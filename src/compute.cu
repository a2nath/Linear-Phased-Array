#include "common.h"
#include "spdlog/spdlog.h"
#include "cuda.cuh"

static int antenna_instance_id = 0;
static cudaStream_t stream = nullptr;

bool network_package::gfx_update = true;
bool network_package::sim_update = true;

static compute_buffers_t buffers[2] = {};
static compute_buffers_t& gfx = buffers[GFX];
static compute_buffers_t& sim = buffers[SIM];

struct sim_compute_buffers_t // extension of the sim buffers
{
	Placements* __device__rx_locations = nullptr; // for numerical numbers
	size_t count = 0;
	// for graphical numbers - you extrapolate the locations using grid size.
} simx;

struct gfx_compute_buffers_t // extension of the gfx buffers
{
	Dimensions<unsigned> grid_dims = { 0, 0 };
} gfxx;

using namespace wificuda;

static int CUDA_INIT(compute_buffers_t& data, const size_t tx_count, const size_t& pxl_count)
{
	if (stream)
	{
		if (tx_count == 0 || pxl_count == 0)
		{
			CUDA_ERROR("cannot initialize CUDA buffers with %zu transmitters and %zu cells", tx_count, pxl_count);
		}

		data.cellcount = pxl_count;
		data.num_tx = tx_count;
		size_t datasize = data.num_tx * data.cellcount * sizeof(double);

		spdlog::info("CUDA INIT: " + str(data.num_tx) \
			+ " transmitters, " + str(data.cellcount) \
			+ " datapoints, " + str(datasize) + " total bytes");

		CUDA_CALL(cudaMalloc(&data.__device__phee_minus_alpha_list, datasize), \
			"allocating memory: __device__phee_minus_alpha_list");
		CUDA_CALL(cudaMalloc(&data.__device__gain_RX_grid, datasize), \
			"allocating memory: __device__gain_RX_grid");
		CUDA_CALL(cudaMalloc(&data.__device__pathloss_list, datasize), \
			"allocating memory: __device__pathloss_list");
		CUDA_CALL(cudaMallocHost(&data.___host___hmatrix, datasize), \
			"allocating memory: ___host___hmatrix");
		CUDA_CALL(cudaMalloc(&data.__device__hmatrix, datasize),
			"allocating memory: __device__hmatrix");
		CUDA_CALL_PRINT(cudaMalloc(&data.__device__polar_data_theta, datasize),
			"allocating memory: __device__polar_data_theta");
		CUDA_CALL_PRINT(cudaMalloc(&data.__device__polar_data_hype, datasize),
			"allocating memory: __device__polar_data_hype");

		spdlog::info("CUDA INIT: allocated " + str(7 * datasize) + " %zu total bytes");

		if (data.coeff_changed_arr.empty())
			data.coeff_changed_arr.resize(tx_count, true);

		return 0;
	}

	CUDA_ERROR("CUDA steam not defined");
	return 1;
}

static int CUDA_INIT(compute_buffers_t& data, const size_t tx_count, const Dimensions<unsigned>& dim)
{
	gfx.type = GFX;
	gfxx.grid_dims = dim;

	return CUDA_INIT(data, tx_count, dim.count());
}

/* Flush the data from the previous allocation. This function does not get called during init, as = {} is enough */
static int CUDA_FREE(compute_buffers_t& data)
{
	CUDA_CALL(cudaFree(data.__device__phee_minus_alpha_list), \
		"releasing memory: __device__phee_minus_alpha_list");
	CUDA_CALL(cudaFree(data.__device__gain_RX_grid), \
		"releasing memory: __device__gain_RX_grid");
	CUDA_CALL(cudaFree(data.__device__pathloss_list), \
		"releasing memory: __device__pathloss_list");
	CUDA_CALL(cudaFreeHost(data.___host___hmatrix), \
		"releasing memory: ___host___hmatrix");
	CUDA_CALL(cudaFree(data.__device__hmatrix), \
		"releasing memory: __device__hmatrix");
	CUDA_CALL(cudaFree(data.__device__polar_data_theta), \
		"releasing memory: __device__polar_data_theta");
	CUDA_CALL(cudaFree(data.__device__polar_data_hype), \
		"releasing memory: __device__polar_data_hype");

	if (data.type == SIM && simx.__device__rx_locations) {
		CUDA_CALL(cudaFree(simx.__device__rx_locations), \
			"releasing memory: __device__rx_locations");
		simx = {};
	}

	if (data.type == GFX) {
		gfxx = {};
	}

	data = {};
	return 0;
}

/* update the antenna array from updated power and scan angle */
static __global__ void antenna_update_kernel(
	const size_t size,
	const double scan_angle,
	const unsigned panel_count,
	double* d_phee_minus_alpha_list,
	double* d_gain_RX_grid,
	double* d_pathloss_list,
	double* d_hmatrix)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

	/* update the antenna gain Gtx */
	if (idx < size)
	{
		double phee = (d_phee_minus_alpha_list[idx] + scan_angle) / 2;

		double sin_term = panel_count * sinf(phee);
		double gain_factor_antenna_system = d_gain_RX_grid[idx]; // xN antennas already

		if (sin_term != 0)
		{
			double pow_base = sinf(panel_count * phee) / sin_term;
			gain_factor_antenna_system *= pow_base * pow_base;
		}

		/* update the channel matrix */
		d_hmatrix[idx] = gain_factor_antenna_system / d_pathloss_list[idx];
	}
}

/* re-calc the signal outs to handsets only (before calling update!) */
static __global__ void antenna_init_kernel(
	const size_t data_size,
	const double current_lambda,
	const double current_spacing,
	const double theta_c,
	const unsigned panel_count,
	const float ant_dim_x,
	const float ant_dim_y,
	double* d_phee_minus_alpha_list,
	double* d_pathloss_list,
	double* d_gain_RX_grid,
	double* d_polar_data_theta,
	double* d_polar_data_hyp)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < data_size)
	{
		const double pioverlambda = M_PIl / current_lambda;
		const double phee_temp = 2 * current_spacing * pioverlambda;
		const double pl_temp_meters = 4 * pioverlambda;
		const double antenna_dim_factor = 10 * ant_dim_x * ant_dim_y / (current_lambda * current_lambda);
		double m_factor = ant_dim_x * pioverlambda;


		double theta_minus_thetaC = d_polar_data_theta[idx] - theta_c;
		double m = m_factor * sinf(theta_minus_thetaC);
		double pow_base = (1 + cos(theta_minus_thetaC)) / 2;
		double singleant_gain = antenna_dim_factor * pow_base * pow_base; // pow(pow_base, 2) equivalent on HOST

		if (m != 0)
		{
			singleant_gain *= pow(sin(m) / m, 2);
		}

		d_phee_minus_alpha_list[idx] = phee_temp * sin(theta_minus_thetaC);
		d_pathloss_list[idx] = pow(pl_temp_meters * d_polar_data_hyp[idx], 2);
		d_gain_RX_grid[idx] = singleant_gain * panel_count;
	}
}

/* for COW and dots for heatmap */
static __global__ void graphics_cart2pol_kernel(
	double* polar_theta,
	double* polar_hype,
	const size_t cellcount,
	const unsigned grid_width,
	const unsigned tx_x,
	const unsigned tx_y)
{
	const size_t pixel_idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (pixel_idx >= cellcount)
	{
		return;
	}

	const unsigned col = static_cast<unsigned>(pixel_idx % grid_width);
	const unsigned row = static_cast<unsigned>(pixel_idx / grid_width);
	const double x = static_cast<double>(col) - static_cast<double>(tx_x);
	const double y = static_cast<double>(row) - static_cast<double>(tx_y);

	polar_theta[pixel_idx] = atan2(y, x);
	polar_hype[pixel_idx] = hypot(x, y);
}

/* for COW and hundreds of handsets only
*
* A __global__ function executes on the GPU and cannot directly access ordinary host globals.
* Also, the name __device__rx_locations is only a name; it does not apply CUDA's __device__ qualifier.
*/

static __global__ void numerical_cart2pol_kernel(
	double* polar_theta,
	double* polar_hype,
	const Placements* rx_locations,
	const size_t rx_station_num,
	const unsigned tx_x,
	const unsigned tx_y)
{
	const size_t rx_idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (rx_idx >= rx_station_num)
	{
		return;
	}

	const double x = (double)rx_locations[rx_idx].x - (double)tx_x;
	const double y = (double)rx_locations[rx_idx].y - (double)tx_y;

	polar_theta[rx_idx] = atan2(y, x);
	polar_hype[rx_idx] = hypot(x, y);
}

/*update the antenna array from updated powerand scan angle */
static void global_update(
	compute_buffers_t& data,
	const size_t tx_id,
	const Settings& current)
{
	if (data.cellcount == 0 || tx_id >= data.num_tx)
	{
		CUDA_ERROR("invalid antenna update layout: transmitter %zu of %zu, cell count %zu",
			tx_id, data.num_tx, data.cellcount);
	}

	int threadsPerBlock = 256;
	int blocksPerGrid = (data.cellcount + threadsPerBlock - 1) / threadsPerBlock;

	const size_t offset = tx_id * data.cellcount;
	const size_t bytes = data.cellcount * sizeof(double);

	antenna_update_kernel <<< blocksPerGrid, threadsPerBlock, 0, stream >>> (
		data.cellcount,
		current.alpha,
		current.panel_count,
		data.__device__phee_minus_alpha_list + offset,
		data.__device__gain_RX_grid + offset,
		data.__device__pathloss_list + offset,
		data.__device__hmatrix + offset
	);

	CUDA_CALL(cudaGetLastError(), "antenna_update_kernel state query");

	CUDA_DEBUG_SYNC(stream);

}

static void global_init(
	compute_buffers_t& data,
	const size_t tx_id,
	const Settings& current
	)
{
	if (data.cellcount == 0)
	{
		CUDA_ERROR("GPU memory not reserved on %s. Call the handset function",
			data.type == SIM ? "SIM" : "GFX");
	}

	int threadsPerBlock = 256;
	int blocksPerGrid = (data.cellcount + threadsPerBlock - 1) / threadsPerBlock;
	const size_t offset = tx_id * data.cellcount;

	antenna_init_kernel <<< blocksPerGrid, threadsPerBlock, 0, stream >>> (
		data.cellcount,
		current.lambda,
		current.spacing,
		current.theta_c,
		current.panel_count,
		current.antenna_dims.x,
		current.antenna_dims.y,
		data.__device__phee_minus_alpha_list + offset,
		data.__device__pathloss_list + offset,
		data.__device__gain_RX_grid + offset,
		data.__device__polar_data_theta + offset,
		data.__device__polar_data_hype + offset
	);

	CUDA_CALL(cudaGetLastError(), "antenna_init_kernel() failed");

	CUDA_DEBUG_SYNC(stream);
}

/* for GUI simulation in the whole grid */
void wificuda::graphical_init(unsigned tx_id, const Dimensions<unsigned>& dims, AAntennaTable& anttable)
{
	if (gfx.cellcount == 0) {
		spdlog::info("Antenna Graphics Init");
		gui_resize(anttable.size(), dims);
	}

	if (gfxx.grid_dims != dims) {
		spdlog::info("Antenna Graphics Re-init");
		gui_resize(anttable.size(), dims);
	}

	global_init(gfx, tx_id, anttable[tx_id].settings());
	network_package::gfx_update = true;
}

/* Simulation outputs signal-to-noise-ratio for each handsets
 * Allocated memory rx_count * tx_count
 * Input: station_id
 *        number of handset
 *        antenna_table
 * Output: none
 * Function will calculate the initial parameters that will be
 * used during the antenna update (power and scan angle)
 */
void wificuda::numerical_init(unsigned tx_id, unsigned ms_num, AAntennaTable& anttable)
{
	global_init(sim, tx_id, anttable[tx_id].settings());
	network_package::sim_update = true;
}

/* update the antenna array from updated power and scan angle */
void wificuda::graphical_update(unsigned tx_id, const Settings& current)
{
	if (network_package::gfx_update && stream)
	{
		spdlog::info("Antenna Graphics update");
		global_update(gfx, tx_id, current);
		network_package::gfx_update = false;
		gfx.coeff_changed_arr[tx_id] = true;
	}
}

void wificuda::numerical_update(unsigned tx_id, const Settings& current)
{
	if (network_package::sim_update && stream)
	{
		spdlog::info("Antenna Numerical update");
		global_update(sim, tx_id, current);
		network_package::sim_update = false;
		sim.coeff_changed_arr[tx_id] = true;
	}
}

void wificuda::update(unsigned tx_id, const Settings& current)
{
	wificuda::graphical_update(tx_id, current);
	wificuda::numerical_update(tx_id, current);
}


void wificuda::graphical_recalc_polar(
	unsigned tx_id,
	const Placements& location,
	const Dimensions<unsigned>& grid_size,
	const Settings& settings)
{
	const size_t cellcount = grid_size.x * grid_size.y;

	if (gfxx.grid_dims != grid_size) // render area changed
	{
		const size_t tx_count = gfx.num_tx;
		CUDA_FREE(gfx);
		CUDA_INIT(gfx, tx_count, grid_size);
	}

	int threadsPerBlock = 256;
	int blocksPerGrid = (cellcount + threadsPerBlock - 1) / threadsPerBlock;

	const size_t offset = tx_id * cellcount;
	graphics_cart2pol_kernel <<< blocksPerGrid, threadsPerBlock, 0, stream >>> (
		gfx.__device__polar_data_theta + offset,
		gfx.__device__polar_data_hype + offset,
		cellcount,
		grid_size.x,
		location.x,
		location.y);
	CUDA_CALL(cudaGetLastError(), "%s", "graphics_cart2pol_kernel launch");

	CUDA_DEBUG_SYNC(stream);

	network_package::gfx_update = true;
}


/* Recalculate transmitter-to-receiver polar data for the simulation. */
void wificuda::numerical_recalc_polar(
	unsigned tx_id,
	const size_t tx_count,
	const Placements& tx_new_location)
{
	if (!stream)
	{
		CUDA_ERROR("numerical_recalc_polar requires a stream, 1+ transmitter, and 1+ receivers");
	}

	if (simx.count == 0)
	{   /* part of init */
		CUDA_ERROR("rx location data is not filled out in the device memory. Call the handset function");
	}

	if (sim.cellcount == 0)
	{   /* part of init */
		spdlog::info("Antenna GPU arrays allocating. CUDA init on sim data");
		sim.type = SIM;
		CUDA_INIT(sim, tx_count, simx.count);
	}
	//ensure_compute_layout(sim, tx_locations.size(), rx_locations.size());
	//
	//if (__device__rx_locations != rx_locations.size())
	//{
	//	if (__device__rx_locations)
	//	{
	//		CUDA_CALL(cudaFree(__device__rx_locations), "%s", "releasing simulation receiver locations");
	//	}
	//
	//	CUDA_CALL(cudaMalloc(
	//		&__device__rx_locations,
	//		rx_locations.size() * sizeof(Placements)),
	//		"%s", "allocating simulation receiver locations");
	//	sim_rx_location_capacity = rx_locations.size();
	//}

	const size_t offset = tx_id * simx.count;

	int threads_per_block = 256;
	int blocks_per_grid = \
		(simx.count + threads_per_block - 1) / threads_per_block;

	numerical_cart2pol_kernel <<< blocks_per_grid, threads_per_block, 0, stream >>> (
		sim.__device__polar_data_theta + offset,
		sim.__device__polar_data_hype + offset,
		simx.__device__rx_locations,
		simx.count,
		tx_new_location.x,
		tx_new_location.y);

	CUDA_CALL_PRINT(cudaGetLastError(), "numerical_cart2pol_kernel launch with \nTX-id:%d"\
		"with %d locations in device memory\nSim array size: %d\nNew location: (%d,%d)\nOffset:%d\n",
		tx_id, simx.count, sim.cellcount, tx_new_location.x, tx_new_location.y, offset);

	CUDA_DEBUG_SYNC(stream);

	network_package::sim_update = true;
}

void wificuda::gui_resize(unsigned tx_count, const Dimensions<unsigned>& dims)
{
	if (gfx.num_tx == tx_count && gfxx.grid_dims == dims)
	{
		return;
	}

	spdlog::info("GUI Resize");
	CUDA_FREE(gfx);
	CUDA_INIT(gfx, tx_count, dims);
}

/* Recalculate the GUI grid for every transmitter. */
//void old_recalc_polar(
//	const Dimensions<unsigned>& grid_size,
//	const Placements& tx_location)
//{
//	if (!stream)
//	{
//		CUDA_ERROR("%s", "old_recalc_polar requires a stream, grid, and transmitters");
//	}
//
//	const size_t cellcount = grid_size.x * grid_size.y;
//
//	sim.grid_width = grid_size.x;
//	sim.grid_height = grid_size.y;
//
//	int threadsPerBlock = 256;
//	int blocksPerGrid = (cellcount + threadsPerBlock - 1) / threadsPerBlock;
//
//	const size_t offset = tx_id * cellcount;
//	graphics_cart2pol_kernel << < blocksPerGrid, threadsPerBlock, 0, stream >> > (
//		sim.__device__polar_data_theta + offset,
//		sim.__device__polar_data_hype + offset,
//		cellcount,
//		grid_size.x,
//		tx_location[tx_id].x,
//		tx_location[tx_id].y);
//	CUDA_CALL(cudaGetLastError(), "%s", "old_recalc_polar launch");
//}
/* hatrix with respect to pixel index (flattened from 2D) */
static inline double* coeff_arr_ptr(compute_buffers_t& data, unsigned tx_id)
{
	if (tx_id >= data.num_tx)
	{
		CUDA_ERROR("coefficients are out of range: "\
			"stations(%u):station_num(%u), cellcount(%u)",
			data.num_tx, tx_id, data.cellcount);
	}

	const size_t offset = tx_id * data.cellcount;
	if (data.coeff_changed_arr[tx_id])
	{
		const size_t bytes = data.cellcount * sizeof(double);

		/* device to host must happen in host side */
		CUDA_CALL(cudaMemcpyAsync(
			data.___host___hmatrix + offset,
			data.__device__hmatrix + offset,
			bytes,
			cudaMemcpyDeviceToHost,
			stream),
			"copy hmatrix memory from GPU to CPU"
		);

		/* required to finish the sync above before CPU reads result */
		CUDA_CALL(cudaStreamSynchronize(stream), "stream sync hmatrix from device -> host");

		data.coeff_changed_arr[tx_id] = false;
	}

	return data.___host___hmatrix + offset;
}


/* sync GUI data from device -> host and then return the pointer */
double* wificuda::gcoeff_ptr(unsigned tx_id)
{
	return coeff_arr_ptr(gfx, tx_id);
}

/* sync SIM data from device -> host and then return the pointer */
double* wificuda::scoeff_ptr(unsigned tx_id)
{
	return coeff_arr_ptr(sim, tx_id);
}

void wificuda::gpu_teardown()
{
	if (stream)
	{
		cudaError_t synerr = cudaStreamSynchronize(stream);
		if (synerr != cudaSuccess)
		{
			fprintf(stderr, "gpu_teardown() synchronize failed: %s\n", cudaGetErrorString(synerr));
			exit(-1);
		}

		CUDA_FREE(gfx);
		CUDA_FREE(sim);

		cudaError_t err = cudaStreamDestroy(stream);
		if (err != cudaSuccess)
		{
			fprintf(stderr, "gpu_teardown() destroy stream failed: %s\n", cudaGetErrorString(err));
			exit(-1);
		}
		stream = nullptr;
	}
}

/* Init function to set the cuda memory for handset placements
   When used- numerical_recalc_polar function will use the positions of the handsets
			  in relation to the TXs.
   numerical_recalc_polar is calculated in __device__
*/
void wificuda::handset_placement_sync(const placement_v& handset_locations)
{
	if (stream)
	{
		if (!simx.__device__rx_locations)
		{
			CUDA_CALL(cudaMalloc(&simx.__device__rx_locations, handset_locations.size() * sizeof(Placements)), \
				"allocating memory: __device__rx_locations");
		}

		CUDA_CALL(cudaMemcpyAsync(
			simx.__device__rx_locations,
			handset_locations.data(),
			handset_locations.size() * sizeof(Placements),
			cudaMemcpyHostToDevice,
			stream),
			"%s", "copying simulation receiver locations");

		simx.count = handset_locations.size(); // set the count if success
		CUDA_DEBUG_SYNC(stream);
		return;
	}

	CUDA_WARN("Stream not initialized when allocating memory for handsets simulation wide");
}


void wificuda::gpu_init()
{
	if (!stream)
	{
		cudaError_t err = cudaStreamCreate(&stream);
		if (err != cudaSuccess)
		{
			fprintf(stderr, "gpu_init() launch failed: %s\n", cudaGetErrorString(err));
			exit(-1);
		}
	}
	//spdlog::info("CUDA stream is now succefully initialized");
}

network_package::AAntenna::AAntenna(
	const unsigned& init_panel_count,
	const double& init_lambda,
	const double& init_antenna_spacing,
	const double& init_antenna_orientation_rads,
	const Placements& init_location,
	const antennadim& init_antdims)
	:
	instance_id(antenna_instance_id++),
	initial{
		0,
		std::numeric_limits<double>::min(),
		init_panel_count,
		init_lambda,
		init_antenna_spacing,
		init_antenna_orientation_rads,
		init_location,
		init_antdims
	},
	current(initial),
	full_re_init(true)
{
}
