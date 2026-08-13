#include "common.h"
#include "spdlog/spdlog.h"
#include "cuda.cuh"

static int antenna_instance_id = 0;
static cudaStream_t stream = nullptr;
static Placements* sim_device_rx_locations = nullptr;
static size_t sim_rx_location_capacity = 0;
using namespace wificuda;

compute_buffers_t gfx, sim;

static int CUDA_INIT(compute_buffers_t& data, const size_t tx_count, const size_t& pxl_count)
{
	if (stream)
	{
		CUDA_CALL(cudaMalloc(&data.__device__phee_minus_alpha_list, tx_count * pxl_count * sizeof(double)), \
			"releasing memory: __device__phee_minus_alpha_list");
		CUDA_CALL(cudaMalloc(&data.__device__gain_RX_grid, tx_count * pxl_count * sizeof(double)), \
			"releasing memory: __device__gain_RX_grid");
		CUDA_CALL(cudaMalloc(&data.__device__pathloss_list, tx_count * pxl_count * sizeof(double)), \
			"releasing memory: __device__pathloss_list");
		CUDA_CALL(cudaMallocHost(&data.___host___hmatrix, tx_count * pxl_count * sizeof(double)), \
			"releasing memory: ___host___hmatrix");
		CUDA_CALL(cudaMalloc(&data.__device__hmatrix, tx_count * pxl_count * sizeof(double)), \
			"releasing memory: __device__hmatrix");
		CUDA_CALL(cudaMalloc(&data.__device__polar_data_theta, tx_count * pxl_count * sizeof(double)), \
			"releasing memory: __device__polar_data_theta");
		CUDA_CALL(cudaMalloc(&data.__device__polar_data_hype, tx_count * pxl_count * sizeof(double)), \
			"releasing memory: __device__polar_data_hype");

		data.cellcount = pxl_count;
		data.num_tx = tx_count;
		return 0;
	}

	CUDA_ERROR("CUDA steam not defined");
	return 1;
}


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
	const double& current_lambda,
	const double& current_spacing,
	const double& theta_c,
	const unsigned& panel_count,
	const float& ant_dim_x,
	const float& ant_dim_y,
	double* d_phee_minus_alpha_list,
	double* d_pathloss_list,
	double* d_gain_RX_grid,
	double* d_polar_data_theta,
	double* d_polar_data_hyp)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < data_size)
	{
		const double& pioverlambda = M_PIl / current_lambda;
		const double& phee_temp = 2 * current_spacing * pioverlambda;
		const double& pl_temp_meters = 4 * pioverlambda;
		const double& antenna_dim_factor = 10 * ant_dim_x * ant_dim_y / (current_lambda * current_lambda);
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

/* for COW and hundreds of handsets only */
static __global__ void numerical_cart2pol_kernel(
	double* polar_theta,
	double* polar_hype,
	const Placements* rx_locations,
	const size_t receiver_count,
	const unsigned tx_x,
	const unsigned tx_y)
{
	const size_t rx_idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (rx_idx >= receiver_count)
	{
		return;
	}

	const double x = static_cast<double>(rx_locations[rx_idx].x)
		- static_cast<double>(tx_x);
	const double y = static_cast<double>(rx_locations[rx_idx].y)
		- static_cast<double>(tx_y);

	polar_theta[rx_idx] = atan2(y, x);
	polar_hype[rx_idx] = hypot(x, y);
}

/*update the antenna array from updated powerand scan angle */
static void global_update(
	compute_buffers_t& data,
	const size_t tx_id,
	const Settings& current)
{
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
	CUDA_CALL(cudaStreamSynchronize(stream), "stream sync");
}

static void global_init(
	compute_buffers_t& data,
	const size_t tx_id,
	const Settings& current
	)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (data.cellcount + threadsPerBlock - 1) / threadsPerBlock;

	antenna_init_kernel <<< blocksPerGrid, threadsPerBlock, 0, stream >>> (
		data.cellcount,
		current.lambda,
		current.spacing,
		current.theta_c,
		current.panel_count,
		current.antenna_dims.x,
		current.antenna_dims.y,
		data.__device__phee_minus_alpha_list,
		data.__device__pathloss_list,
		data.__device__gain_RX_grid,
		data.__device__polar_data_theta,
		data.__device__polar_data_hype
	);

	CUDA_CALL(cudaGetLastError(), "antenna_init_kernel state query");

	data.modified = true;
}

/* for GUI simulation in the whole grid */
void wificuda::graphics_init(unsigned tx_id, const size_t polar_data_size, const Settings& settings)
{
	spdlog::info("Antenna Graphics Re-init");
	if (sim.cellcount == 0 || sim.cellcount != polar_data_size)
	{
		CUDA_FREE(gfx);
		CUDA_INIT(gfx, tx_id, polar_data_size);
	}

	global_init(gfx, tx_id, settings);
}

/* for bare-minimum numerical calculations needed at the mobile_stations only */
void wificuda::numerical_init(unsigned tx_id, const size_t polar_data_size, const Settings& settings)
{
	spdlog::info("Antenna Numerical Re-init");
	if (sim.cellcount == 0 || sim.cellcount != polar_data_size)
	{
		CUDA_FREE(sim);
		CUDA_INIT(sim, tx_id, polar_data_size);
	}

	//simulation.host_hmatrix = host_hmatrix_sim;
	//host_polar_sim = polar_data.data_ptr;

	//cudaStream_t stream;
	//cudaStreamCreate(&stream);

	//cudaMemcpyAsync(d_polar_data_sim, host_polar_sim, polar_data.array_size * sizeof(Polar_Coordinates), cudaMemcpyHostToDevice, stream);

	//cudaStreamSynchronize(stream);
	//cudaStreamDestroy(stream);
	//cudaMemcpy(d_polar_data_sim, polar_data.data_ptr, polar_data.array_size * sizeof(Polar_Coordinates), cudaMemcpyHostToDevice);

	global_init(sim, tx_id, settings);
}

/* update the antenna array from updated power and scan angle */
void wificuda::graphics_update(unsigned tx_id, const Settings& current)
{
	if (gfx.modified && stream)
	{
		spdlog::info("Antenna Graphics update");
		global_update(gfx, tx_id, current);
		gfx.modified = false;
	}
}

void wificuda::numerical_update(unsigned tx_id, const Settings& current)
{
	if (sim.modified && stream)
	{
		spdlog::info("Antenna Numerical update");
		global_update(sim, tx_id, current);
		sim.modified = false;
	}
}

void wificuda::update(unsigned tx_id, const Settings& current)
{
	wificuda::graphics_update(tx_id, current);
	wificuda::numerical_update(tx_id, current);
}

void wificuda::graphics_recalc_polar(
	unsigned tx_id,
	const Placements& location,
	const Dimensions<unsigned>& grid_size,
	const Settings& settings)
{
	const size_t cellcount = grid_size.x * grid_size.y;

	if (gfx.cellcount != cellcount) // render area changed
	{
		CUDA_FREE(gfx);
		CUDA_INIT(gfx, gfx.num_tx, cellcount);
	}

	//ensure_compute_layout(gfx, tx_location.size(), cellcount);
	gfx.grid_width = grid_size.x;
	gfx.grid_height = grid_size.y;

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

	gfx.modified = true;
}


/* Recalculate transmitter-to-receiver polar data for the simulation. */
void wificuda::numerical_recalc_polar(
	unsigned tx_id,
	const Placements& tx_location,
	const placement_v& rx_locations)
{
	if (!stream || rx_locations.empty())
	{
		CUDA_ERROR("numerical_recalc_polar requires a stream, a transmitter, and receivers");
	}

	//ensure_compute_layout(sim, tx_locations.size(), rx_locations.size());

	if (sim_rx_location_capacity != rx_locations.size())
	{
		if (sim_device_rx_locations)
		{
			CUDA_CALL(cudaFree(sim_device_rx_locations), "%s", "releasing simulation receiver locations");
		}

		CUDA_CALL(cudaMalloc(
			&sim_device_rx_locations,
			rx_locations.size() * sizeof(Placements)),
			"%s", "allocating simulation receiver locations");
		sim_rx_location_capacity = rx_locations.size();
	}

	CUDA_CALL(cudaMemcpyAsync(
		sim_device_rx_locations,
		rx_locations.data(),
		rx_locations.size() * sizeof(Placements),
		cudaMemcpyHostToDevice,
		stream),
		"%s", "copying simulation receiver locations");

	int threads_per_block = 256;
	int blocks_per_grid = \
		(rx_locations.size() + threads_per_block - 1) / threads_per_block;


	const size_t offset = tx_id * rx_locations.size();
	numerical_cart2pol_kernel <<<blocks_per_grid, threads_per_block, 0, stream>>> (
		sim.__device__polar_data_theta + offset,
		sim.__device__polar_data_hype + offset,
		sim_device_rx_locations,
		rx_locations.size(),
		tx_location.x,
		tx_location.y);
	CUDA_CALL(cudaGetLastError(), "%s", "numerical_cart2pol_kernel launch");


	sim.modified = true;
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
double wificuda::coeff(compute_buffers_t& data, unsigned tx_id, const size_t data_idx)
{
	if (tx_id >= data.num_tx || data_idx >= data.cellcount)
	{
		CUDA_ERROR("coefficients are out of range: "\
			"stations(%u):station_num(%u), cellcount(%u):pixel_idx(%u)",
			data.num_tx, tx_id, data.cellcount, data_idx);
	}
	const size_t index = tx_id * data.cellcount + data_idx;
	return data.___host___hmatrix[index];
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
		if (sim_device_rx_locations)
		{
			CUDA_CALL(cudaFree(sim_device_rx_locations), "%s", "releasing simulation receiver locations");
			sim_device_rx_locations = nullptr;
			sim_rx_location_capacity = 0;
		}

		cudaError_t err = cudaStreamDestroy(stream);
		if (err != cudaSuccess)
		{
			fprintf(stderr, "gpu_teardown() destroy stream failed: %s\n", cudaGetErrorString(err));
			exit(-1);
		}
		stream = nullptr;
	}
}

void wificuda::gpu_init()
{
	cudaError_t err = cudaStreamCreate(&stream);
	if (err != cudaSuccess)
	{
		fprintf(stderr, "gpu_init() launch failed: %s\n", cudaGetErrorString(err));
		exit(-1);
	}

	//spdlog::info("CUDA stream is now succefully initialized");
}

network_package::AAntenna::AAntenna(
	const unsigned& init_panel_count,
	const double& init_lambda,
	const double& init_antenna_spacing,
	const double& init_antenna_orientation_rads,
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
		init_antdims
	},
dummy(nullptr) // constant initial setup
{

}