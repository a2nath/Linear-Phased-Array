#if defined(__CUDACC__) || !defined(__device__)
#include <math_constants.h>
#include "cuda_runtime.h"
#include "network.h"

using namespace network_package;

#ifdef __INTELLISENSE__
#define __CUDA_ARCH__ 800
#endif
#define CUDA_CALL(call)                                                         \
	{                                                                           \
		cudaError_t err = call;                                                 \
		if (err != cudaSuccess) {                                               \
			fprintf(stderr, "CUDA error in file '%s' in line %i: %s.\n",        \
					__FILE__, __LINE__, cudaGetErrorString(err));               \
			exit(EXIT_FAILURE);                                                 \
		}                                                                       \
	}

__constant__ double CONSTANT_PI = 3.1415926535897931;

static double* d_phee_minus_alpha_list_gfx = nullptr;
static double* d_gain_RX_grid_gfx = nullptr;
static double* d_pathloss_list_gfx = nullptr;

static double* host_hmatrix_gfx = nullptr;
static double* d_hmatrix_gfx = nullptr;

static Polar_Coordinates* d_polar_data_gfx = nullptr;
static Polar_Coordinates* host_polar_gfx = nullptr;

static size_t d_malloc_size_gfx = 0;

static double* d_phee_minus_alpha_list_sim = nullptr;
static double* d_gain_RX_grid_sim = nullptr;
static double* d_pathloss_list_sim = nullptr;

static double* host_hmatrix_sim = nullptr;
static double* d_hmatrix_sim = nullptr;

static Polar_Coordinates* d_polar_data_sim = nullptr;
static Polar_Coordinates* host_polar_sim = nullptr;

static size_t d_malloc_size_sim = 0;

__host__ void CUDA_GINIT(const size_t& malloc_size)
{
	cudaMalloc(&d_phee_minus_alpha_list_gfx, malloc_size * sizeof(double));
	cudaMalloc(&d_gain_RX_grid_gfx, malloc_size * sizeof(double));
	cudaMalloc(&d_pathloss_list_gfx, malloc_size * sizeof(double));
	cudaMallocHost(&host_hmatrix_gfx, malloc_size * sizeof(double));
	cudaMalloc(&d_hmatrix_gfx, malloc_size * sizeof(double));
	cudaMalloc(&d_polar_data_gfx, malloc_size * sizeof(Polar_Coordinates));
	d_malloc_size_gfx = malloc_size;
}

__host__ void CUDA_INIT(const size_t& malloc_size)
{
	cudaMalloc(&d_phee_minus_alpha_list_sim, malloc_size * sizeof(double));
	cudaMalloc(&d_gain_RX_grid_sim, malloc_size * sizeof(double));
	cudaMalloc(&d_pathloss_list_sim, malloc_size * sizeof(double));
	cudaMallocHost(&host_hmatrix_sim, malloc_size * sizeof(double));
	cudaMalloc(&d_hmatrix_sim, malloc_size * sizeof(Polar_Coordinates));
	cudaMalloc(&d_polar_data_sim, malloc_size * sizeof(Polar_Coordinates));
	d_malloc_size_sim = malloc_size;
}

__host__ void CUDA_GFREE()
{
	cudaFree(d_phee_minus_alpha_list_gfx);
	cudaFree(d_gain_RX_grid_gfx);
	cudaFree(d_pathloss_list_gfx);
	cudaFreeHost(host_hmatrix_gfx);
	cudaFree(d_hmatrix_gfx);
	cudaFree(d_polar_data_gfx);
}

__host__ void CUDA_FREE()
{
	cudaFree(d_phee_minus_alpha_list_sim);
	cudaFree(d_gain_RX_grid_sim);
	cudaFree(d_pathloss_list_sim);
	cudaFreeHost(host_hmatrix_sim);
	cudaFree(d_hmatrix_sim);
	cudaFree(d_polar_data_sim);
}

/* update the antenna array from updated power and scan angle */
__global__ void antenna_update_kernel(
	const size_t size,
	const double& scan_angle,
	const unsigned& panel_count,
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

/*update the antenna array from updated powerand scan angle */
__host__ void AAntenna::update(
	const size_t& malloc_size,
	double* phee_minus_alpha_list,
	double* gain_RX_grid,
	double* pathloss_list,
	double* gpu_hmatrix,
	double* host_hmatrix)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (malloc_size + threadsPerBlock - 1) / threadsPerBlock;

	cudaStream_t stream;
	cudaStreamCreate(&stream);

	antenna_update_kernel <<< blocksPerGrid, threadsPerBlock >>> (
		malloc_size, current.alpha, current.panel_count, phee_minus_alpha_list, gain_RX_grid, pathloss_list, gpu_hmatrix);

	cudaMemcpyAsync(host_hmatrix, gpu_hmatrix, malloc_size * sizeof(double), cudaMemcpyDeviceToHost, stream);
	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);

	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "antenna_update_kernel() launch failed: %s\n", cudaGetErrorString(cudaStatus));
		exit(-1);
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		exit(-1);
	}
}



/* re-calc the signal outs to handsets only (before calling update!) */
__global__ void antenna_init_kernel(
	const size_t data_size,
	const double& current_lambda,
	const double& current_spacing,
	const double& theta_c,
	const unsigned& panel_count,
	const antennadim& antenna_dims,
	double* d_phee_minus_alpha_list,
	double* d_pathloss_list,
	double* d_gain_RX_grid,
	Polar_Coordinates* d_polar_data)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < data_size)
	{
		const double& pioverlambda = CONSTANT_PI / current_lambda;
		const double& phee_temp = 2 * current_spacing * pioverlambda;
		const double& pl_temp_meters = 4 * pioverlambda;
		const double& antenna_dim_factor = 10 * antenna_dims.x * antenna_dims.y / (current_lambda * current_lambda);
		double m_factor = antenna_dims.x * pioverlambda;

		auto& cell_polar_data = d_polar_data[idx];
		double theta_minus_thetaC = cell_polar_data.theta - theta_c;
		double m = m_factor * sinf(theta_minus_thetaC);
		double pow_base = (1 + cached::cos(theta_minus_thetaC)) / 2;
		double singleant_gain = antenna_dim_factor * pow_base * pow_base; // pow(pow_base, 2) equivalent on HOST

		if (m != 0)
		{
			singleant_gain *= pow(cached::sin(m) / m, 2);
		}

		d_phee_minus_alpha_list[idx] = phee_temp * cached::sin(theta_minus_thetaC);
		d_pathloss_list[idx] = pow(pl_temp_meters * cell_polar_data.hype, 2);
		d_gain_RX_grid[idx] = singleant_gain * panel_count;
	}
}

__host__ void AAntenna::init(
	const size_t& d_malloc_size,
	double* d_phee_minus_alpha_list,
	double* d_pathloss_list,
	double* d_gain_RX_grid,
	const Polar_Coordinates* d_polar_data)
{
	int threadsPerBlock = 256;
	int blocksPerGrid = (d_malloc_size + threadsPerBlock - 1) / threadsPerBlock;

	antenna_init_kernel <<< blocksPerGrid, threadsPerBlock >>> (
		d_malloc_size, current.lambda, current.spacing, current.theta_c, current.panel_count, current.antenna_dims,
		d_phee_minus_alpha_list, d_pathloss_list, d_gain_RX_grid, d_polar_data);

	cudaError_t cudaStatus = cudaGetLastError();

	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "antenna_init_kernel() launch failed: %s\n", cudaGetErrorString(cudaStatus));
		exit(-1);
	}

	/* cudaDeviceSynchronize waits for the kernel to finish, and returns
	 any errors encountered during the launch. */
	cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		exit(-1);
	}
}

/* update the antenna array from updated power and scan angle */
__host__ void AAntenna::graphics_update()
{
	if (graphic.modified)
	{
		update(d_malloc_size_gfx, d_phee_minus_alpha_list_gfx, d_gain_RX_grid_gfx, d_pathloss_list_gfx, d_hmatrix_gfx, host_hmatrix_gfx);
		graphic.modified = false;
	}
}

__host__ void AAntenna::numerical_update()
{
	if (simulation.modified)
	{
		update(d_malloc_size_sim, d_phee_minus_alpha_list_sim, d_gain_RX_grid_sim, d_pathloss_list_sim, d_hmatrix_sim, host_hmatrix_sim);
		simulation.modified = false;
	}
}


/* for GUI simulation in the whole grid */
__host__ void AAntenna::graphics_init(PolarArray& polar_data)
{
	if (d_malloc_size_sim == 0 || d_malloc_size_sim != polar_data.array_size)
	{
		CUDA_GFREE();
		CUDA_GINIT(polar_data.array_size);
	}

	graphic.host_hmatrix = host_hmatrix_gfx;
	host_polar_gfx = polar_data.data_ptr;

	cudaStream_t stream;
	cudaStreamCreate(&stream);

	cudaMemcpyAsync(d_polar_data_gfx, host_polar_gfx, polar_data.array_size * sizeof(Polar_Coordinates), cudaMemcpyHostToDevice, stream);

	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);
	//cudaMemcpy(d_polar_data_gfx, host_polar_gfx, polar_data.array_size * sizeof(Polar_Coordinates), cudaMemcpyHostToDevice);

	init(d_malloc_size_gfx, d_phee_minus_alpha_list_gfx, d_pathloss_list_gfx, d_gain_RX_grid_gfx, d_polar_data_gfx);
	graphic.modified = true;
}


/* for bare-minimum numerical calculations needed at the mobile_stations only */
__host__ void AAntenna::numerical_init(PolarArray& polar_data)
{
	if (d_malloc_size_sim == 0 || d_malloc_size_sim != polar_data.array_size)
	{
		CUDA_FREE();
		CUDA_INIT(polar_data.array_size);
	}

	simulation.host_hmatrix = host_hmatrix_sim;
	host_polar_sim = polar_data.data_ptr;

	cudaStream_t stream;
	cudaStreamCreate(&stream);

	cudaMemcpyAsync(d_polar_data_sim, host_polar_sim, polar_data.array_size * sizeof(Polar_Coordinates), cudaMemcpyHostToDevice, stream);

	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);
	//cudaMemcpy(d_polar_data_sim, polar_data.data_ptr, polar_data.array_size * sizeof(Polar_Coordinates), cudaMemcpyHostToDevice);

	init(d_malloc_size_sim, d_phee_minus_alpha_list_sim, d_pathloss_list_sim, d_gain_RX_grid_sim, d_polar_data_sim);
	simulation.modified = true;
}


__host__ AAntenna::~AAntenna()
{
	cudaError_t cudaStatus = cudaDeviceSynchronize();

	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d in Destructor of Antenna object!\n", cudaStatus);
		exit(-1);
	}

	CUDA_GFREE();
	CUDA_FREE();
}

/* hatrix with respect to pixel index (flattened from 2D) */
const double& AAntenna::gcoeff(const unsigned& pixel_idx) const
{
	return graphic.host_hmatrix[pixel_idx];
}

const double& AAntenna::coeff(const unsigned& rx_sta) const
{
	return simulation.host_hmatrix[rx_sta];
}

#endif
