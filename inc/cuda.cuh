#pragma once

#if defined(__CUDACC__) || !defined(__device__)
#include <math_constants.h>
#include "cuda_runtime.h"
#include "network.cuh"

//#include "spdlog/spdlog.h"

#define FMT_USE_INLINE_VARIABLES 0
#ifdef __INTELLISENSE__
#define __CUDA_ARCH__ 800
#endif

#include <cstdio>
#include <stdexcept>

#define CUDA_ERROR(message, ...) do { \
	char cuda_message[1024]; \
	std::snprintf(cuda_message, sizeof(cuda_message), \
		message, __VA_ARGS__); \
	fprintf(stderr, "ERROR [%s:%d]: %s\n", __FILE__, __LINE__, cuda_message); \
	exit(-1); \
} while(0)

#define CUDA_WARN(message, ...) do { \
	char cuda_message[1024]; \
	std::snprintf(cuda_message, sizeof(cuda_message), \
		message, __VA_ARGS__); \
	printf("WARNING [%s:%d]: %s\n", __FILE__, __LINE__, cuda_message); \
} while(0)

#define CUDA_CALL(x, message, ...) do { \
	const cudaError_t cuda_call_status = (x); \
	if(cuda_call_status != cudaSuccess) { \
		fprintf(stderr, "%s\n", cudaGetErrorString(cuda_call_status)); \
		CUDA_ERROR(message, __VA_ARGS__); \
	} \
} while(0)

/* synchronize after critical functions if debug is defined */
#if defined(_DEBUG) || !defined(NDEBUG)
#define CUDA_DEBUG_SYNC(s) do { \
	CUDA_CALL(cudaStreamSynchronize(s), "debug sync"); \
} while(0)

#define CUDA_CALL_PRINT(x, message, ...) do { \
	const cudaError_t cuda_call_status = (x); \
	printf(message, __VA_ARGS__); \
	printf("\n"); \
	if(cuda_call_status != cudaSuccess) { \
		fprintf(stderr, "%s\n", cudaGetErrorString(cuda_call_status)); \
		CUDA_ERROR("L>cuda error call print"); \
	} \
} while(0)
#else
#define CUDA_DEBUG_SYNC(s) do { } while (0)
#define CUDA_CALL_PRINT CUDA_CALL
#endif


/* CUDA datasize = datasize x num_tx*/
static enum BUFFER_TYPE { GFX = 0, SIM, NONE };

struct compute_buffers_t
{
	BUFFER_TYPE type = NONE;
	size_t cellcount = 0; // pixel or recievers
	size_t num_tx = 0;

	double* __device__phee_minus_alpha_list = nullptr;
	double* __device__gain_RX_grid = nullptr;
	double* __device__pathloss_list = nullptr;

	double* ___host___hmatrix = nullptr;
	double* __device__hmatrix = nullptr;

	double* __device__polar_data_theta = nullptr;
	double* __device__polar_data_hype = nullptr;
};


using AAntennaTable = std::vector < network_package::AAntenna >;
namespace wificuda {
	void graphical_init(unsigned tx_id, const Dimensions<unsigned>& dim, AAntennaTable& anttable);
	void numerical_init(unsigned tx_id, unsigned ms_num, AAntennaTable& anttable);
	void gui_resize(unsigned tx_count, const Dimensions<unsigned>& dims);

	void graphical_update(unsigned tx_id, const Settings& current);
	void numerical_update(unsigned tx_id, const Settings& current);
	void update(unsigned tx_id, const Settings& current);

	/* Recalculate one transmitter slice without reallocating the grid buffers. */
	void graphical_recalc_polar(unsigned tx_id, const Placements& placement, const Dimensions<unsigned>& grid_size, const Settings& tx_locations);
	void numerical_recalc_polar(unsigned tx_id, const size_t tx_count, const Placements& tx_location);

	double gcoeff(unsigned tx_id, const size_t rx_sta);
	double scoeff(unsigned tx_id, const size_t rx_sta);

	//update_cow_rfs
	void handset_placement_sync(const placement_v& handset_locations);

	void gpu_teardown();
	void gpu_init();

	//void old_recalc_polar(const Dimensions<unsigned>& grid_size, const Placements& tx_location);
};

#endif
