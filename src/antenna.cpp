#include "network.h"
using namespace network_package;

static std::vector<unsigned> indices_with_inf;
static std::vector<unsigned> indices_with_z;


/*update the antenna array from updated powerand scan angle */
void AAntenna::update(
	const size_t& malloc_size,
	double* phee_minus_alpha_list,
	double* gain_RX_grid,
	double* pathloss_list,
	double* dummy,
	double* host_hmatrix)
{
	/* update the antenna gain Gtx */
	indices_with_inf.clear();
	indices_with_z.clear();

	for (long long idx = 0; idx < malloc_size; ++idx)
	{
		if (pathloss_list[idx] == 0)
		{ // this is going to be a inf
			indices_with_inf.emplace_back(idx);
			continue;
		}
		else if (gain_RX_grid[idx] == 0)
		{ // this is going to be a zero
			indices_with_z.emplace_back(idx);
			continue;
		}

		double phee = (phee_minus_alpha_list[idx] + current.alpha) / 2;

		double sin_term = current.panel_count * sin(phee);
		double gain_factor_antenna_system = gain_RX_grid[idx]; // xN antennas already

		if (sin_term != 0)
		{
			gain_factor_antenna_system *= pow_2(sin(current.panel_count * phee) / sin_term);
		}

		/* update the channel matrix */
		host_hmatrix[idx] = gain_factor_antenna_system / pathloss_list[idx];
	}

	for (long i = 0; i < indices_with_inf.size(); ++i)
	{
		auto& problem_index = indices_with_inf[i];

		if (0 <= problem_index - 1)
		{
			host_hmatrix[problem_index] = host_hmatrix[problem_index - 1];
		}
		else if (problem_index + 1 < malloc_size)
		{
			host_hmatrix[problem_index] = host_hmatrix[problem_index + 1];
		}
		// else all of them are infinity
	}

	for (long i = 0; i < indices_with_z.size(); ++i)
	{
		auto& problem_index = indices_with_z[i];

		if (0 <= problem_index - 1)
		{
			host_hmatrix[problem_index] = host_hmatrix[problem_index - 1];
		}
		else if (problem_index + 1 < malloc_size)
		{
			host_hmatrix[problem_index] = host_hmatrix[problem_index + 1];
		}
	}
}

/* re-calc the signal outs to handsets only (before calling update!) */
void AAntenna::init(
	const size_t& malloc_size,
	double* phee_minus_alpha_list,
	double* pathloss_list,
	double* gain_RX_grid,
	const Polar_Coordinates* polar_data)
{
	const double& pioverlambda = M_PIl / current.lambda;
	const double& phee_temp = 2 * current.spacing * pioverlambda;
	const double& pl_temp_meters = 4 * pioverlambda;
	const double& antenna_dim_factor = 10 * current.antenna_dims.x * current.antenna_dims.y / pow_2(current.lambda);
	double m_factor = current.antenna_dims.x * pioverlambda;

	for (size_t idx = 0; idx < malloc_size; ++idx)
	{
		auto& cell_polar_data = polar_data[idx];

		double theta_minus_thetaC = cell_polar_data.theta - current.theta_c;
		double m = m_factor * sin(theta_minus_thetaC);
		double singleant_gain = antenna_dim_factor * pow_2((1 + cos(theta_minus_thetaC)) / 2);

		if (m != 0)
		{
			singleant_gain *= pow_2(sin(m) / m);
		}

		phee_minus_alpha_list[idx] = phee_temp * sin(theta_minus_thetaC);
		pathloss_list[idx] = pow_2(pl_temp_meters * cell_polar_data.hype);
		gain_RX_grid[idx] = singleant_gain * current.panel_count;

		if (pathloss_list[idx] == 0)
		{
			spdlog::warn("[index:" + str(idx) + "] " + str(pathloss_list[idx]) + " pathloss_list value is 0. " + \
				"Possible reason, cell_polar_data {theta,hype}:{" + str(cell_polar_data.theta) + "," + str(cell_polar_data.hype) + "} hype is zero");

		}
		else if (std::isinf(pathloss_list[idx]))
		{
			spdlog::warn("[index:" + str(idx) + "] " + str(pathloss_list[idx]) + " pathloss_list value is inf. ");
		}

		if (gain_RX_grid[idx] == 0)
		{
			spdlog::warn("[index:" + str(idx) + "] " + str(gain_RX_grid[idx]) + " gain_RX_grid value is 0. " + \
				"Possible reason, cos(theta_minus_thetaC):" + str(cos(theta_minus_thetaC)) + " + 1 = zero");

		}
		else if (std::isinf(gain_RX_grid[idx]))
		{
			spdlog::warn("[index:" + str(idx) + "] " + str(gain_RX_grid[idx]) + " gain_RX_grid value is inf. ");// +\
						//"Possible reasons cos(theta_minus_thetaC):" + str(cos(theta_minus_thetaC)) + " + 1 = zero");
		}
	}
}

void AAntenna::graphics_update()
{
	if (graphic.modified)
	{
		spdlog::info("Antenna Graphics update");
		update(graphic.hmatrix.size(), &graphic.phee_minus_alpha_list[0], &graphic.gain_RX_grid[0], &graphic.pathloss_list[0], dummy, &graphic.hmatrix[0]);
		graphic.modified = false;
	}
}

void AAntenna::numerical_update()
{
	if (simulation.modified)
	{
		spdlog::info("Antenna Numerical update");
		update(simulation.hmatrix.size(), &simulation.phee_minus_alpha_list[0], &simulation.gain_RX_grid[0], &simulation.pathloss_list[0], dummy, &simulation.hmatrix[0]);
		simulation.modified = false;
	}
}

/* for GUI simulation in the whole grid */
void AAntenna::graphics_init(PolarArray& polar_data)
{
	spdlog::info("Antenna Graphics Re-init");
	graphic.resize(polar_data.array_size);
	init(polar_data.array_size, &graphic.phee_minus_alpha_list[0], &graphic.pathloss_list[0], &graphic.gain_RX_grid[0], polar_data.data_ptr);
	graphic.modified = true;
}

/* for bare-minimum numerical calculations needed at the mobile_stations only */
void AAntenna::numerical_init(PolarArray& polar_data)
{
	spdlog::info("Antenna Numerical Re-init");
	simulation.resize(polar_data.array_size);
	init(polar_data.array_size, &simulation.phee_minus_alpha_list[0], &simulation.pathloss_list[0], &simulation.gain_RX_grid[0], polar_data.data_ptr);
	simulation.modified = true;
}

AAntenna::~AAntenna()
{
}

/* hatrix with respect to pixel index (flattened from 2D) */
const double& AAntenna::gcoeff(const unsigned& pixel_idx) const
{
	return graphic.hmatrix[pixel_idx];
}

const double& AAntenna::coeff(const unsigned& rx_sta) const
{
	return simulation.hmatrix[rx_sta];
}