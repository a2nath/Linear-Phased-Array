#pragma once
#include <cmath>
#include <smmintrin.h>
#include <immintrin.h>
#include "coordinates.h"
#include "common.h"

namespace network_package
{
	using antennadim = Dimensions<double>;

	inline double rad2deg(double rad)
	{
		return rad * 180.0 / M_PIl;
	}
	inline double deg2rad(double deg)
	{
		return deg * M_PIl / 180.0;
	}
	inline double log2lin(double log)
	{
		return pow(10, log / 10);
	}
	inline double lin2dB(double lin)
	{
		return 10 * log10(lin);
	}
	inline double dBm2watt(double dBm)
	{
		return log2lin(dBm - 30);
	}

	inline double getLambda(double frequency)
	{
		if (frequency > 0)
		{
			return C_SPEED / frequency;
		}

		throw std::invalid_argument("Divide by zero error from passing 0 frequency in getLambda call");
	}

	/* Get system noise figure in [dBm] with thermal noise by passing B/W in [Mhz], and system noise in [dB] */
	inline double getThermalSystemNoise(const double& bandwidth, const double& system_noise)
	{
		if (bandwidth > 0)
		{
			return -174 + round(10 * log10(bandwidth)) + system_noise;
		}

		throw std::invalid_argument("Bandwidth cannot be zero or when setting system noise factor");
	}

	/* Linear Phase Array Antenna */
	class AAntenna
	{
		//typedef void (*FuncPtr)(double);
		double     power;              // array power
		double     alpha;              // scan angle
		double     ms_grx_linear;      // gain factor for a client station
		unsigned   panel_count;        // number of panels
		double     lambda;             // wavelength of the RF signal
		double     spacing;            // wavelength of the RF signal
		double     theta_c;            // antenna array direction
		antennadim antenna_dims;

		struct Calculations
		{
			/* Gtx gain for the whole grid */
			std::vector<double> gain_RX_grid;
			std::vector<double> phee_minus_alpha_list;
			std::vector<double> pathloss_list;
			std::vector<double> hmatrix;                 // hmatrix from BS pov

			void resize(const size_t& size)
			{
				gain_RX_grid.resize(size);
				phee_minus_alpha_list.resize(size);
				pathloss_list.resize(size);
				hmatrix.resize(size);
			}
		};

		/* Calculations needed to create/update the H-matrix coefficient table */
		Calculations simulation, graphic;

	public:
		/* accepts a RX station ID with respect to THIS station
		   and returns the coefficient
		*/
		const double& coeff(const unsigned& rx_sta) const
		{
			return simulation.hmatrix[rx_sta];
		}

		/* set power in watts */
		void set_power(const double& power_watts)
		{
			power = power_watts;
		}

		/* get power in watts */
		const double& get_power() const
		{
			return power;
		}

		/* get alpha in rads */
		const double& getAlpha() const
		{
			return alpha;
		}

		/* set station grx in linear_factor */
		void set_grx_gain(const unsigned& linear_factor)
		{
			ms_grx_linear = linear_factor;
		}

		/* get station grx in linear_factor */
		const double& grx_gain() const
		{
			return ms_grx_linear;
		}
		/* sett panel count in the antenna array */
		void set_antpanelcount(const unsigned& count)
		{
			panel_count = count;
		}

		/* get panel count in the antenna array */
		const unsigned& antpanelcount() const
		{
			return panel_count;
		}

		/* set wavelength in meters */
		void set_antlambda(const double& meters_lambda)
		{
			lambda = meters_lambda;
		}

		/* get wavelength in meters */
		const double& antlambda() const
		{
			return lambda;
		}

		/* set the physical antenna panel spacing in meters */
		void set_antspacing(const double& meters_separation)
		{
			spacing = meters_separation;
		}

		/* get the physical antenna panel spacing in meters */
		const double& antspacing() const
		{
			return spacing;
		}

		/* get the physical antenna direction in rads */
		void set_beamdir(const double& rads_direction)
		{
			theta_c = rads_direction;
		}

		/* get the physical antenna direction */
		const double& beamdir() const
		{
			return theta_c;
		}

		/* set the physical size in meters */
		void set_antdim(const antennadim& meters_dim)
		{
			antenna_dims = meters_dim;
		}

		/* get the physical size in meters */
		const antennadim antdim() const
		{
			return antenna_dims;
		}

		/* update the antenna array from updated power and scan angle */
		inline void update(const double& new_alpha, Calculations& calculations)
		{
			if (new_alpha != alpha)
			{
				/* update the antenna gain Gtx */
				for (unsigned idx = 0; idx < calculations.phee_minus_alpha_list.size(); ++idx)
				{
					double phee = (calculations.phee_minus_alpha_list[idx] + new_alpha) / 2;

					double sin_term = panel_count * sin(phee);
					double gain_factor_antenna_system = calculations.gain_RX_grid[idx]; // xN antennas already

					if (sin_term != 0)
					{
						gain_factor_antenna_system *= pow(cached::sin(panel_count * phee) / sin_term, 2);
					}

					/* update the channel matrix */
					calculations.hmatrix[idx] = ms_grx_linear * gain_factor_antenna_system / calculations.pathloss_list[idx];
				}

				/* update the scan angle of the antenna array */
				alpha = new_alpha;
			}
		}

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void graphics_update(const double& new_alpha)
		{
			update(new_alpha, graphic);
		}

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void numerical_update(const double& new_alpha)
		{
			update(new_alpha, simulation);
		}

		/* re-calc the signal outs to handsets only (before calling update!) */
		inline void init(const std::vector<Polar_Coordinates>& polar_data, Calculations& calculations)
		{
			const double& pioverlambda = M_PIl / lambda;
			const double& phee_temp = 2 * spacing * pioverlambda;
			const double& pl_temp_meters = 4 * pioverlambda;
			const double& antenna_dim_factor = 10 * antenna_dims.x * antenna_dims.y / pow(lambda, 2);
			double m_factor = antenna_dims.x * pioverlambda;

			for (size_t idx = 0; idx < polar_data.size(); ++idx)
			{
				auto& cell_polar_data = polar_data[idx];

				double theta_minus_thetaC = cell_polar_data.theta - theta_c;
				double m = m_factor * cached::sin(theta_minus_thetaC);
				double singleant_gain = antenna_dim_factor * pow((1 + cached::cos(theta_minus_thetaC)) / 2, 2);

				if (m != 0)
				{
					singleant_gain *= pow(cached::sin(m) / m, 2);
				}

				calculations.phee_minus_alpha_list[idx] = phee_temp * cached::sin(theta_minus_thetaC);
				calculations.pathloss_list[idx]         = pow(pl_temp_meters * cell_polar_data.hype, 2);
				calculations.gain_RX_grid[idx]          = singleant_gain * panel_count;
			}
		}

		/* for GUI simulation in the whole grid */
		void graphics_init(const std::vector<Polar_Coordinates>& polar_data)
		{
			graphic.resize(polar_data.size());
			init(polar_data, graphic);
		}

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void numerical_init(const std::vector<Polar_Coordinates>& polar_data)
		{
			simulation.resize(polar_data.size());
			init(polar_data, simulation);
		}

		AAntenna(
			const double& init_ms_grx_linear,
			const unsigned& init_panel_count,
			const double& init_lambda,
			const double& init_antenna_spacing,
			const double& init_antenna_orientation_rads,
			const antennadim& init_antdims)
			:
			power(0),
			alpha(0),
			ms_grx_linear(init_ms_grx_linear),
			panel_count(init_panel_count),
			lambda(init_lambda),
			spacing(init_antenna_spacing),
			theta_c(init_antenna_orientation_rads),
			antenna_dims(init_antdims)
		{}
	};
};
