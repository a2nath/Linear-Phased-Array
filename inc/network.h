#pragma once
#include <cmath>
#include <smmintrin.h>
#include <immintrin.h>
#include "coordinates.h"
#include "common.h"

namespace network_package
{
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
	inline double watt2dBm(double lin)
	{
		return lin2dB(lin) + 30;
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

	/* keep track of sinr and note the configurations and bindings in decreasing order of SINR */
	struct dataitem_t
	{
		unsigned tslot_id;
		double power;
		double alpha;
		Placements pos;
		unsigned cow_id;
		unsigned sta_id;
		double sinr;

		friend std::ostream& operator<<(std::ostream& out, const dataitem_t& item)
		{
			out
				<< std::fixed
				<< std::setprecision(2)
				<< std::setprecision(2)
				<< std::fixed
				<< std::fixed
				<< std::fixed
				<< std::fixed
				<< std::setprecision(2);
			out
				<< "timeslot: " << item.tslot_id << "\t"
				<< "power: " << watt2dBm(item.power) << "\t"
				<< "alpha: " << rad2deg(item.alpha) << "\t"
				<< "placement: " << item.pos.x << ","
				<< item.pos.y << "\t"
				<< "cow_id: " << item.cow_id << "\t"
				<< "sta_id: " << item.sta_id << "\t"
				<< "sinr: " << lin2dB(item.sinr);

			return out;
		}

		// Deleted copy constructor and copy assignment operator
		dataitem_t(const dataitem_t&) = delete;
		dataitem_t& operator=(const dataitem_t&) = delete;

		// Defaulted move constructor and move assignment operator
		dataitem_t(dataitem_t&&) noexcept = default;
		dataitem_t& operator=(dataitem_t&& other) noexcept = default;

		dataitem_t(
			const unsigned& itslot_id,
			const double& ipower,
			const double& ialpha,
			const Placements& ipos,
			const unsigned& icow_id,
			const unsigned& ista_id,
			const double& isinr)
			:
			tslot_id(itslot_id),
			power(ipower),
			alpha(ialpha),
			pos(ipos),
			cow_id(icow_id),
			sta_id(ista_id),
			sinr(isinr)
		{}
	};

	struct data_comparator
	{
		bool operator()(const dataitem_t& p1, const dataitem_t& p2) const
		{
			return p1.sinr < p2.sinr; // descending order
		}
	};


	/* Linear Phase Array Antenna */
	class AAntenna
	{
		const Settings initial;
		Settings prev, current;

		struct Calculations
		{
			/* Gtx gain for the whole grid */
			std::vector<double> gain_RX_grid;
			std::vector<double> phee_minus_alpha_list;
			std::vector<double> pathloss_list;
			std::vector<double> hmatrix;                 // hmatrix from BS pov
			bool modified;

			void resize(const size_t& size)
			{
				gain_RX_grid.resize(size);
				phee_minus_alpha_list.resize(size);
				pathloss_list.resize(size);
				hmatrix.resize(size);

				modified = true;
			}

			Calculations() : modified(false) {}
		};

		/* Calculations needed to create/update the H-matrix coefficient table */
		Calculations simulation, graphic;

		std::vector<unsigned> indices_with_inf;
		std::vector<unsigned> indices_with_z;

	public:

		const Settings& settings() const
		{
			return current;
		}

		/* accepts a RX station ID with respect to THIS station
		   and returns the coefficient
		*/
		const double& coeff(const unsigned& rx_sta) const
		{
			return simulation.hmatrix[rx_sta];
		}

		/* hatrix with respect to pixel index (flattened from 2D) */
		const double& gcoeff(const unsigned& pixel_idx) const
		{
			return graphic.hmatrix[pixel_idx];
		}

		/* set power in watts */
		void set_power(const double& power_watts)
		{
			prev.power = current.power;
			current.power = power_watts;

			graphic.modified = true;
		}

		/* get power in watts */
		const double& get_power() const
		{
			return current.power;
		}

		/* get alpha in rads */
		const double& getAlpha() const
		{
			return current.alpha;
		}

		/* sett panel count in the antenna array */
		void set_antpanelcount(const unsigned& count)
		{
			prev.panel_count = current.panel_count;
			current.panel_count = count;
		}

		/* get panel count in the antenna array */
		const unsigned& antpanelcount() const
		{
			return current.panel_count;
		}

		/* set wavelength in meters */
		void set_antlambda(const double& meters_lambda)
		{
			prev.lambda = current.lambda;
			current.lambda = meters_lambda;
		}

		/* get wavelength in meters */
		const double& antlambda() const
		{
			return current.lambda;
		}

		/* set the physical antenna panel spacing in meters */
		void set_antspacing(const double& meters_separation)
		{
			prev.spacing = current.spacing;
			current.spacing = meters_separation;
		}

		/* get the physical antenna panel spacing in meters */
		const double& antspacing() const
		{
			return current.spacing;
		}

		/* set the physical antenna direction in rads */
		void rotate_cow_at(const double& rads_direction)
		{
			prev.theta_c = current.theta_c;
			current.theta_c = rads_direction;
		}

		/* get the physical antenna direction */
		const double& beamdir() const
		{
			return current.theta_c;
		}

		/* set the physical size in meters */
		void set_antdim(const antennadim& meters_dim)
		{
			prev.antenna_dims = current.antenna_dims;
			current.antenna_dims = meters_dim;
		}

		/* get the physical size in meters */
		const antennadim antdim() const
		{
			return current.antenna_dims;
		}

		void set_alpha(const double& dir_rads)
		{
			prev.alpha = current.alpha;
			current.alpha = dir_rads;

			simulation.modified = true;
			graphic.modified = true;
		}

		/*update the antenna array from updated powerand scan angle */
		inline void update(Calculations& calculations)
		{
			/* update the antenna gain Gtx */
			for (long long idx = 0; idx < calculations.phee_minus_alpha_list.size(); ++idx)
			{
				if (calculations.pathloss_list[idx] == 0)
				{ // this is going to be a inf
					indices_with_inf.emplace_back(idx);
					continue;
				}
				else if (calculations.gain_RX_grid[idx] == 0)
				{ // this is going to be a zero
					indices_with_z.emplace_back(idx);
					continue;
				}

				double phee = (calculations.phee_minus_alpha_list[idx] + current.alpha) / 2;

				double sin_term = current.panel_count * sin(phee);
				double gain_factor_antenna_system = calculations.gain_RX_grid[idx]; // xN antennas already

				if (sin_term != 0)
				{
					gain_factor_antenna_system *= cached::pow_2(cached::sin(current.panel_count * phee) / sin_term);
				}

				/* update the channel matrix */
				calculations.hmatrix[idx] = gain_factor_antenna_system / calculations.pathloss_list[idx];
			}

			for (long i = 0; i < indices_with_inf.size(); ++i)
			{
				auto& problem_index = indices_with_inf[i];

				if (0 <= problem_index - 1)
				{
					calculations.hmatrix[problem_index] = calculations.hmatrix[problem_index - 1];
				}
				else if (problem_index + 1 < calculations.hmatrix.size())
				{
					calculations.hmatrix[problem_index] = calculations.hmatrix[problem_index + 1];
				}
				// else all of them are infinity
			}

			for (long i = 0; i < indices_with_z.size(); ++i)
			{
				auto& problem_index = indices_with_z[i];

				if (0 <= problem_index - 1)
				{
					calculations.hmatrix[problem_index] = calculations.hmatrix[problem_index - 1];
				}
				else if (problem_index + 1 < calculations.hmatrix.size())
				{
					calculations.hmatrix[problem_index] = calculations.hmatrix[problem_index + 1];
				}
			}
		}

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void graphics_update()
		{
			if (graphic.modified)
			{
				update(graphic);
				graphic.modified = false;
			}
		}

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void numerical_update()
		{
			if (simulation.modified)
			{
				update(simulation);
				simulation.modified = false;
			}
		}

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void numerical_update(const double& new_alpha)
		{
			prev.alpha = current.alpha;
			current.alpha = new_alpha;

			numerical_update();
		}

		/* re-calc the signal outs to handsets only (before calling update!) */
		inline void init(const std::vector<Polar_Coordinates>& polar_data, Calculations& calculations)
		{
			//if (current.theta_c != new_theta_c)
			//{
			const double& pioverlambda = M_PIl / current.lambda;
			const double& phee_temp = 2 * current.spacing * pioverlambda;
			const double& pl_temp_meters = 4 * pioverlambda;
			const double& antenna_dim_factor = 10 * current.antenna_dims.x * current.antenna_dims.y / cached::pow_2(current.lambda);
			double m_factor = current.antenna_dims.x * pioverlambda;

			for (size_t idx = 0; idx < polar_data.size(); ++idx)
			{
				auto& cell_polar_data = polar_data[idx];

				double theta_minus_thetaC = cell_polar_data.theta - current.theta_c;
				double m = m_factor * cached::sin(theta_minus_thetaC);
				double singleant_gain = antenna_dim_factor * cached::pow_2((1 + cached::cos(theta_minus_thetaC)) / 2);

				if (m != 0)
				{
					singleant_gain *= cached::pow_2(cached::sin(m) / m);
				}

				calculations.phee_minus_alpha_list[idx] = phee_temp * cached::sin(theta_minus_thetaC);
				calculations.pathloss_list[idx] = cached::pow_2(pl_temp_meters * cell_polar_data.hype);
				calculations.gain_RX_grid[idx] = singleant_gain * current.panel_count;

				if (calculations.pathloss_list[idx] == 0)
				{
					spdlog::warn("[index:" + str(idx) + "] " + str(calculations.pathloss_list[idx]) + " pathloss_list value is 0. " + \
						"Possible reason, cell_polar_data {theta,hype}:{" + str(cell_polar_data.theta) + "," + str(cell_polar_data.hype) + "} hype is zero");

				}
				else if (std::isinf(calculations.pathloss_list[idx]))
				{
					spdlog::warn("[index:" + str(idx) + "] " + str(calculations.pathloss_list[idx]) + " pathloss_list value is inf. ");
				}

				if (calculations.gain_RX_grid[idx] == 0)
				{
					spdlog::warn("[index:" + str(idx) + "] " + str(calculations.gain_RX_grid[idx]) + " gain_RX_grid value is 0. " + \
						"Possible reason, cos(theta_minus_thetaC):" + str(cos(theta_minus_thetaC)) + " + 1 = zero");

				}
				else if (std::isinf(calculations.gain_RX_grid[idx]))
				{
					spdlog::warn("[index:" + str(idx) + "] " + str(calculations.gain_RX_grid[idx]) + " gain_RX_grid value is inf. ");// +\
						//"Possible reasons cos(theta_minus_thetaC):" + str(cos(theta_minus_thetaC)) + " + 1 = zero");
				}
			}

			//current.theta_c = new_theta_c;
		//}
		}

		/* for GUI simulation in the whole grid */
		void graphics_init(const std::vector<Polar_Coordinates>& polar_data)
		{
			graphic.resize(polar_data.size());
			init(polar_data, graphic);
			graphic.modified = true;
		}

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void numerical_init(const std::vector<Polar_Coordinates>& polar_data)
		{
			simulation.resize(polar_data.size());
			init(polar_data, simulation);
			simulation.modified = true;
		}

		/* return if antenna phy parameters changed */
		const bool modified() const
		{
			return simulation.modified || graphic.modified;
		}

		/* change the state back to the previous one */
		void undo()
		{
			std::swap(prev, current);
		}

		void reset()
		{
			prev = current;
			current = initial;
		}

		AAntenna(
			const unsigned& init_panel_count,
			const double& init_lambda,
			const double& init_antenna_spacing,
			const double& init_antenna_orientation_rads,
			const antennadim& init_antdims)
			:
			initial{ 0, std::numeric_limits<double>::min(), init_panel_count, init_lambda, init_antenna_spacing, init_antenna_orientation_rads, init_antdims } // constant initial setup
		{
		}
	};

};
