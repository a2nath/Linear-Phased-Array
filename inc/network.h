#pragma once
#include <cmath>
#include <smmintrin.h>
#include <immintrin.h>
#include "coordinates.h"
#include "common.h"

namespace network_package
{
	using namespace rf_math;

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
			double* host_hmatrix;

			bool modified;

			void resize(const size_t& size)
			{
				gain_RX_grid.resize(size);
				phee_minus_alpha_list.resize(size);
				pathloss_list.resize(size);
				hmatrix.resize(size);
			}

			Calculations() : modified(false), host_hmatrix(nullptr) {}
		};

		/* Calculations needed to create/update the H-matrix coefficient table */
		Calculations simulation, graphic;
		double* dummy;

	public:

		const Settings& settings() const
		{
			return current;
		}

		/* accepts a RX station ID with respect to THIS station
		   and returns the coefficient
		*/
		const double& coeff(const unsigned& rx_sta) const;

		/* hatrix with respect to pixel index (flattened from 2D) */
		const double& gcoeff(const unsigned& pixel_idx) const;

		/* set power in watts */
		void set_power(const double& power_watts)
		{
			prev.power = current.power;
			current.power = power_watts;

			simulation.modified = true;
			graphic.modified = true;
		}

		/* get power in watts */
		float& get_power()
		{
			return current.power;
		}

		/* get alpha in rads */
		float& getAlpha()
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
		float& antlambda()
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
		float& antspacing()
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
		float& beamdir()
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
		inline void update(
			const size_t& malloc_size,
			double* phee_minus_alpha_list,
			double* gain_RX_grid,
			double* pathloss_list,
			double* gpu_hmatrix,
			double* host_hmatrix);

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void graphics_update();

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void numerical_update();

		/* re-calc the signal outs to handsets only (before calling update!) */
		inline void init(
			const size_t& malloc_size,
			double* d_phee_minus_alpha_list,
			double* d_pathloss_list,
			double* d_gain_RX_grid,
			const Polar_Coordinates* d_polar_data);

		/* for GUI simulation in the whole grid */
		void graphics_init(PolarArray& polar_info);

		/* for bare-minimum numerical calculations needed at the mobile_stations only */
		void numerical_init(PolarArray& polar_info);

		/* change the state back to the previous one */
		void undo()
		{
			std::swap(prev, current);
		}

		void reset()
		{
			prev = initial;
			current = initial;
		}

		~AAntenna();

		AAntenna(
			const unsigned& init_panel_count,
			const double& init_lambda,
			const double& init_antenna_spacing,
			const double& init_antenna_orientation_rads,
			const antennadim& init_antdims)
			:
			initial{ 0,
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
	};

};
