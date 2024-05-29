#pragma once
#include <numeric>
#include "common.h"
#include "network.h"
#include "random.h"

using namespace network_package;

class Cow
{
	const unsigned station_id;
	unsigned power_idx;

	/* antenna parameters for all calculations */
	AAntenna   antenna;

	/* more antenna tracking for gui reset */
	const struct Init_Antenna
	{
		const Placements& location;
		const double& ms_grx_linear;
		const unsigned& panel_count;
		const double& lambda;
		const double& antenna_spacing;
		const double& antenna_orientation;
		const network_package::antennadim& antenna_dim;

		Init_Antenna(
			const Placements& init_bs_location,
			const double& init_ms_grx_linear,
			const unsigned& init_panel_count,
			const double& init_lambda,
			const double& init_antenna_spacing,
			const double& init_antenna_orientation,
			const network_package::antennadim& init_antenna_dim)
			:
			location(init_bs_location),
			ms_grx_linear(init_ms_grx_linear),
			panel_count(init_panel_count),
			lambda(init_lambda),
			antenna_spacing(init_antenna_spacing),
			antenna_orientation(init_antenna_orientation),
			antenna_dim(init_antenna_dim)
		{}

	} init;

	Placements location;
	Dimensions<unsigned> gui_grid_size;

	const std::vector<Placements>& ms_station_loc;
	const unsigned& ms_stations;

	/* each cell is a mobile station */
	std::vector<Polar_Coordinates> polar_data;

	/* each cell is a gui grid, hence much larger */
	std::vector<Polar_Coordinates> gui_polar_data;

public:

	/* set Gtx power in linear */
	void setPower(double& input_power)
	{
		antenna.set_power(input_power);
	}

	/* return Gtx power in linear */
	const double& getxPower() const
	{
		return antenna.get_power();
	}

	const double& getAlpha() const
	{
		return antenna.getAlpha();
	}

	const Placements& getPosition() const
	{
		return location;
	}

	const unsigned& sid() const
	{
		return station_id;
	}

	/* update the parameters of the Base Station and get channel state to each station */
	void antennaUpdate(const double& alpha)
	{
		antenna.numerical_update(alpha);
	}

	/* set the signal level in Watts (linear): second parameter */
	inline void getSignalLevel(const unsigned& node_id, double& signal_level_lin) const
	{
		signal_level_lin = antenna.coeff(node_id) * antenna.get_power();
	}


	/* polar data for gui grid rendering */
	void set_polar_data(const size_t& width, const size_t& height, const Placements& new_location, std::vector<Polar_Coordinates>& output)
	{
		size_t idx = 0;
		output.resize(width * height);

		for (int row = 0; row < width; ++row)
		{
			for (int col = 0; col < height; ++col)
			{
				long int diffx = col - new_location.x;
				long int diffy = row - new_location.y;
				output[idx++] = cart2pol(diffx, diffy);
			}
		}
	}

	/* polar data for simulation */
	void set_polar_data(const Placements& new_location, std::vector<Polar_Coordinates>& output)
	{
		size_t idx = 0;
		output.resize(ms_stations);

		for (auto& mstation : ms_station_loc)
		{
			long int diffx = mstation.x - new_location.x;
			long int diffy = mstation.y - new_location.y;
			output[idx++] = cart2pol(diffx, diffy);
		}
	}

	/* when user resets all the changes in the simulation */
	void reset()
	{
		//antenna.set_power(0);
		//antenna.set_grx_gain(init.ms_grx_linear);
		//antenna.set_antpanelcount(init.panel_count);
		//antenna.set_antlambda(init.lambda);
		//antenna.set_antspacing(init.antenna_spacing);
		//antenna.set_beamdir(init.antenna_orientation);
		//antenna.set_antdim(init.antenna_dim);
		//set_polar_data(init.location, polar_data);
		//antenna.numerical_init(polar_data);
	}

	//antennadim dim_meters, double theta, double spacing, int antenna_count,
	Cow(
		unsigned& id,
		const Placements& bs_location,
		const std::vector<Placements>& ms_pos_list,
		const double& ms_grx_linear,
		const unsigned& panel_count,
		const double& lambda,
		const double& antenna_spacing,
		const double& antenna_orientation,
		const network_package::antennadim& antenna_dim)
		:
		station_id(id),
		init(bs_location, ms_grx_linear, panel_count, lambda, antenna_spacing, antenna_orientation, antenna_dim),
		location(init.location),
		ms_station_loc(ms_pos_list),
		ms_stations(ms_pos_list.size()),
		antenna(ms_grx_linear, panel_count, lambda, antenna_spacing, antenna_orientation, antenna_dim),
		power_idx(0)
	{
		set_polar_data(init.location, polar_data);
		antenna.numerical_init(polar_data);
	}
};

/* client or handseet */
class Station
{
	unsigned station_id;
	std::vector<Cow>& cow_data; // cow that's associated to this station at any given timeslot
	const double& tnf_watt;           // thermal noise floor
	double sinr;
public:

	const unsigned& sid() const
	{
		return station_id;
	}

	/* get signal power in watts */
	inline double get_rx_signal_power(unsigned cow_id)
	{
		double power;
		cow_data[cow_id].getSignalLevel(station_id, power);
		return power;
	}

	void setRX(unsigned bs_id)
	{
		double signal = get_rx_signal_power(bs_id);

		double interference = 0;

		for (unsigned b = 0; b < cow_data.size(); ++b)
		{
			if (cow_data[b].sid() != bs_id)
			{
				auto power = get_rx_signal_power(b);
				interference += power;
			}
		}

		sinr = signal / (interference + tnf_watt);
	}

	/* get SINR in linear by passing signal and inteference/noise in dB */
	double getSINR() const
	{
		return lin2dB(sinr);
	}

	Station(unsigned id, std::vector<Cow>& cows, const double& noise) :
		station_id(id),
		cow_data(cows),
		tnf_watt(noise),
		sinr(0)
	{
		/* update the stations */
	}

};
