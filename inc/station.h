#pragma once
#include <numeric>
#include "common.h"
#include "network.h"
#include "random.h"

using namespace network_package;

///* to update the GUI */
//struct station_details_t;

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
	void set_power(const double& input_power)
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

	/* get signal level in Watts (linear): second parameter, needs to be reshaped to M x N */
	inline void getHeatMapData(std::vector<double>& signal)
	{
		signal.assign(polar_data.size(), 0);
		for (size_t id = 0; id < signal.size(); ++id)
		{
			getSignalLevel(id, signal[id]);
		}
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

	/* signal heat map for cells in the GUI for one "cow" */
	void set_gui_matrix(
		Placements&  new_location,
		double&      new_antpower_watts,
		double&      new_gain_grx_linear,
		unsigned&    new_panel_count,
		double&      new_frequency, // for lambda
		double&      new_ant_spacing,
		double&      new_theta_c_rads,
		const antennadim&  new_antdim_meters,
		unsigned&    rows,
		unsigned&    cols)
	{
		/* recalculate entire antenna + grid based on the exact config change(s) */
		if (rows != gui_grid_size.x || cols != gui_grid_size.y)
		{
			set_polar_data(rows, cols, new_location, gui_polar_data);
			gui_grid_size = { rows, cols };
		}

		if (location != new_location)
		{
			set_polar_data(new_location, polar_data);
		}

		if (antenna.get_power() != new_antpower_watts)
		{
			antenna.set_power(new_antpower_watts);
		}

		if (antenna.grx_gain() != new_gain_grx_linear)
		{
			antenna.set_grx_gain(new_gain_grx_linear);
		}

		if (antenna.antpanelcount() != new_panel_count)
		{
			antenna.set_antpanelcount(new_panel_count);
		}

		double new_lambda = getLambda(new_frequency);
		if (antenna.antlambda() != new_lambda)
		{
			antenna.set_antlambda(new_lambda);
		}

		if (antenna.antspacing() != new_ant_spacing)
		{
			antenna.set_antspacing(new_ant_spacing);
		}

		if (antenna.beamdir() != new_theta_c_rads)
		{
			antenna.set_beamdir(new_theta_c_rads);
		}

		auto& dim = antenna.antdim();
		if (dim.x != new_antdim_meters.x || dim.y != new_antdim_meters.y)
		{
			antenna.set_antdim(new_antdim_meters);
		}

		antenna.graphics_init(gui_polar_data); // always calculates based on the above
	}

	/* when user resets all the changes in the simulation */
	void reset()
	{
		antenna.set_power(0);
		antenna.set_grx_gain(init.ms_grx_linear);
		antenna.set_antpanelcount(init.panel_count);
		antenna.set_antlambda(init.lambda);
		antenna.set_antspacing(init.antenna_spacing);
		antenna.set_beamdir(init.antenna_orientation);
		antenna.set_antdim(init.antenna_dim);
		set_polar_data(init.location, polar_data);
		antenna.numerical_init(polar_data);
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
	const unsigned station_id;  // id of the station in the list of stations
	std::vector<Cow>& bs_station_list;
	const double& tnf_watt;           // thermal noise floor
	double sinr;
public:

	const unsigned& sid() const
	{
		return station_id;
	}

	/* get signal power in watts */
	inline void get_rx_signal_power(const unsigned& bs_id, double& power_watts)
	{
		bs_station_list[bs_id].getSignalLevel(station_id, power_watts);
	}

	/* calcalate the rx signal based on the station-base station binding */
	void set_rx(const unsigned& associated_bs_id)
	{
		double signal, interference = 0, power;

		for (size_t bs_id = 0; bs_id < bs_station_list.size(); ++bs_id)
		{
			if (bs_id != associated_bs_id)
			{
				get_rx_signal_power(bs_id, power);
				interference += power;
			}
			else
			{
				get_rx_signal_power(associated_bs_id, power);
				signal = power;
			}
		}

		sinr = signal / (interference + tnf_watt);
	}

	/* get SINR in linear by passing signal and inteference/noise in linear factor */
	double get_sinr() const
	{
		return sinr;
	}

	Station(unsigned id, std::vector<Cow>& icow, const double& inoise) :
		station_id(id),
		bs_station_list(icow),
		tnf_watt(inoise),
		sinr(0)
	{}
};
