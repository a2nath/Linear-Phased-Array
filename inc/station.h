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
	const Placements& init_location;
	Placements prev_location, location;
	Dimensions<size_t> init_gui_grid_size, prev_gui_grid_size, gui_grid_size;

	const std::vector<Placements>& ms_station_loc;
	const unsigned& ms_stations;

	/* each cell is a mobile station */
	std::vector<Polar_Coordinates> polar_data;

	/* each cell is a gui grid, hence much larger */
	std::vector<Polar_Coordinates> gui_polar_data;

	/* helps with gui changes and recalculate based on what the exact change */
	//struct Last_Antenna_State
	//{
	//	double& antpower_watts;
	//	double& grx_gain;
	//	double& gain_grx_linear;
	//	unsigned& panel_count;
	//	double& lambda;
	//	double& ant_spacing;
	//	double& theta_c_rads;
	//	antennadim& antdim_meters;
	//};
public:
	//void details(double& ant_spacing,
	//	unsigned& ant_panel_count,
	//	double& ant_tx_power,
	//	double& ant_scan_angle,
	//	Placements& position
	//) const
	//{
	//	ant_spacing = antenna.interAntSpacing();
	//	ant_panel_count = antenna.getPanelCount();
	//	ant_tx_power = antenna.get_power();
	//	ant_scan_angle = antenna.alpha();
	//	position = location;
	//}

	const Placements& where() const
	{
		return location;
	}

	const AAntenna& antennaDetails() const
	{
		return antenna;
	}

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

	const double& alpha() const
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
	void antenna_update(const double& alpha)
	{
		antenna.numerical_update(alpha);
	}

	/* get signal level in Watts (linear): second parameter, needs to be reshaped to M x N */
	void gui_udpate(const double& alpha)
	{
		antenna.graphics_update(alpha);
	}

	/* set the signal level in Watts (linear): second parameter */
	inline void signal_power(const unsigned& node_id, double& signal_level_lin) const
	{
		signal_level_lin = antenna.coeff(node_id) * antenna.get_power();
	}

	/* set the signal level in Watts (linear): second parameter */
	inline void g_signal_power(const size_t& pidx, double& signal_level_lin) const
	{
		signal_level_lin = antenna.gcoeff(pidx) * antenna.get_power();
	}

	void heatmap(std::vector<double>& output)
	{
		if (antenna.state_changed())
		{
			for (size_t pixel_idx = 0; pixel_idx < gui_polar_data.size(); ++pixel_idx)
			{
				g_signal_power(pixel_idx, output[pixel_idx]);
			}
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
				long int diffx = col - new_location.x; // diff with respect to pixel (think of col, row has location of rx)
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
			long int diffx = mstation.x - new_location.x; // diff with respect to rx
			long int diffy = mstation.y - new_location.y;
			output[idx++] = cart2pol(diffx, diffy);
		}
	}


	/* checks if [antenna.graphics_init] and [antenna.numerical_init] is needed based on what changed */
	void reset_gui(
		const size_t& rows,
		const size_t& cols,
		const double& new_theta_c,
		const Placements& new_location)
	{
		gui_grid_size = { rows, cols };
		location = new_location;
		antenna.rotate_cow_at(new_theta_c);

		if (location != prev_location && antenna.state_changed())
		{
			set_polar_data(rows, cols, location, gui_polar_data);
			antenna.graphics_init(gui_polar_data);

			set_polar_data(location, polar_data);   // before sinr calc. need to setup polar data on new_loc
			antenna.numerical_init(polar_data);
		}
		else if (antenna.state_changed())
		{
			antenna.graphics_init(gui_polar_data);
			antenna.numerical_init(polar_data);
		}
		else if (gui_grid_size != prev_gui_grid_size)
		{
			set_polar_data(rows, cols, location, gui_polar_data);
			antenna.graphics_init(gui_polar_data);
		}

		prev_location = location;
		prev_gui_grid_size = gui_grid_size;
	}

	/* return state is [true]=init done, else [false]=not called "init_gui" yet */
	const bool gui_state() const
	{
		return gui_polar_data.size() > 0;
	}

	/* initialize */
	void init_gui(const size_t& rows, const size_t& cols)
	{
		antenna.reset(); // antenna reset does not reset the location. See reset()

		init_gui_grid_size = { rows, cols };
		prev_gui_grid_size = init_gui_grid_size;
		gui_grid_size = init_gui_grid_size;

		set_polar_data(rows, cols, location, gui_polar_data);
		antenna.graphics_init(gui_polar_data); // always calcu
	}

	/* when user resets all the changes in the simulation */
	void reset()
	{
		antenna.reset();
		prev_location = init_location;
		location = init_location;

		set_polar_data(location, polar_data);
		antenna.numerical_init(polar_data);

		gui_polar_data.clear();
	}

	//antennadim dim_meters, double theta, double spacing, int antenna_count,
	Cow(
		unsigned& id,
		const Placements& bs_location,
		const std::vector<Placements>& ms_pos_list,
		const unsigned& panel_count,
		const double& lambda,
		const double& antenna_spacing,
		const double& antenna_orientation,
		const network_package::antennadim& antenna_dim)
		:
		station_id(id),
		init_location(bs_location),
		ms_station_loc(ms_pos_list),
		ms_stations(ms_pos_list.size()),
		antenna(panel_count, lambda, antenna_spacing, antenna_orientation, antenna_dim),
		power_idx(0)
	{
		reset();
	}
};

/* client or handseet */
class Station
{
	const unsigned  station_id;      // station id
	const double    tnf_watt;        // thermal noise floor
	const double    gain_rx;
	double sinr;
public:

	const unsigned& sid() const
	{
		return station_id;
	}

	void set_sinr(const double& calculation)
	{
		sinr = calculation;
	}

	/* get SINR in linear by passing signal and inteference/noise in linear factor */
	double get_sinr() const
	{
		return sinr;
	}

	/* get receiver gain watts */
	const double& get_grx() const
	{
		return gain_rx;
	}

	/* get system noise figure in watts */
	const double& get_nf() const
	{
		return tnf_watt;
	}

	Station(unsigned id, const double& inoise, const double& grx) :
		station_id(id),
		tnf_watt(inoise),
		gain_rx(grx),
		sinr(0)
	{}
};
using cow_v = std::vector<Cow>;
using sta_v = std::vector<Station>;
