#pragma once
#include <numeric>
#include "common.h"
#include "network.h"
#include "random.h"

using namespace network_package;
using AAntennaTable = std::vector<AAntenna>;

///* to update the GUI */
//struct station_details_t;

class Cow
{
	const unsigned station_id;
	unsigned power_idx;

	/* antenna parameters for all calculations */
	AAntennaTable& antennatable;
	AAntenna& antenna;

	/* more antenna tracking for gui reset */
	const Placements init_location;
	Placements prev_location, location;
	Dimensions<unsigned> init_gui_grid_size, gui_grid_size;

	const std::vector<Placements>& ms_station_loc;
	const unsigned ms_stations_num;

public:

	/* only gui calls this so "update" both [sim] and [gui] components */
	void update(const Settings& new_settings, const Placements& new_location)
	{
		auto& current = antenna.settings();
		bool ant_reinit = false;
		bool ant_update = false;
		//bool gui_reinit = false;

		if (current.antenna_dims != new_settings.antenna_dims)
		{
			antenna.set_antdim(new_settings.antenna_dims);
			ant_reinit = true;
		}

		if (location != new_location)
		{
			antenna.set_location(new_settings.location);
			ant_reinit = true;
		}

		if (current.lambda != new_settings.lambda)
		{
			antenna.set_antlambda(new_settings.lambda);
			ant_reinit = true;
		}

		if (current.panel_count != new_settings.panel_count)
		{
			antenna.set_antpanelcount(new_settings.panel_count);
			ant_reinit = true;
		}

		if (current.spacing != new_settings.spacing)
		{
			antenna.set_antspacing(new_settings.spacing);
			ant_reinit = true;
		}

		if (current.theta_c != new_settings.theta_c)
		{
			antenna.rotate_cow_at(new_settings.theta_c);
			ant_reinit = true;
		}

		if (current.alpha != new_settings.alpha)
		{
			antenna.set_alpha(new_settings.alpha);
			ant_update = true;
		}

		if (current.power != new_settings.power)
		{
			antenna.set_power(new_settings.power);
		}

		if (ant_reinit)
		{
			antenna.numerical_init(ms_station_loc);

			antenna.graphics_init(gui_grid_size.x, gui_grid_size.y);
		}

		if (ant_reinit || ant_update)
		{
			antenna.numerical_update();
			antenna.graphics_update();
		}
	}

	void update_minimal(const double& power, const double& alpha)
	{
		antenna.set_power(power);
		antenna.set_alpha(alpha);
		antenna.numerical_update();

		if (gui_ready()) // minimal update during start of sim or visualization
		{
			antenna.graphics_update();
		}
	}

	/* set Gtx power in linear */
	void set_power(const double& input_power)
	{
		antenna.set_power(input_power);
	}

	/* return Gtx power in linear */
	float& getxPower()
	{
		return antenna.get_power();
	}

	float& alpha()
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

	/* set the signal level in Watts (linear): second parameter */
	inline void signal_power(const unsigned& node_id, double& signal_level_lin)
	{
		signal_level_lin = antenna.coeff(node_id) * antenna.get_power();
	}

	inline void heatmap(std::vector<double>& output, const bool& debug = false) const
	{
		antenna.get_signal_power(output, debug);
	}

	const std::string str() const
	{
		return "cow:" + std::to_string(station_id) + ", location:" + location.str() + ", antenna settings:" + antenna.settings().str();
	}

	/* state of the TX station */
	graphics::State get_state()
	{
		graphics::State state(station_id);
		state.settings = antenna.settings();

		return state;
	}

	/* return state is [true]=init done, else [false]=not called "init_gui" yet */
	bool gui_ready() const
	{
		return gui_grid_size.count() > 0;
	}

	/* resize the gui window;
	NOTE: the [second] part of this function cannot be DELAYED further */
	void resize_gui(const Dimensions<unsigned>& new_dimension)
	{
		if (gui_ready() && !new_dimension.is_zero() && gui_grid_size != new_dimension)
		{
			gui_grid_size = new_dimension;
			antenna.graphics_init(gui_grid_size.x, gui_grid_size.y);
			antenna.graphics_update();
		}
		else if (!gui_ready())
		{
			spdlog::error("GUI not ready when resizing");
			throw std::runtime_error("GUI needs initializing");
		}
		else if (new_dimension.is_zero())
		{
			spdlog::error("New size is not valid");
			throw std::runtime_error("Enter valid size");
		}
		else
		{
			spdlog::warn("performance concern: calling resize() more than once with the same parameters");
		}
	}

	/* when user resets or initializes data structures GUI */
	void init_gui(const unsigned& rows, const unsigned& cols)
	{
		if (rows == 0 || cols == 0)
		{
			spdlog::error("New size is not valid");
			throw std::runtime_error("Enter valid size");
		}

		init_gui_grid_size = { rows, cols };
		gui_grid_size = init_gui_grid_size;

		antenna.graphics_init(gui_grid_size.x, gui_grid_size.y);
	}

	void undo()
	{
		antenna.undo();
		std::swap(prev_location, location);

		antenna.numerical_init(ms_station_loc);
	}

	/* when user resets all the changes in the simulation */
	void init_sim()
	{
		antenna.reset();
		prev_location = init_location;
		location = init_location;

		antenna.numerical_init(ms_station_loc);
	}

	//antennadim dim_meters, double theta, double spacing, int antenna_count,
	Cow(
		unsigned& id,
		const std::vector<Placements>& ms_pos_list,
		AAntennaTable& antenna_table)
		:
		station_id(id),
		power_idx(0),
		antennatable(antenna_table),
		antenna(antenna_table.at(id)),
		init_location(antenna.settings().location),
		ms_station_loc(ms_pos_list),
		ms_stations_num(ms_station_loc.size())
	{
		init_sim();
	}
};

/* client or handseet */
class Station
{
	const unsigned  station_id;      // station id
	const double&   tnf_watt;        // thermal noise floor
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

	/* get system noise factor in watts */
	const double& get_nf() const
	{
		return tnf_watt;
	}

	Station(unsigned id, const double& inoise, const double& grx) :
		station_id(id),
		tnf_watt(inoise),
		gain_rx(grx),
		sinr(0)
	{
	}
};
using cow_v = std::vector<Cow>;
using sta_v = std::vector<Station>;
