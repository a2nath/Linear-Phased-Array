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
	Dimensions<unsigned> init_gui_grid_size, gui_grid_size;

	const std::vector<Placements>& ms_station_loc;
	const unsigned ms_stations;

	/* each cell is a mobile station or a gui drid */
	PolarArray polar_data, gui_polar_data;

public:

	const Placements& where() const
	{
		return location;
	}

	/* only gui calls this so "update" both [sim] and [gui] components */
	void update(const Settings& new_settings, const Placements& new_location)
	{
		auto& current = antenna.settings();
		bool ant_reinit = false;
		bool ant_update = false;
		bool gui_reinit = false;

		if (current.antenna_dims != new_settings.antenna_dims)
		{
			antenna.set_antdim(new_settings.antenna_dims);
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

		if (location != new_location)
		{
			location = new_location;
			ant_reinit = true;
		}

		if (ant_reinit)
		{
			set_polar_data(location);
			antenna.numerical_init(polar_data);

			set_polar_data(gui_grid_size, location);
			antenna.graphics_init(gui_polar_data);
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

	/* set the signal level in Watts (linear): second parameter */
	inline void g_signal_power(const size_t& pidx, double& signal_level_lin, const bool& debug)
	{
		signal_level_lin = antenna.gcoeff(pidx) * antenna.get_power();
	}

	void heatmap(std::vector<double>& output, const bool& debug = false)
	{
		for (size_t pixel_idx = 0; pixel_idx < gui_polar_data.array_size; ++pixel_idx)
		{
			g_signal_power(pixel_idx, output[pixel_idx], debug);
		}
	}

	/* polar data for gui grid rendering */
	void set_polar_data(const Dimensions<unsigned>& new_size, const Placements& new_location)
	{
		size_t idx = 0;
		gui_polar_data.set(new_size.x * new_size.y);

		for (int row = 0; row < new_size.y; ++row)
		{
			for (int col = 0; col < new_size.x; ++col)
			{
				long int diffx = col - new_location.x; // diff with respect to pixel (think of col, row has location of rx)
				long int diffy = row - new_location.y;
				gui_polar_data.data_ptr[idx++] = cart2pol(diffx, diffy);
			}
		}
	}

	/* polar data for simulation */
	void set_polar_data(const Placements& new_location)
	{
		size_t idx = 0;
		polar_data.set(ms_stations);

		for (auto& mstation : ms_station_loc)
		{
			long int diffx = mstation.x - new_location.x; // diff with respect to rx
			long int diffy = mstation.y - new_location.y;
			polar_data.data_ptr[idx++] = cart2pol(diffx, diffy);
		}
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
		state.location = location;

		return state;
	}

	/* return state is [true]=init done, else [false]=not called "init_gui" yet */
	const bool gui_ready() const
	{
		return gui_polar_data.array_size > 0;
	}

	/* resize the gui window;
	NOTE: the [second] part of this function cannot be DELAYED further */
	void resize_gui(const Dimensions<unsigned>& new_dimension)
	{
		if (gui_ready() && !new_dimension.is_zero() && gui_grid_size != new_dimension)
		{
			gui_grid_size = new_dimension;
			set_polar_data(gui_grid_size, location);
			antenna.graphics_init(gui_polar_data);
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

		set_polar_data(gui_grid_size, location);
		antenna.graphics_init(gui_polar_data);
	}

	void undo()
	{
		antenna.undo();
		std::swap(prev_location, location);

		set_polar_data(location);
		antenna.numerical_init(polar_data);
	}

	/* when user resets all the changes in the simulation */
	void init_sim()
	{
		antenna.reset();
		prev_location = init_location;
		location = init_location;

		set_polar_data(location);
		antenna.numerical_init(polar_data);
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
		const antennadim& antenna_dim)
		:
		station_id(id),
		init_location(bs_location),
		ms_station_loc(ms_pos_list),
		ms_stations(ms_pos_list.size()),
		antenna(panel_count, lambda, antenna_spacing, antenna_orientation, antenna_dim),
		power_idx(0)
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
	{}
};
using cow_v = std::vector<Cow>;
using sta_v = std::vector<Station>;
