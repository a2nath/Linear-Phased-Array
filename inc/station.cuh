#pragma once
#include <numeric>
#include "common.h"
#include "network.cuh"
#include "random.h"
#include "cuda.cuh"
#include "undo.h"

using namespace network_package;
#define UNDO_BUFFER_SIZE     50

class Cow
{
	const unsigned station_id;
	unsigned power_idx;

	/* antenna parameters for all calculations */
	AAntennaTable& antennatable;
	AAntenna& antenna;

	/* more antenna tracking for gui reset */
	const Placements init_location;
	Placements location;
	UndoRedoBuffer<Dimensions<unsigned>> undoredo_guisize_tracker;

	const std::vector<Placements>& ms_station_loc;
	const unsigned ms_stations_num;

	/* Support undo operations for undo 50 actions
	* and store undo actions and redo */
	UndoRedoBuffer<Settings> undoredo_settings_tracker;

	double* host_simdata_coefficients_ptr;
	double* host_guidata_coefficients_ptr;
public:
	/* only gui calls this so "update" both [sim] and [gui] components */
	void update(const Settings& new_settings)
	{
		const Settings current = antenna.settings();
		const Settings& next = antenna.settings();

		bool ant_reinit = false;
		bool ant_update = false;
		//bool gui_reinit = false;

		if (next.antenna_dims != new_settings.antenna_dims)
		{
			antenna.set_antdim(new_settings.antenna_dims);
		}

		if (next.location != new_settings.location)
		{
			antenna.set_location(new_settings.location);
		}

		if (next.lambda != new_settings.lambda)
		{
			antenna.set_antlambda(new_settings.lambda);
		}

		if (next.panel_count != new_settings.panel_count)
		{
			antenna.set_antpanelcount(new_settings.panel_count);
		}

		if (next.spacing != new_settings.spacing)
		{
			antenna.set_antspacing(new_settings.spacing);
		}

		if (next.theta_c != new_settings.theta_c)
		{
			antenna.rotate_cow_at(new_settings.theta_c);
		}

		if (next.alpha != new_settings.alpha)
		{
			antenna.set_alpha(new_settings.alpha);
			ant_update = true;

			sim_update = true;
			gfx_update = true;
		}

		if (next.power != new_settings.power)
		{
			antenna.set_power(new_settings.power);

			sim_update = true;
			gfx_update = true;
		}

		update_changes(current, next);
		undoredo_settings_tracker.emplace(antenna.settings());
	}

	/* If COW is shifted anyone on the map or graphics need update */
	void update_changes(const Settings& current, const Settings& next)
	{
		if (current != next)
		{
			if (antenna.full_re_init)
			{
				/* the location of the COW has been changed */
				auto& gui_grid_size = undoredo_guisize_tracker.get_current();

				wificuda::numerical_recalc_polar(
					station_id, antennatable.size(), next.location);
				wificuda::numerical_init(station_id, ms_stations_num, antennatable);
				wificuda::graphical_recalc_polar(station_id, next.location, gui_grid_size, next);
				wificuda::graphical_init(station_id, gui_grid_size, antennatable);

				antenna.full_re_init = false;
				location = next.location;
			}

			// compare any other characteristics.
			wificuda::numerical_update(station_id, next);
			wificuda::graphical_update(station_id, next);
		}
	}

	/* Inputs {antenna-power(watt), antenna-scan-angle(deg)} */
	void update_antenna_rf(const double& power, const double& alpha)
	{
		antenna.set_power(power);
		antenna.set_alpha(alpha);

		sim_update = true;
		gfx_update = true;

		wificuda::numerical_update(station_id, antenna.settings());
	}

	/* set Gtx power in linear */
	void set_alpha(const double& angle_rad)
	{
		antenna.set_alpha(angle_rad);

		sim_update = true;
		gfx_update = true;

		wificuda::update(station_id, antenna.settings());
	}


	/* set Gtx power in linear */
	void set_power(const double& input_power)
	{
		antenna.set_power(input_power);

		sim_update = true;
		gfx_update = true;

		wificuda::update(station_id, antenna.settings());
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
		host_simdata_coefficients_ptr = wificuda::scoeff_ptr(station_id);

		if (!host_simdata_coefficients_ptr)
		{
			spdlog::error("Did not call CUDA INIT first. Data pointer is null");
		}

		signal_level_lin = host_simdata_coefficients_ptr[node_id] * antenna.get_power();
	}

	void heatmap(std::vector<double>& output, const bool& debug = false)
	{
		/* get the hmatrix or coefficients from the device to host */
		host_guidata_coefficients_ptr = wificuda::gcoeff_ptr(station_id);

		if (!host_guidata_coefficients_ptr)
		{
			spdlog::error("Did not call CUDA INIT first. Data pointer is null");
		}

		/* get the cell count from the current setting */
		const auto& cells_num = undoredo_guisize_tracker.get_current().count();

		for (size_t pixel_idx = 0; pixel_idx < cells_num; ++pixel_idx)
		{
			output[pixel_idx] = host_guidata_coefficients_ptr[pixel_idx] * antenna.get_power();
		}
	}

	const std::string str() const
	{
		return "cow:" + std::to_string(station_id) + ", location:" + (location).str() + ", antenna settings:" + antenna.settings().str();
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
		return undoredo_guisize_tracker.get_current().count() > 0;
	}

	/* resize the gui window;
	NOTE: the [second] part of this function cannot be DELAYED further */
	void resize_gui(const Dimensions<unsigned>& new_dimension) // -> resize_event should be more high level
	{
		auto gui_grid_size = undoredo_guisize_tracker.get_current();

		if (!new_dimension.is_zero() && gui_grid_size != new_dimension)
		{
			undoredo_guisize_tracker.emplace(new_dimension);
			gui_grid_size = undoredo_guisize_tracker.get_current();

			wificuda::gui_resize(antennatable.size(), gui_grid_size);
			wificuda::graphical_recalc_polar(station_id, location, gui_grid_size, antenna.settings());
			wificuda::graphical_init(station_id, gui_grid_size, antennatable);
			wificuda::graphical_update(station_id, antenna.settings());

			return;
		}

		if (new_dimension.is_zero())
		{
			spdlog::error("New size is not valid");
			throw std::runtime_error("Enter valid size");
		}
		else // grid_size = new_dimension
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

		undoredo_guisize_tracker.emplace({ rows, cols });
		const auto& gui_grid_size = undoredo_guisize_tracker.get_current();

		wificuda::gui_resize(antennatable.size(), gui_grid_size);
		wificuda::graphical_recalc_polar(station_id, location, gui_grid_size, antenna.settings());
		wificuda::graphical_init(station_id, gui_grid_size, antennatable);
	}


	void undo()
	{
		if (undoredo_settings_tracker.has_changed())
		{
			Settings& current = undoredo_settings_tracker.get_current();
			undoredo_settings_tracker.undo();

			Settings& previous = undoredo_settings_tracker.get_current();
			update_changes(current, previous);
		}

		if (undoredo_guisize_tracker.has_changed())
		{
			Dimensions<unsigned>& current = undoredo_guisize_tracker.get_current();
			undoredo_guisize_tracker.undo();

			Dimensions<unsigned>& previous = undoredo_guisize_tracker.get_current();
			wificuda::gui_resize(antennatable.size(), previous);
			wificuda::graphical_recalc_polar(station_id, location, previous, antenna.settings());
			wificuda::graphical_init(station_id, previous, antennatable);
			wificuda::graphical_update(station_id, antenna.settings());
		}
	}

	void redo()
	{
		if (undoredo_settings_tracker.has_changed())
		{
			Settings& current = undoredo_settings_tracker.get_current();
			undoredo_settings_tracker.redo();

			Settings& previous = undoredo_settings_tracker.get_current();
			update_changes(current, previous);
		}

		if (undoredo_guisize_tracker.has_changed())
		{
			Dimensions<unsigned>& current = undoredo_guisize_tracker.get_current();
			undoredo_guisize_tracker.redo();

			Dimensions<unsigned>& previous = undoredo_guisize_tracker.get_current();
			wificuda::gui_resize(antennatable.size(), previous);
			wificuda::graphical_recalc_polar(station_id, location, previous, antenna.settings());
			wificuda::graphical_init(station_id, previous, antennatable);
			wificuda::graphical_update(station_id, antenna.settings());
		}

	}

	/* when user resets all the changes in the simulation */
	void init_sim()
	{
		antenna.reset();
		location = init_location;

		wificuda::numerical_recalc_polar(station_id, antennatable.size(), location);
		wificuda::numerical_init(station_id, ms_stations_num, antennatable);
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
		antenna(antenna_table[id]),
		init_location(antenna.settings().location),
		ms_station_loc(ms_pos_list),
		ms_stations_num(ms_station_loc.size()),
		host_simdata_coefficients_ptr(nullptr),
		host_guidata_coefficients_ptr(nullptr)
	{
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
