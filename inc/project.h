#pragma once
#include <unordered_set>
#include <vector>
#include <cassert>
#include <any>
#include <ostream>
#include <limits>
#include <queue>
#include "common.h"
#include "network.h"
#include "station.h"
#include "random.h"
#include "args.h"

#ifdef GRAPHICS
#include "visuals.h"

#include <future>
#include <thread>

struct GraphicsHelper
{
	Logger& logger;
	unsigned num_tx, num_rx;
	unsigned known_height;
	unsigned known_width;
	double noise_factor;
	std::vector<double_v> raw_dbg_lin_data; // raw split data that shows signal strength of each COW
	std::vector<double_v> raw_dbg_dBm_data; // merged data that shows SINR based on strongest signal
	std::vector<double_v> raw_mrg_data; // merged data that shows SINR based on strongest signal
	std::vector<size_t>   cow_sigids;   // max signal ids for each pixel, to be used with merged
	std::vector<size_t>   rx_ids_p_idx; // rx index to find the SINR value from the mrg lut data
	state_v curr_states;
	Bench btime;

	bool verbose;
	bool is_rendering;

	inline bool ready() const
	{
		return raw_dbg_lin_data.size() > 0;
	}

	//void rqst_secondary_mem()
	//{
		 // resize the part with merged cow heat data (2.)
		//cow_sigids.resize(rows * cols);   // cow id to use as signal source (3.) - this can be a local variable
		//rx_ids_p_idx.resize(num_rx);      // indices inside the merged cow heat data to represent receiver SINR. See (2.) (4.)
	//}

	/* GUI setup for simulation changes */
	void update_tx(Cow& cow, const graphics::State& state)
	{
		spdlog::info("Reconfigering #" + str(cow.sid()) + " transmitter in visuals:\n \
			power(" + str(state.settings.power, 2) + \
			") antenna_dir(" + str(state.settings.theta_c, 2) + \
			") scan_angle(" + str(state.settings.alpha, 2) + \
			") rows(" + str(known_height) + \
			") cols(" + str(known_width) + ")");
			//") placement(" + str(placement) + ")");

		cow.update(state.settings, state.location);

		/* update cow heat data */
		cow.heatmap(raw_dbg_lin_data[cow.sid()]);
	}

	/* GUI setup for all cows together, inputs: rows, cols */
	void setup_tx(cow_v& cows, const double_v& lut_power_list, const double_v& lut_scan_angle_list)
	{
		raw_dbg_lin_data.assign(num_tx, double_v(known_width * known_height));
		spdlog::info("Setting up " + str(cows.size()) + " transmitters in visuals");

		for (auto& cow : cows)
		{
			spdlog::info("Setting up TX " + str(cow.sid()) + " for GUI");

			cow.init_gui(known_width, known_height);
			cow.update_minimal(lut_power_list[cow.sid()], lut_scan_angle_list[cow.sid()]);
			cow.heatmap(raw_dbg_lin_data[cow.sid()]);
		}
	}

	void generate_data_debug(cow_v& txlist)
	{
		spdlog::info("Translate linear to dBm values in debug");

		for (int i = 0; i < txlist.size(); ++i)
		{
			size_t index = 0;

			for (size_t row = 0; row < known_height; ++row)
			{
				//logger.write("cow " + str(cow_idx) + '\n');

				for (size_t col = 0; col < known_width; ++col)
				{
					raw_dbg_dBm_data[i][index] = cached::watt2dBm(raw_dbg_lin_data[i][index]); // already filled with heat from "setup_tx" now convert
					++index;
				}

			}
		}
	};

	void generate_data_intf(const cow_v& txlist)
	{
		spdlog::info("Create interference data");
		raw_mrg_data.assign(num_tx, double_v(known_height * known_width));

		size_t pxl_idx = 0;
		double num = 0;

		for (unsigned row = 0; row < known_height; ++row)
		{
			for (unsigned col = 0; col < known_width; ++col)
			{
				for (auto& cow : txlist)
				{
					auto& cow_idx = cow.sid();

					const double& signal = raw_dbg_lin_data[cow_idx][pxl_idx];
					double interference = 0;

					for (unsigned c = (cow_idx + 1) % txlist.size(); c != cow_idx;)
					{
						interference += raw_dbg_lin_data[c][pxl_idx];
						c = (c + 1) % txlist.size();
					}

					raw_mrg_data[cow_idx][pxl_idx] = lin2dB(signal / (interference + noise_factor));
				}

				++pxl_idx;
			}
		}
	}

	inline std::future<void> async_generate_data_intf(const cow_v& txlist)
	{
		spdlog::info("- Async - ");
		return std::async(std::launch::async, [&]()
			{
				generate_data_intf(txlist);
			});

		// Optionally wait for the future if synchronization is needed later
		// future.wait();
	}

	std::thread gui_change_detect(cow_v& txlist, graphics::DataSync& sync)
	{
		return std::thread([&]()
			{
				while (is_rendering)
				{
					std::unique_lock<std::mutex> lock_compute(graphics::compute_sim_mutex);  // Lock the mutex
					graphics::consig.wait(lock_compute, [&]()
						{
							return sync.is_computing > 0 || sync.resize_event || !is_rendering;
						}
					);  // Wait for new data or rendering to stop

					if (!is_rendering)
					{
						break;
					}

					if (verbose)
					{
						btime.bench_start();
					}

					if (sync.is_debugging && sync.mainq.empty())
					{
						generate_data_debug(txlist);

						if (verbose)
						{
							btime.mark();
							spdlog::info("Debug data rebuilt in " + str(btime.get()) + " milliseconds");
						}
					}

					while (sync.mainq.size())
					{
						auto& state = *sync.mainq.front();

						update_tx(txlist[state.tx_idx], state);

						if (sync.is_debugging)
						{
							generate_data_debug(txlist);
							async_generate_data_intf(txlist); // no rush when viewing in debug mode
						}
						else
						{
							generate_data_intf(txlist);
						}

						sync.pop();

						if (verbose)
						{
							btime.mark();
							spdlog::info("Heat Checker Thread: Update Done on TX-" + str(state.tx_idx)\
								+ " in " + str(btime.get()) + " milliseconds");
						}
					}

					if (verbose)
					{
						btime.bench_stop();
						spdlog::info("Total time to finish thread: " + str(btime.get()) + " milliseconds");
					}

					sync.is_computing = 0;
					sync.render_tx_id = sync.def_render_tx_id;

					graphics::consig.notify_one();
				}
			});
	}

	void debug_mode_max_signal_map(cow_v& txlist,
		const placement_v& mobile_stations_loc)
	{
		//if (raw_dbg_lin_data[0].empty())
		//{
		//	spdlog::critical("Visualization not enabled, exiting render() function");
		//	return;
		//}
		//
		//
		//rqst_secondary_mem();
		//
		//setup_tx(txlist, bs_txpower, scan_alpha_list, rows, cols);


		size_t index = 0;

		for (size_t row = 0; row < known_height; ++row)
		{
			for (size_t col = 0; col < known_width; ++col)
			{
				/* find max and the corresponding cow-id (that becomes SIGNAL) */
				double max_value = 0;// sep_channels.first;
				int cow_idx = -1;

				for (auto& cow : txlist)
				{
					if (max_value < raw_dbg_lin_data[cow.sid()][index])
					{
						max_value = raw_dbg_lin_data[cow.sid()][index];
						cow_idx = cow.sid();
					}
				}

				cow_sigids[index] = cow_idx;
				++index;

			}
		}

		size_t px_index = 0;
		size_t rx_index = 0;
		double num = 0;


		auto& indices_ms = rx_ids_p_idx;

		std::unordered_map<size_t, std::unordered_set<size_t>> rx_coords_hash_lut;

		for (auto& loc : mobile_stations_loc)
		{
			rx_coords_hash_lut[loc.x].emplace(loc.y);
		}

		for (size_t row = 0; row < known_height; ++row)
		{
			for (size_t col = 0; col < known_width; ++col)
			{
				const unsigned& cow_idx = cow_sigids[px_index];
				const double& signal = raw_dbg_lin_data[cow_idx][px_index];

				double interference = 0;

				for (unsigned c = (cow_idx + 1) % txlist.size(); c != cow_idx;)
				{
					interference += raw_dbg_lin_data[c][px_index];
					c = (c + 1) % txlist.size();
				}

				num = lin2dB(signal / (interference + noise_factor));
				raw_mrg_data[cow_idx][px_index] = num;

				// it is not one of the RX station coordinates
				if (!(rx_coords_hash_lut.find(row) == rx_coords_hash_lut.end() || rx_coords_hash_lut[row].find(col) == rx_coords_hash_lut[row].end()))
				{
					indices_ms[rx_index++] = px_index;
				}

				++px_index;
			}
		}

	}

	void render(cow_v& txlist,
		graphics::DataSync& sync,
		const placement_v& mobile_stations_loc,
		const double_v& lut_power_list,
		const double_v& lut_scan_angle_list)
	{
		state_v init_states;

		if (!ready())
		{
			/* fill the TX indivisual simulation data across the entire grid */
			setup_tx(txlist, lut_power_list, lut_scan_angle_list);

			/* fill the overall interference simulation across the entire grid */
			generate_data_intf(txlist);

			for (auto& tx : txlist)
			{
				init_states.emplace_back(tx.get_state());
			}
		}

		spdlog::info("Beginning Visualization from Initial TX states");

		graphics::render(logger,
			mobile_stations_loc,
			init_states,
			curr_states,
			raw_dbg_lin_data,
			raw_dbg_dBm_data,
			raw_mrg_data,
			known_height,
			known_width,
			sync,
			is_rendering);
	}

	void plot(cow_v& txlist,
		const placement_v& mobile_stations_loc,
		const placement_v& base_stations_loc,
		const double_v& bs_txpower,
		const double_v& bs_theta_c,
		const double_v& scan_alpha_list)
	{
		if (!ready())
		{
			setup_tx(txlist, bs_txpower, scan_alpha_list);
			generate_data_intf(txlist);

			curr_states.clear();
			for (auto& tx : txlist)
			{
				curr_states.emplace_back(tx.get_state());
			}
		}
		// else capture plot as it is

		spdlog::info("Capturing plots\n\tDebug Plots");

		graphics::capture_plot(logger,
			"transmitter_dbg_",
			curr_states,
			mobile_stations_loc,
			raw_dbg_lin_data,
			known_height,
			known_width,
			true);

		spdlog::info("\tInterference Plots");

		graphics::capture_plot(logger,
			"transmitter_int_",
			curr_states,
			mobile_stations_loc,
			raw_mrg_data,
			known_height,
			known_width);
	}

	GraphicsHelper(Logger& ilogger,
		const size_t num_transmitters,
		const size_t num_receivers,
		const unsigned& pixel_rows,
		const unsigned& pixel_cols,
		const double& inoise_factor,
		const bool& iverbose = true)
		:
		logger(ilogger),
		num_tx(num_transmitters),
		num_rx(num_receivers),
		known_height(pixel_rows),
		known_width(pixel_cols),
		noise_factor(inoise_factor),
		raw_dbg_lin_data(num_tx),
		raw_mrg_data(num_tx),
		is_rendering(true),
		verbose(iverbose)
	{
		raw_dbg_lin_data.clear();
		raw_dbg_dBm_data.clear();
		raw_mrg_data.clear();
	}
};
#endif


struct SimulationHelper
{
	cow_v& cows;
	const unsigned& cow_count, ms_stations;
	const std::vector<double_v>& powers_lut;                 // TX power level from base station in integer dBm lut
	const std::vector<double_v>& alphas_lut;                 // Antenna array directions in integer rads lut
	const std::vector<std::vector<unsigned>>& binding_station_ids_lut;  // base_station to station id binding lut
	std::priority_queue<dataitem_t, std::vector<dataitem_t>, data_comparator> pqueue;

	unsigned timeslot_idx;   // timeslot in the schedule
	unsigned timeslots;      // total timeslots
	unsigned nodes_rem;

	void inc_timeslot()
	{
		++timeslot_idx;
	}

	const unsigned& get_timeslots() const
	{
		return timeslots;
	}

	/* gets called to get BS to MS bindings for a single timeslot. Output size: #bs */
	inline void get_bindings(std::vector<unsigned>& output) const
	{
		output = binding_station_ids_lut[timeslot_idx % timeslots];
	}

	/* get antenna array power level in each timeslot */
	inline void get_lut_tx_power(double_v& output) const
	{
		output = powers_lut[timeslot_idx % timeslots];
	}

	/* get permutations of some simulation property */
	template<class T>
	inline bool get_perm(std::vector<T>& list) const
	{
		return std::next_permutation(list.begin(), list.end());
	}

	/* set antenna array directivity before starting each simulation */
	inline void get_lut_scan_angle(double_v& output) const
	{
		output = alphas_lut[timeslot_idx % timeslots];
	}

	/* output into console */
	void printout()
	{
		while (pqueue.size())
		{
			std::cout << pqueue.top() << std::endl;
			pqueue.pop();
		}
	}

	/* need to update all coefficients before getting rx power in any 1 station */
	void setup_tx()
	{
		spdlog::info("Initializating " + str(cows.size()) + " transmitters for simulation");
		double_v scan_angles, power_nums;
		get_lut_scan_angle(scan_angles);
		get_lut_tx_power(power_nums);

		for (unsigned c = 0; c < cows.size(); ++c)
		{
			cows[c].update_minimal(power_nums[c], scan_angles[c]);
		}
	}

	/* performance logging and analysis latter */
	void update_results(const unsigned& bs_id, const unsigned& ms_id, const double& sinr, const Placements& mx_location)
	{
		pqueue.emplace(
			timeslot_idx,
			cows[bs_id].getxPower(),
			cows[bs_id].alpha(),
			mx_location,
			cows[bs_id].sid(),
			ms_id,
			sinr);
	}

	SimulationHelper(cow_v& cowlist,
		const unsigned& timeslot_num,
		const unsigned& timeslot_count,
		const unsigned& mobile_stations,
		const std::vector<double_v>& power_bindings,
		const std::vector<double_v>& alpha_bindings,
		const std::vector<std::vector<unsigned>>& station_bindings)
		:
		cows(cowlist),
		cow_count(cowlist.size()),
		ms_stations(mobile_stations),
		binding_station_ids_lut(station_bindings),
		powers_lut(power_bindings),
		alphas_lut(alpha_bindings),
		timeslot_idx(timeslot_num),
		timeslots(timeslot_count),
		nodes_rem(mobile_stations)
	{
	}
};

class Simulator
{
	Logger& logger;
	std::string sim_error;
	cow_v cows;
	sta_v stations;
	SimulationHelper* simhelper;

	const unsigned&    timeslot;
	const double&      frequency;
	const double       lambda;
	const double&      bandwidth;
	const double&      symrate;
	const double&      blockspersym;
	const double&      antenna_height;
	const unsigned&    mobile_station_count;
	const unsigned&    base_station_count;
	const unsigned&    timeslot_count;
	const double       sinr_limit_linear;
	const double_v     bs_theta_c;
	const placement_v& base_stations_loc;
	const placement_v& mobile_stations_loc;
	const unsigned_v&  bs_antenna_counts;
	const double_v&    power_range_dBm;
	const double_v&    scan_angle_range;
	const double_v&    antenna_spacing;
	const double_v&    antenna_dims;
	const double       system_noise_factor_w;
	GraphicsHelper     visuals;

	const std::vector<double_v>& bs_tx_requested_power_watts;
	const std::vector<double_v>& bs_requested_scan_alpha_rad;
	const std::vector<std::vector<unsigned>>& ms2bs_requested_bindings;

	void setup(const std::vector<double>& grxlist)
	{
		for (unsigned i = 0; i < mobile_station_count; ++i)
		{
			stations.emplace_back(i, system_noise_factor_w, grxlist[i]);
		}

		/* setup the system */
		for (unsigned bs_id = 0; bs_id < base_station_count; ++bs_id)
		{
			cows.emplace_back(bs_id,
				base_stations_loc[bs_id],
				mobile_stations_loc,
				//ms_gain_gtrx_lin,
				bs_antenna_counts[bs_id],
				lambda,
				antenna_spacing[bs_id],
				bs_theta_c[bs_id],
				antennadim(antenna_dims[0], antenna_dims[1])
				);
		}

		/* initialize the simulation helper */
		simhelper = new SimulationHelper(
			cows,
			timeslot,
			timeslot_count,
			mobile_station_count,
			bs_tx_requested_power_watts,
			bs_requested_scan_alpha_rad,
			ms2bs_requested_bindings
		);
	}

	/* calcalate the rx signal based on all base stations */
	double calculate(Station& station, const unsigned& associated_bs_id)
	{
		double signal = 0, interference = 0, power;

		for (auto& cow : cows)
		{
			cow.signal_power(station.sid(), power);

			if (cow.sid() != associated_bs_id)
			{
				interference += power * station.get_grx();
			}
			else
			{
				signal = power * station.get_grx();
			}
		}

		return signal / (interference + station.get_nf());;
	}

	void run_single_timelot()
	{
		/* IMPORTANT need to set this before iterating over cows again */
		simhelper->setup_tx();

		/* set the selected stations as per bindings */
		std::vector<unsigned> select_stations;
		simhelper->get_bindings(select_stations);

		/* update the SINR for each scenario (base-station-count factorial permutations) */
		bool permutation_state_bindings = true;
		while (permutation_state_bindings)
		{
			/* create a new list of power nums and run calculations on its permutations */
			std::vector<double> power_list;
			simhelper->get_lut_tx_power(power_list);

			bool permutation_state_power = true;
			while (permutation_state_power)
			{
				for (auto& cow : cows)
				{
					auto& bs_id = cow.sid();
					auto& rx_idx = select_stations[bs_id];
					auto& rx_station = stations[rx_idx];
					auto& rx_loc = mobile_stations_loc[rx_idx];

					auto sinr = calculate(rx_station, bs_id); // reassociate with a new tx if possible
					simhelper->update_results(bs_id, rx_idx, sinr, mobile_stations_loc[rx_idx]);
				}

				permutation_state_power = simhelper->get_perm(power_list);
			}

			permutation_state_bindings = simhelper->get_perm(select_stations);
		}
	}

public:

	/* run the simulation */
	void run()
	{
		/* setup the cows first such that the simulation is set with alphas and TX powers */
		for (unsigned slot = 0; slot < timeslot_count; ++slot)
		{
			run_single_timelot();
			simhelper->inc_timeslot(); // increment the timeslot to get the next parameters
		}
	}

	void print()
	{
		simhelper->printout();
	}

	/* render a scene in a new window and show the simulation */
	void gui_run()
	{
#ifdef GRAPHICS
		double_v init_ant_power, init_scan_angle;

		simhelper->get_lut_tx_power(init_ant_power);
		simhelper->get_lut_scan_angle(init_scan_angle);
		graphics::DataSync sync(cows.size());

		auto queue_thread = visuals.gui_change_detect(cows, sync);

		visuals.render(cows, sync, mobile_stations_loc, init_ant_power, init_scan_angle);

		queue_thread.join();
#endif
	}

	/* drop of a plot in png format in the output directory */
	void gui_print()
	{
#ifdef GRAPHICS
		double_v antenna_power, scan_angles;

		simhelper->get_lut_tx_power(antenna_power);
		simhelper->get_lut_scan_angle(scan_angles);

		visuals.plot(cows, mobile_stations_loc, base_stations_loc, antenna_power, bs_theta_c, scan_angles);
#endif
	}

	~Simulator()
	{
		delete simhelper;
	}

	Simulator(const MyArgs& args, Logger& ilogger)
		:
		logger(ilogger),
		sim_error(""),
		simhelper(nullptr),
		timeslot(args.timeslot),
		frequency(args.frequency),
		lambda(getLambda(frequency)),
		bandwidth(args.bandwidth),
		symrate(args.symrate),
		blockspersym(args.blockspersym),
		antenna_height(args.antenna_height),
		mobile_station_count(args.mobile_station_count),
		base_station_count(args.base_station_count),
		timeslot_count(args.timeslots),
		sinr_limit_linear(cached::log2lin(args.sinr_limit_dB)),
		bs_theta_c(cached::deg2rad(args.bs_theta_c)),
		base_stations_loc(args.tx_loc.data),
		mobile_stations_loc(args.rx_loc.data),
		bs_antenna_counts(args.bs_antenna_count),
		power_range_dBm(args.power_range_dBm),
		scan_angle_range(args.scan_angle_range),
		antenna_spacing(args.antenna_spacing),
		antenna_dims(args.antenna_dims),
		system_noise_factor_w(dBm2watt(getThermalSystemNoise(args.bandwidth, args.system_noise))),
		visuals(ilogger, args.base_station_count, args.mobile_station_count,
			args.field_size[0], args.field_size[1], system_noise_factor_w),
		bs_tx_requested_power_watts(args.tx_powerlist().data),
		bs_requested_scan_alpha_rad(args.tx_alphalist().data),
		ms2bs_requested_bindings(args.ms_id_selections.binding_data)
	{
		setup(double_v(mobile_station_count, cached::log2lin(args.gain_gtrx)));

		if (!args.nogui)
		{
			/* initialize the visual helper */
			spdlog::info("Visualization enabled for field size "\
				+ str(args.field_size[0]) + "x" + str(args.field_size[1]) + " meters");
		}
	}
};
