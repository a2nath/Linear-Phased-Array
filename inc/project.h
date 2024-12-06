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

struct GraphicsHelper
{
	Logger& logger;
	unsigned num_tx, num_rx;
	size_t rows;
	size_t cols;
	graphics::DataSync       synced_state;
	std::vector<double_v>    raw_cow_data; // raw split data that shows signal strength of each COW
	double_v&    raw_mrg_data; // merged data that shows SINR based on strongest signal
	std::vector<size_t>  cow_sigids;   // max signal ids for each pixel, to be used with merged
	std::vector<size_t>      rx_ids_p_idx; // rx index to find the SINR value from the mrg lut data

	bool is_rendering;
	double noise_factor;

	void init(const double& inoise_factor)
	{
		raw_cow_data.assign(raw_cow_data.size(), double_v(rows * cols)); // resize the part with individual cow heat data (1.)
		noise_factor = inoise_factor;
	}

	void init_late()
	{
		raw_mrg_data.resize(rows * cols); // resize the part with merged cow heat data (2.)
		cow_sigids.resize(rows * cols);   // cow id to use as signal source (3.) - this can be a local variable
		rx_ids_p_idx.resize(num_rx);      // indices inside the merged cow heat data to represent receiver SINR. See (2.) (4.)
	}

	/* debug with logger */
	void setup_cow_heat(unsigned cow_idx)
	{
		size_t index = 0;

		for (size_t row = 0; row < rows; ++row)
		{
			//logger.write("cow " + str(cow_idx) + '\n');

			for (size_t col = 0; col < cols; ++col)
			{
				double num = lin2dB(raw_cow_data[cow_idx][index]); // already filled with heat from "setup_tx" now convert
				raw_cow_data[cow_idx][index] = num;
				++index;
			}
		}
	}

	/* GUI setup for simulation changes */
	void update_tx(Cow& cow, const double& power, const double& antenna_dir, const double& scan_angle, const size_t& new_rows, const size_t& new_cols, const Placements& placement)
	{
		if (rows != new_rows || rows != new_cols)
		{
			rows = new_rows;
			cols = new_cols;
			raw_cow_data.resize(raw_cow_data.size(), double_v(rows * cols));
			init_late();
		}

		spdlog::info("Reconfigering #" + str(cow.sid()) + " transmitter in visuals: \
			power(" + str(power) + \
			") antenna_dir (" + str(antenna_dir) + \
			") scan_angle(" + str(scan_angle) + \
			") rows(" + str(rows) + \
			") cols(" + str(cols));
			//") placement(" + str(placement) + ")");


		cow.reset_gui(rows, cols, antenna_dir, placement);
		cow.gui_udpate(scan_angle);
		cow.set_power(power);

		/* update cow heat data */
		cow.heatmap(raw_cow_data[cow.sid()]);

		setup_cow_heat(cow.sid()); // <-- this is need to convert to logarithmic, or else all graphs will be plain yellow
	}

	/* GUI setup for all cows together, inputs: rows, cols */
	void setup_tx(cow_v& cows, const double_v& antenna_power_list, const double_v& scan_alpha_list, const size_t& rows, const size_t& cols)
	{
		spdlog::info("Setting up " + str(cows.size()) + " transmitters in visuals");

		for (unsigned c = 0; c < cows.size(); ++c)
		{
			cows[c].init_gui(rows, cols);
			cows[c].gui_udpate(scan_alpha_list[c]);
			cows[c].set_power(antenna_power_list[c]);

			/* update cow heat data */
			cows[c].heatmap(raw_cow_data[c]);

			setup_cow_heat(c); // <-- this is need to convert to logarithmic, or else all graphs will be plain yellow
		}
	}

	/* returns min and max sinr for transmitter [tx_id] */
	inline std::pair<double, double> compute_colorspan(unsigned tx_id)
	{
		auto& cow_raw_data = raw_cow_data[tx_id];
		auto [fmin_ite, fmax_ite] = std::minmax_element(cow_raw_data.begin(), cow_raw_data.end());

		return { *fmin_ite, *fmax_ite };
	}

	/* returns relative min and max sinr between transmitters */
	std::pair<double, double> compute_colorspan(cow_v& txlist)
	{
		auto double_min = std::numeric_limits<double>::max(), double_max = std::numeric_limits<double>::lowest();
		for (auto& cow : txlist)
		{
			auto [fmin, fmax] = compute_colorspan(cow.sid());

			double_min = std::min(fmin, double_min);
			double_max = std::max(fmax, double_max);
		}

		return { double_min, double_max };

	}

	void render(cow_v& txlist,
		const placement_v& mobile_stations_loc,
		const placement_v& base_stations_loc,
		const double_v& bs_txpower,
		const double_v& bs_theta_c,
		const double_v& scan_alpha_list)
	{
		if (raw_cow_data[0].empty())
		{
			spdlog::critical("Visualization not enabled, exiting render() function");
			return;
		}

		init_late();

		setup_tx(txlist, bs_txpower, scan_alpha_list, rows, cols);

		std::pair<double, double> sep_channels = compute_colorspan(txlist);


		/*auto& signal_table = raw_cow_data;
		auto& strongest_sigid = cow_sigids;*/


		size_t index = 0;

		for (size_t row = 0; row < rows; ++row)
		{
			for (size_t col = 0; col < cols; ++col)
			{
				/* find max and the corresponding cow-id (that becomes SIGNAL) */
				double max_value = sep_channels.first;
				int cow_idx = -1;

				for (auto& cow : txlist)
				{
					if (max_value < raw_cow_data[cow.sid()][index])
					{
						max_value = raw_cow_data[cow.sid()][index];
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


		auto& signal_output = raw_mrg_data;
		auto& indices_ms = rx_ids_p_idx;

		std::unordered_map<size_t, std::unordered_set<size_t>> rx_coords_hash_lut;

		for (auto& loc : mobile_stations_loc)
		{
			rx_coords_hash_lut[loc.x].emplace(loc.y);
		}

		for (size_t row = 0; row < rows; ++row)
		{
			for (size_t col = 0; col < cols; ++col)
			{
				const unsigned& cow_idx = cow_sigids[px_index];
				const double& signal = raw_cow_data[cow_idx][px_index];

				double interference = 0;

				for (unsigned c = (cow_idx + 1) % txlist.size(); c != cow_idx;)
				{
					interference += raw_cow_data[c][px_index];
					c = (c + 1) % txlist.size();
				}

				num = lin2dB(signal / (interference + noise_factor));
				signal_output[px_index] = num;

				// it is not one of the RX station coordinates
				if (!(rx_coords_hash_lut.find(row) == rx_coords_hash_lut.end() || rx_coords_hash_lut[row].find(col) == rx_coords_hash_lut[row].end()))
				{
					indices_ms[rx_index++] = px_index;
				}

				++px_index;
			}
		}

		std::pair<double, double> com_channels = compute_colorspan(raw_cow_data.size() - 1);

		graphics::render(logger,
			mobile_stations_loc,
			base_stations_loc,
			raw_cow_data,
			cow_sigids,
			bs_txpower,
			bs_theta_c,
			scan_alpha_list,
			rows,
			cols,
			double_v{ sep_channels.first, com_channels.first },
			double_v{ sep_channels.second, com_channels.second },
			synced_state,
			is_rendering);
	}

	void plot(cow_v& txlist,
		const placement_v& mobile_stations_loc,
		const placement_v& base_stations_loc,
		const double_v& bs_txpower,
		const double_v& bs_theta_c,
		const double_v& scan_alpha_list)
	{
		if (raw_cow_data[0].empty())
		{
			spdlog::critical("Visualization not enabled, exiting plot() function");
			return;
		}

		if (txlist[0].gui_state() == false)
			setup_tx(txlist, bs_txpower, scan_alpha_list, rows, cols);

		// else capture plot as it is

		std::vector<unsigned> tx_ids(txlist.size() + 1);
		std::iota(tx_ids.begin(), tx_ids.end(), 0);

		auto sep_ch = compute_colorspan(txlist);
		auto com_ch = compute_colorspan(txlist.size() - 1);

		for (auto& id : tx_ids)
		{

			std::string num = id == tx_ids.back() ? "combined" : str(id);

			graphics::capture_plot(logger,
				"transmitter_" + num + ".png",
				mobile_stations_loc,
				base_stations_loc,
				raw_cow_data[id],
				bs_txpower,
				bs_theta_c,
				scan_alpha_list,
				rows,
				cols,
				id == tx_ids.back() ? com_ch.first : sep_ch.first,
				id == tx_ids.back() ? com_ch.second : sep_ch.second);
		}
	}

	GraphicsHelper(const size_t num_transmitters, const size_t num_receivers, const size_t& pixel_rows, const size_t& pixel_cols, Logger& ilogger)
		:
		num_tx(num_transmitters),
		num_rx(num_receivers),
		rows(pixel_rows),
		cols(pixel_cols),
		synced_state(num_transmitters),
		is_rendering(true),
		raw_cow_data(num_tx + 1),
		raw_mrg_data(raw_cow_data.back()),
		logger(ilogger)
	{
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
			cows[c].antenna_update(scan_angles[c]);
			cows[c].set_power(power_nums[c]);
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

	const unsigned& timeslot;
	const double&   frequency;
	const double    lambda;
	const double&   bandwidth;
	const double&   symrate;
	const double&   blockspersym;
	const double&   antenna_height;
	const unsigned& mobile_station_count;
	const unsigned& base_station_count;
	const unsigned& timeslot_count;
	const double    sinr_limit_linear;
	const double_v bs_theta_c;
	const placement_v& base_stations_loc;
	const placement_v& mobile_stations_loc;
	const unsigned_v& bs_antenna_counts;
	const double_v& power_range_dBm;
	const double_v& scan_angle_range;
	const double_v& antenna_spacing;
	const double_v& antenna_dims;

	GraphicsHelper visuals;

	const std::vector<double_v>& bs_tx_requested_power_watts;
	const std::vector<double_v>& bs_requested_scan_alpha_rad;
	const std::vector<std::vector<unsigned>>& ms2bs_requested_bindings;

	void setup(const std::vector<double>& grxlist, const double& system_noise_lin)
	{
		for (unsigned i = 0; i < mobile_station_count; ++i)
		{
			stations.emplace_back(i, system_noise_lin, grxlist[i]);
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
		double_v antenna_power, scan_angles;

		simhelper->get_lut_tx_power(antenna_power);
		simhelper->get_lut_scan_angle(scan_angles);
		visuals.render(cows, mobile_stations_loc, base_stations_loc, antenna_power, bs_theta_c, scan_angles);
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
		visuals(args.base_station_count, args.mobile_station_count, args.field_size[0], args.field_size[1], ilogger),
		bs_tx_requested_power_watts(args.tx_powerlist().data),
		bs_requested_scan_alpha_rad(args.tx_alphalist().data),
		ms2bs_requested_bindings(args.ms_id_selections.binding_data)
	{
		setup(double_v(args.mobile_station_count, cached::log2lin(args.gain_gtrx)),
			cached::log2lin(getThermalSystemNoise(bandwidth, args.system_noise)));

		if (!args.nogui)
		{
			/* initialize the visual helper */
			spdlog::info("Visualization enabled");
			visuals.init(stations[0].get_nf());
		}
	}
};


class Simulator
{
protected:
    Logger& logger;
    std::string sim_error;
    std::vector<Cow> cows;
    std::vector<Station> stations;
    PerfMon* monitor;

    //std::vector<double> dBm2watts;              // in watts
    //std::vector<double> deg2rads;  // scan angle
    const double&   frequency;
    const double    lambda;
    const double&   bandwidth;
    const double&   symrate;
    const double&   blockspersym;
    const double&   antenna_height;
    const double&   ms_gain_gtrx_lin;
    const double&   system_noise_lin;
    const unsigned& mobile_station_count;
    const unsigned& base_station_count;
    const unsigned& timeslot_count;
    const double&   sinr_limit_linear;
    const std::vector<double>& bs_theta_c;
    const std::vector<Placements>& base_stations_loc;
    const std::vector<Placements>& mobile_stations_loc;
    const std::vector<unsigned>& bs_antenna_counts;
    const std::vector<double>& power_range_dBm;
    const std::vector<double>& scan_angle_range;
    const std::vector<double>& antenna_spacing;
    const std::vector<double>& antenna_dims;
    const bool&     showgui;
    const bool&     debug;

    const std::vector<std::vector<double>>& bs_tx_requested_power_dBm;
    const std::vector<std::vector<double>>& bs_requested_scan_alpha_deg;
    const std::vector<std::vector<unsigned>>& ms2bs_requested_bindings;

    void setup()
    {
        //for (unsigned power = power_range_dBm.front(); power <= power_range_dBm.back(); ++power)
        //{
        //    dBm2watts.emplace_back(cached::dBm2watt(power));
        //}
        //
        //for (unsigned alpha = power_range_dBm.front(); alpha <= power_range_dBm.back(); ++alpha)
        //{
        //    deg2rads.emplace_back(cached::dBm2watt(alpha));
        //}

        //for (auto& theta : bs_theta_c)
        //{
        //    theta = deg2rad(theta);
        //}

        for (unsigned i = 0; i < mobile_stations_loc.size(); ++i)
        {
            stations.emplace_back(i, cows, system_noise_lin);
        }

        /* setup the system */
        for (unsigned bs_id = 0; bs_id < base_station_count; ++bs_id)
        {
            cows.emplace_back(bs_id,
                base_stations_loc[bs_id],
                mobile_stations_loc,
                ms_gain_gtrx_lin,
                bs_antenna_counts[bs_id],
                lambda,
                antenna_spacing[bs_id],
                bs_theta_c[bs_id],
                antennadim(antenna_dims[0], antenna_dims[1])
                );
        }
    }

    void run_single_timelot()
    {
        auto& simhelper = monitor->params;

        /* set the selected stations as per bindings */
        std::vector<unsigned> select_stations;
        simhelper.get_bindings(select_stations);

        /* create a new list of power nums and run calculations on its permutations */
        std::vector<double> power_list;
        simhelper.get_power(power_list);

        /* update the SINR for each scenario (base-station-count factorial permutations) */
        bool permutation_state = true;
        while (permutation_state)
        {
            /* change the antenna power across all base stations */
            simhelper.set_power(power_list);

            for (size_t bs_id = 0; bs_id < base_station_count; ++bs_id)
            {
                auto& station = stations[select_stations[bs_id]];
                station.set_rx(bs_id);

                monitor->update(bs_id, station.get_sinr(), select_stations);
            }

            permutation_state = simhelper.get_perm_power(power_list);
        }
    }

    void run_multi_timeslot()
    {
        auto& simhelper = monitor->params;

        /* setup the cows first such that the simulation is set with alphas and TX powers */
        for (unsigned slot = 0; slot < simhelper.timeslots; ++slot)
        {
            run_single_timelot();
            simhelper.inc_timeslot(); // increment the timeslot to get the next parameters
        }
    }

public:

    /* run the simulation */
    void run()
    {
        SimulationHelper simparams(cows,
            timeslot_count,
            mobile_station_count,
            bs_tx_requested_power_dBm,
            bs_requested_scan_alpha_deg,
            ms2bs_requested_bindings
        );

        monitor = new PerfMon(simparams);

        if (timeslot_count > 1)
            run_multi_timeslot();
        else
            run_single_timelot();
    }


    void print()
    {
        monitor->print();
    }

    ~Simulator()
    {
        delete monitor;
    }

    Simulator(const MyArgs& args, Logger& ilogger)
        :
        logger(ilogger),
        sim_error(""),
        monitor(nullptr),
        frequency(args.frequency),
        lambda(getLambda(frequency)),
        bandwidth(args.bandwidth),
        symrate(args.symrate),
        blockspersym(args.blockspersym),
        antenna_height(args.antenna_height),
        ms_gain_gtrx_lin(cached::log2lin(args.gain_gtrx)),
        system_noise_lin(cached::log2lin(getThermalSystemNoise(bandwidth, args.system_noise))),
        mobile_station_count(args.mobile_station_count),
        base_station_count(args.base_station_count),
        timeslot_count(args.timeslots),
        sinr_limit_linear(cached::log2lin(args.sinr_limit_dB)),
        bs_theta_c(cached::deg2rad(args.bs_theta_c)),
        base_stations_loc(args.base_stations_loc.data),
        mobile_stations_loc(args.mobile_stations_loc.data),
        bs_antenna_counts(args.bs_antenna_count),
        power_range_dBm(args.power_range_dBm),
        scan_angle_range(args.scan_angle_range),
        antenna_spacing(args.antenna_spacing),
        antenna_dims(args.antenna_dims),
        showgui(args.showgui),
        debug(args.debug),
        bs_tx_requested_power_dBm(args.bs_tx_power_dBm.value().data),
        bs_requested_scan_alpha_deg(args.bs_scan_alpha_deg.value().data),
        ms2bs_requested_bindings(args.ms_id_selections.binding_data)
    {
        setup();
    }
};
