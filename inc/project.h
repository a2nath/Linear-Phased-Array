#pragma once
#include <unordered_set>
#include <vector>
#include <cassert>
#include <any>
#include <ostream>
#include "common.h"
#include "network.h"
#include "station.h"
#include "random.h"
#include "args.h"
#include <optional>

//using cow_id_distribution = std::uniform_int_distribution<int>;

struct SimulationHelper
{
    std::vector<Cow>& cows;
    const unsigned& cow_count, ms_stations;
    const std::vector<std::vector<double>>& powers_lut;                 // TX power level from base station in integer dBm lut
    const std::vector<std::vector<double>>& alphas_lut;                 // Antenna array directions in integer rads lut
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

    const std::vector<Cow>& cowlist() const // performance monitor needs a list
    {
        return cows;
    }

    /* gets called to get BS to MS bindings for a single timeslot. Output size: #bs */
    inline void get_bindings(std::vector<unsigned>& output) const
    {
        output = binding_station_ids_lut[timeslot_idx % timeslots];
    }

    /* set antenna array power level in each timeslot */
    inline void set_power(const std::vector<double>& power_list) const
    {
        for (unsigned i = 0; i < power_list.size(); ++i)
        {
            cows[i].set_power(power_list[i]);
        }
    }

    /* set antenna array power level in each timeslot */
    inline void reset_power() const
    {
        set_power(powers_lut[timeslot_idx % timeslots]);
    }

    /* get antenna array power level in each timeslot */
    inline void get_power(std::vector<double>& output) const
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
    inline void get_scana(std::vector<double>& output) const
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

    /* performance logging and analysis latter */
    void update_results(
        const unsigned& bs_id,
        const double& sinr,
        const std::vector<Placements>& rx_placement,
        const std::vector<unsigned>& rx_selected_idxlist)
    {
        auto& cow = cows[bs_id];

        pqueue.emplace(
            timeslot_idx,
            cow.getxPower(),
            cow.alpha(),
            rx_placement[rx_selected_idxlist[bs_id]],
            cow.sid(),
            rx_selected_idxlist[bs_id],
            sinr);
    }

    SimulationHelper(std::vector<Cow>& cowlist,
        const unsigned& timeslot_num,
        const unsigned& timeslot_count,
        const unsigned& mobile_stations,
        const std::vector<std::vector<double>>& power_bindings,
        const std::vector<std::vector<double>>& alpha_bindings,
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
protected:
    Logger& logger;
    std::string sim_error;
    std::vector<Cow> cows;
    std::vector<Station> stations;
    SimulationHelper* simhelper;

    const unsigned& timeslot;
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

    const std::vector<std::vector<double>>& bs_tx_requested_power_watts;
    const std::vector<std::vector<double>>& bs_requested_scan_alpha_rad;
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

    void run_single_timelot()
    {
        /* set the scan angles for each transmitter */
        std::vector<double> scan_angles;
        simhelper->get_scana(scan_angles);

        /* set the selected stations as per bindings */
        std::vector<unsigned> select_stations;
        simhelper->get_bindings(select_stations);

        /* update the SINR for each scenario (base-station-count factorial permutations) */
        bool permutation_state_bindings = true;
        while (permutation_state_bindings)
        {
            /* create a new list of power nums and run calculations on its permutations */
            std::vector<double> power_list;
            simhelper->get_power(power_list);

            bool permutation_state_power = true;
            while (permutation_state_power)
            {
                /* change the antenna power across all base stations */
                simhelper->set_power(power_list);

                for (size_t bs_id = 0; bs_id < base_station_count; ++bs_id)
                {
                    auto& rx_idx = select_stations[bs_id];
                    auto& rx_station = stations[rx_idx];

                    cows[bs_id].antenna_update(scan_angles[bs_id]);  // update hmatrix
                    cows[bs_id].set_power(power_list[bs_id]);        // update power

                    rx_station.set_tx_idx(bs_id); // reassociate with a new tx if possible
                    simhelper->update_results(bs_id, rx_station.get_sinr(), mobile_stations_loc, select_stations);
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
        for (unsigned slot = 0; slot < simhelper->timeslots; ++slot)
        {
            run_single_timelot();
            simhelper->inc_timeslot(); // increment the timeslot to get the next parameters
        }
    }

    void print()
    {
        simhelper->printout();
    }

    void guistat(const size_t& pixel_rows, const size_t& pixel_cols)
    {
        /* create a new list of power nums and run calculations on its permutations */
        std::vector<double> c_power_list;
        simhelper->get_power(c_power_list);

        /* set the selected stations as per bindings */
        std::vector<unsigned> select_stations;
        simhelper->get_bindings(select_stations);

        //std::unordered_map<unsigned, std::unordered_map<unsigned, bool>> visited;

        /* update the SINR for each scenario (base-station-count factorial permutations) */
        std::vector<double> power_list = c_power_list;

        for (unsigned timeslot = 0; timeslot < timeslot_count; ++timeslot)
        {
            //monitor->params.set_scan_dir();
            for (unsigned c = 0; c < simhelper->cows.size(); ++c)
            {
                auto& cow = simhelper->cows[c];
                auto& alpha = bs_requested_scan_alpha_rad[timeslot][c];

                cow.set_gui_matrix(pixel_rows, pixel_cols, cow.where());

                std::vector<double> signal_power_lin(pixel_rows * pixel_cols);
                cow.heatmap_udpate(alpha, signal_power_lin);

                size_t index = 0;
                for (size_t r = 0; r = pixel_rows; ++r)
                {
                    logger.write("cow " + str(cow.sid()));
                    for (size_t c = 0; c = pixel_cols; ++c)
                    {
                        logger.write(' ' + str(watt2dBm(signal_power_lin[index++])));
                    }
                    logger.write("\n");
                }
            }
        }
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
        timeslot(args.timeslot.value()),
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
