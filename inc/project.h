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

using cow_id_distribution = std::uniform_int_distribution<int>;


struct SimulationHelper
{
    std::vector<Cow>& cows;
    const unsigned& cow_count, ms_stations;
    const std::vector<std::vector<double>>& powers_lut;                 // TX power level from base station in integer dBm lut
    const std::vector<std::vector<double>>& alphas_lut;                 // Antenna array directions in integer rads lut
    const std::vector<std::vector<unsigned>>& binding_station_ids_lut;  // base_station to station id binding lut

    int      timeslot_idx;            // timeslot in the schedule
    unsigned timeslots;           // total timeslots
    unsigned nodes_rem;

    void incrementTimeslot()
    {
        ++timeslot_idx;
    }

    const unsigned& getTimeslots() const
    {
        return timeslots;
    }

    const std::vector<Cow>& getCows() const // this is for performance monitor
    {
        return cows;
    }

    /* gets called to get BS to MS bindings for a single timeslot. Output size: #bs */
    const std::vector<unsigned>& getBindings()
    {
        timeslot_idx = (timeslot_idx + 1) % timeslots;
        return binding_station_ids_lut[timeslot_idx];
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
    void resetPower()
    {
        timeslot_idx = (timeslot_idx + 1) % timeslots;
        set_power(powers_lut[timeslot_idx]);
    }

    /* get antenna array power level in each timeslot */
    void get_power(std::vector<double>& output)
    {
        timeslot_idx = (timeslot_idx + 1) % timeslots;
        output = powers_lut[timeslot_idx];
    }

    /* get permutations of power applied across base stations */
    bool get_perm_power(std::vector<double>& power_list)
    {
        return std::next_permutation(power_list.begin(), power_list.end());
    }

    /* set antenna array directivity before starting each simulation */
    void set_scan_dir()
    {
        timeslot_idx = (timeslot_idx + 1) % timeslots;
        auto& source = alphas_lut[timeslot_idx];

        for (unsigned c = 0; c < cows.size(); ++c)
        {
            cows[c].antennaUpdate(source[c]);
        }
    }

    SimulationHelper(std::vector<Cow>& cowlist,
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
        timeslot_idx(-1),
        timeslots(timeslot_count),
        nodes_rem(mobile_stations)
    {
        /* change the scan angle of the base stations as per the bindings */
        set_scan_dir();
    }
    //SimulationHelper(Input* args, std::vector<Cow>& c) :
    //    cows(c),
    //    cow_count(args->base_stations),
    //    binding_station_ids(station_ids),
    //    powers(args->ant.array_power_wtts),
    //    alphas(args->ant.array_scan_angle),
    //    timeslot(0)
    //{
    //
    //}
};

struct telemetry_t
{
    const double power;
    const double alpha;
    const Placements pos;
    const unsigned cow_id;
    const unsigned timeslot_idx;
    const unsigned sta_id;
    double sinr;

    telemetry_t(const Cow& cow, const unsigned& timeslot, const unsigned& sta, double& _sinr)
        :
        power(cow.getxPower()),
        alpha(cow.getAlpha()),
        pos(cow.getPosition()),
        cow_id(cow.sid()),
        timeslot_idx(timeslot),
        sta_id(sta),
        sinr(_sinr)
    {}
    telemetry_t(const Cow& cow, const unsigned& sta, double& _sinr)
        :
        power(cow.getxPower()),
        alpha(cow.getAlpha()),
        pos(cow.getPosition()),
        cow_id(cow.sid()),
        timeslot_idx(0),
        sta_id(sta),
        sinr(_sinr)
    {}
    telemetry_t() :
        power(0),
        alpha(0),
        pos(),
        cow_id(0),
        timeslot_idx(0),
        sta_id(0),
        sinr(0) {}
};

std::ostream& operator<<(std::ostream& out, telemetry_t const& data)
{
    out << std::fixed << std::fixed << std::setprecision(3);
    out << "power:\t" << data.power << " "
        << "alpha:\t" << data.alpha << " "
        << "placement:\t" << data.pos.x << "," << data.pos.y << " "
        << "cow_id:\t" << data.cow_id << " "
        << "timeslot:\t" << data.timeslot_idx << " "
        << "sta_id:\t" << data.sta_id;

    return out;
}

/* performance logging and analysis latter */
struct PerfMon
{
    const SimulationHelper& params;
    const std::vector<unsigned>& bindings;

    /* keep track of sinr and note the configurations and bindings in decreasing order of SINR */
    std::map<double, std::vector<telemetry_t>, std::greater<double>> sinrlist;

    const std::map<double, std::vector<telemetry_t>, std::greater<double>>& get_data() const
    {
        return sinrlist;
    }

    void update_timeslot(double sinr)
    {
        auto& cows = params.getCows();
        for (unsigned c = 0; c < params.cow_count; ++c)
        {
            sinrlist[sinr].emplace_back(telemetry_t(cows[c], params.timeslot_idx, bindings[c], sinr));
        }
    }

    void update(const unsigned& cow_id, double sinr)
    {
        sinrlist[sinr].emplace_back(telemetry_t(params.getCows()[cow_id], bindings[cow_id], sinr));
    }

    PerfMon(const SimulationHelper& parameters, std::vector<unsigned>& cowbindings) :
        params(parameters),
        bindings(cowbindings)
    {}
};


class Simulator
{
protected:
    Logger& logger;
    std::string sim_error;
    std::vector<Cow> cows;
    std::vector<Station> stations;
    PerfMon* perfmonitor;
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

    void run_multi_timeslot(SimulationHelper& simhelper)
    {
        using slot_t = unsigned;
        using cow_t = unsigned;
        using sta_t = unsigned;

        std::vector<unsigned> select_stations;
        perfmonitor = new PerfMon(simhelper, select_stations);

        unsigned served = 0;
        std::unordered_map<sta_t, std::unordered_map<slot_t, cow_t>> station_binding;

        std::vector<double> power_list;
        while (served < mobile_station_count)
        {
            simhelper.timeslot_idx = 0;
            /* setup the cows first such that the simulation is set with alphas and TX powers */
            for (unsigned slot = 0; slot < simhelper.timeslots; ++slot)
            {
                /* retrive the coefficients from base stations after antenna reconfig (alpha change) */

                select_stations = simhelper.getBindings();

                /* change the antenna power and recalculate the Rx signal power + get SINR */
                simhelper.set_power(power_list); // random power in the range [-30, 30]

                /* change the scan angle of the base stations as per the bindings */
                simhelper.set_scan_dir();

                /* check sinr for each receiving stations */
                for (unsigned bs_id = 0; bs_id < select_stations.size(); ++bs_id)
                {
                    /* setup the stations with the COW bindings to calculate RX signal and interference */
                    auto& station_id = select_stations[bs_id];

                    auto& station = stations[station_id];
                    station.set_rx(bs_id);

                    auto sinr = station.get_sinr();
                    if (sinr >= sinr_limit_linear) // check this settings
                    {
                        perfmonitor->update(bs_id, sinr);
                        ++served;
                    }

                } //end of cow_id

                ++simhelper.timeslot_idx; // increment the timeslot to get the next parameters

            } //end of sta_id
        } //end of sta_id
    }

    void run_single_timelot(SimulationHelper& simhelper)
    {
        using slot_t = unsigned;
        using cow_t = unsigned;
        using sta_t = unsigned;

        auto ms_station_list = simhelper.getBindings();

        perfmonitor = new PerfMon(simhelper, ms_station_list);

        /* create a new list of power nums and run calculations on its permutations */
        std::vector<double> power_list;
        simhelper.get_power(power_list);

        /* update the SINR for each scenario (base-station-count factorial permutations) */
        bool permutation_state = true;
        while (permutation_state)
        {
            /* change the antenna power across all base stations */
            simhelper.set_power(power_list);

            for (size_t bs = 0; bs < base_station_count; ++bs)
            {
                auto& station = stations[ms_station_list[bs]];
                station.set_rx(bs);

                perfmonitor->update(bs, station.get_sinr());
            }

            permutation_state = simhelper.get_perm_power(power_list);
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

        if (timeslot_count > 1)
            run_multi_timeslot(simparams);
        else
            run_single_timelot(simparams);

    }


    const PerfMon& get_perf() const
    {
        return *perfmonitor;
    }

    Simulator(const MyArgs& args, Logger& ilogger)
        :
        logger(ilogger),
        sim_error(""),
        perfmonitor(nullptr),
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