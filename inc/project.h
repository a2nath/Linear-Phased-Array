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


using cow_id_distribution = std::uniform_int_distribution<int>;

struct SimulationHelper
{
    std::vector<Cow>& cows;
    unsigned cow_count;
    std::vector<std::vector<unsigned>>& binding_station_ids_lut;  // base_station to station id binding lut
    std::vector<unsigned>& binding_station_ids;                   // base_station to station id binding
    std::vector<std::vector<double>>& powers_lut;                 // TX power level from base station in integer dBm lut
    std::vector<double>& powers;                                  // TX power level from base station in integer dBm
    std::vector<std::vector<double>>& alphas_lut;                 // Antenna array directions in integer rads lut
    std::vector<double>& alphas;                                  // Antenna array directions in integer rads

    unsigned timeslot;            // timeslot in the schedule
    unsigned timeslots;           // total timeslots
    unsigned nodes_rem;

    const std::vector<Cow>& getCows() const // this is for performance monitor
    {
        return cows;
    }

    /* gets called to get BS to MS bindings for a single timeslot. Output size: #bs */
    std::vector<unsigned>& getBindings()
    {
        return binding_station_ids_lut.size() ? binding_station_ids_lut[timeslot] : binding_station_ids;
    }

    /* set antenna array power level in each timeslot */
    inline void setPower(std::vector<double>& power_list) const
    {
        for (unsigned i = 0; i < power_list.size(); ++i)
        {
            cows[i].setPower(power_list[i]);
        }
    }

    /* set antenna array power level in each timeslot */
    void resetPower() const
    {
        setPower(powers_lut.size() ? powers_lut[timeslot] : powers);
    }

    /* get antenna array power level in each timeslot */
    void getPower(std::vector<double>& output) const
    {
        output = powers_lut.size() ? powers_lut[timeslot] : powers;
    }

    /* get permutations of power applied across base stations */
    bool get_perm_power(std::vector<double>& power_list)
    {
        return std::next_permutation(power_list.begin(), power_list.end());
    }

    /* set antenna array directivity before starting each simulation */
    void setAlpha()
    {
        auto& source = alphas_lut.size() ? alphas_lut[timeslot] : alphas;
        for (unsigned c = 0; c < cows.size(); ++c)
        {
            cows[c].antennaUpdate(source[c]);
        }
    }

    SimulationHelper(Input* args, std::vector<Cow>& c) :
        cows(c),
        cow_count(args->base_stations),
        binding_station_ids_lut(station_ids_lut),
        binding_station_ids(station_ids),
        powers_lut(args->ant.array_power_wtts_lut),
        powers(args->ant.array_power_wtts),
        alphas_lut(args->ant.array_scan_angle_lut),
        alphas(args->ant.array_scan_angle),
        timeslot(0),
        timeslots(std::ceil(args->mobile_stations / cow_count)),
        nodes_rem(args->mobile_stations)
    {
        /* change the scan angle of the base stations as per the bindings */
        setAlpha();
    }
};

struct telemetry_t
{
    const double power;
    const double alpha;
    const Placements pos;
    const unsigned cow_id;
    const unsigned timeslot;
    const unsigned sta_id;
    double sinr;

    telemetry_t(const Cow& cow, const unsigned& timeslot, const unsigned& sta, double& _sinr)
        :
        power(cow.getxPower()),
        alpha(cow.getAlpha()),
        pos(cow.getPosition()),
        cow_id(cow.sid()),
        timeslot(timeslot),
        sta_id(sta),
        sinr(_sinr)
    {}
    telemetry_t(const Cow& cow, const unsigned& sta, double& _sinr)
        :
        power(cow.getxPower()),
        alpha(cow.getAlpha()),
        pos(cow.getPosition()),
        cow_id(cow.sid()),
        timeslot(0),
        sta_id(sta),
        sinr(_sinr)
    {}
    telemetry_t() :
        power(0),
        alpha(0),
        pos(),
        cow_id(0),
        timeslot(0),
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
        << "timeslot:\t" << data.timeslot << " "
        << "sta_id:\t" << data.sta_id;

    return out;
}

/* performance logging and analysis latter */
struct PerfMon
{
    const SimulationHelper& params;
    const std::vector<unsigned>& bindings;
    std::vector<telemetry_t> output;

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
            sinrlist[sinr].emplace_back(telemetry_t(cows[c], params.timeslot, bindings[c], sinr));
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
