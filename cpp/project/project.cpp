// project.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include "project.h"
#include "matplot.h"

using namespace std;

/* in case of errors, check this */
string sim_error = "";

vector<Cow> cows;
vector<Station> stations;

network_package::Input* sys;

vector<double> power_list;          // in watts
vector<vector<double>> alphatable;  // scan angle
string input_file;
vector<vector<Polar_Coordinates>> polar_data;

void setup()
{
    /* fill containers with test info */
    get_parameters(input_file);
    get_mobile_station_cords(mobile_station_pos);
    get_base_station_cords(base_station_pos);
    get_base_station_power(sys->ant.array_power_wtts);
    get_base_station_scan_alpha(sys->ant.array_scan_angle);
    get_allocation_per_bs(sys->ant.antcount_per_base);
    get_base_station_orientation(sys->ant.antenna_orientation);
    get_base_station_antspace(sys->ant.antenna_spacing);
    get_base_station_anntena_dims(sys->ant.antenna_dim_mtrs);

    for (unsigned power = sys->bs_Ptx_min; power <= sys->bs_Ptx_max; ++power)
    {
        power_list.emplace_back(dBm2watt(power));
    }

    for (unsigned i = 0; i < mobile_station_pos.size(); ++i)
    {
        stations.emplace_back(i, cows, sys->NF_watt);
    }

    auto bs_Ptx_min = dBm2watt(sys->bs_Ptx_min);
    auto bs_Ptx_max = dBm2watt(sys->bs_Ptx_max);
    auto bs_scan_min = deg2rad(sys->bs_scan_min);
    auto bs_scan_max = deg2rad(sys->bs_scan_max);

    /* setup the system */
    polar_data.resize(sys->base_stations);
    for (unsigned bs_id = 0; bs_id < sys->base_stations; ++bs_id)
    {

        unsigned idx = 0;
        polar_data[bs_id].resize(sys->mobile_stations);

        for (auto& mstation : mobile_station_pos)
        {
            long int diffx = mstation.x - base_station_pos[bs_id].x;
            long int diffy = mstation.y - base_station_pos[bs_id].y;
            polar_data[bs_id][idx++] = cart2pol(diffx, diffy);
        }

        cows.emplace_back(bs_id, sys, power_list, polar_data[bs_id]);
    }
}


/*
sinr_timeslot
--------------------------------------------------------
STA      SLOT0   SLOT1   SLOT2   SLOT3   SLOT4   SLOT5
--------------------------------------------------------
1
--------------------------------------------------------
2
--------------------------------------------------------
3
--------------------------------------------------------
4
--------------------------------------------------------
5
--------------------------------------------------------
6
--------------------------------------------------------
7
--------------------------------------------------------
..
..

*/

void run(PerfMon*& perf)
{
    using slot_t = unsigned;
    using cow_t = unsigned;
    using sta_t = unsigned;

    double SINR_CUT_OFF_VALID = log2lin(24); // dBm;
    //unordered_set<SimulationHelper, SimulationHelper::combo_hash> combocheck;

    /* run the simulation */
    SimulationHelper simparams(cows);
    vector<unsigned> select_stations;
    perf = new PerfMon(simparams, select_stations);

    int served = 0;
    unordered_map<sta_t, unordered_map<slot_t, cow_t>> station_binding;


    while (served < sys->mobile_stations)
    {
        simparams.timeslot = 0;
        /* setup the cows first such that the simulation is set with alphas and TX powers */
        for (int slot = 0; slot < simparams.timeslots; ++slot)
        {
            /* retrive the coefficients from base stations after antenna reconfig (alpha change) */

            select_stations = simparams.getBindings();

            /* change the antenna power and recalculate the Rx signal power + get SINR */
            simparams.setPower(); // random power in the range [-30, 30]

            /* change the scan angle of the base stations as per the bindings */
            simparams.setAlpha();

            /* check sinr for each receiving stations */
            for (unsigned bs_id = 0; bs_id < select_stations.size(); ++bs_id)
            {
                /* setup the stations with the COW bindings to calculate RX signal and interference */
                auto& station_id = select_stations[bs_id];

                auto& station = stations[station_id];
                station.setRX(bs_id);

                auto sinr = station.getSINR();
                if (sinr >= SINR_CUT_OFF_VALID) // check this settings
                {
                    perf->update(sinr);
                    ++served;
                }

            } //end of cow_id

            ++simparams.timeslot; // increment the timeslot to get the next parameters

        } //end of sta_id
    } //end of sta_id
}


int main(int argc, char** argv, char** envp)
{
    /* setup the simulation runtime parameters */
    setup();

    /* run the simulation and get the SNR table */
    PerfMon* perfmonitor;

    run(perfmonitor);
    delete perfmonitor;

    plot(mobile_station_pos);

    return 0;
}
