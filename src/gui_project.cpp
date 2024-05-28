#include "visuals.h"
// project.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include "project.h"
//#include "matplot.h"
using namespace std;

/* in case of errors, check this */
string sim_error = "";

vector<Cow> cows;
vector<Station> stations;

Input* args = nullptr;
PerfMon* perfmonitor = nullptr;
string jason_file = "";

vector<double> power_list;          // in watts
vector<vector<double>> alphatable;  // scan angle

void setup(string config_file = "")
{
    /* if the program is called using a test file or a directory of .json files */
    if (config_file.size())
    {
        /* populate a container with args data */
        JasonHelper config(config_file);

        /* take the same args from above and set the program config */
        config.parse_data(args);
    }

    for (unsigned power = args->bs_Ptx_min; power <= args->bs_Ptx_max; ++power)
    {
        power_list.emplace_back(dBm2watt(power));
    }

    for (unsigned i = 0; i < mobile_station_pos.size(); ++i)
    {
        stations.emplace_back(i, cows[i], args->NF_watt);
    }

    auto bs_Ptx_min = dBm2watt(args->bs_Ptx_min);
    auto bs_Ptx_max = dBm2watt(args->bs_Ptx_max);
    auto bs_scan_min = deg2rad(args->bs_scan_min);
    auto bs_scan_max = deg2rad(args->bs_scan_max);

    /* setup the system location with respect to each base */
    for (unsigned bs_id = 0; bs_id < args->base_stations; ++bs_id)
    {
        vector<Polar_Coordinates> polar_data(args->mobile_stations);

        for (unsigned ms_id = 0; ms_id = args->mobile_stations; ++ms_id)
        {
            auto& client = mobile_station_pos[ms_id];
            long int diffx = client.x - base_station_pos[bs_id].x;
            long int diffy = client.y - base_station_pos[bs_id].y;
            polar_data[ms_id] = cart2pol(diffx, diffy);
        }

        cows.emplace_back(bs_id, args, power_list, move(polar_data));
    }
}

void run(SimulationHelper& helper)
{
    using slot_t = unsigned;
    using cow_t = unsigned;
    using sta_t = unsigned;

    double SINR_CUT_OFF_VALID = log2lin(24); // dBm;

    /* run the simulation */

    auto ms_station_list = helper.getBindings();
    helper.incrementTimeslot(); // update the timeslot

    perfmonitor = new PerfMon(helper, ms_station_list);

    /* create a new list of power nums and run calculations on its permutations */
    vector<double> power_list;
    helper.getPower(power_list);

    bool ok = power_list.size() > 0;

    /* update the SINR for each scenario (base-station-count factorial permutations) */
    while (ok)
    {
        /* change the antenna power across all base stations */
        helper.setPower(power_list);

        for (unsigned bs = 0; bs < args->base_stations; ++bs)
        {
            auto& station = stations[ms_station_list[bs]];
            station.setRX(cows);

            auto sinr = station.getSINR();
            perfmonitor->update(bs, sinr);
        }

        ok = helper.get_perm_power(power_list);
    }
}

void run_across_timslots()
{
    using slot_t = unsigned;
    using cow_t = unsigned;
    using sta_t = unsigned;

    double SINR_CUT_OFF_VALID = log2lin(24); // dBm;
    //unordered_set<SimulationHelper, SimulationHelper::combo_hash> combocheck;

    /* run the simulation */
    SimulationHelper simparams(args, cows);
    vector<unsigned> select_stations;
    perfmonitor = new PerfMon(simparams, select_stations);

    unsigned served = 0;
    unordered_map<sta_t, unordered_map<slot_t, cow_t>> station_binding;

    vector<double> power_list;
    while (served < args->mobile_stations)
    {
        simparams.timeslot = 0;
        /* setup the cows first such that the simulation is set with alphas and TX powers */
        for (unsigned slot = 0; slot < simparams.timeslots; ++slot)
        {
            /* retrive the coefficients from base stations after antenna reconfig (alpha change) */

            select_stations = simparams.getBindings();

            /* change the antenna power and recalculate the Rx signal power + get SINR */
            simparams.setPower(power_list); // random power in the range [-30, 30]

            /* change the scan angle of the base stations as per the bindings */
            simparams.setAlpha();

            /* check sinr for each receiving stations */
            for (unsigned bs_id = 0; bs_id < select_stations.size(); ++bs_id)
            {
                /* setup the stations with the COW bindings to calculate RX signal and interference */
                auto& station_id = select_stations[bs_id];

                auto& station = stations[station_id];
                station.setRX(cows);

                auto sinr = station.getSINR();
                if (sinr >= SINR_CUT_OFF_VALID) // check this settings
                {
                    perfmonitor->update(bs_id, sinr);
                    ++served;
                }

            } //end of cow_id

            ++simparams.timeslot; // increment the timeslot to get the next parameters

        } //end of sta_id
    } //end of sta_id
}

void close()
{
    delete perfmonitor;
    delete args;
}

void plot(vector<Signal2D>& signal_strength)
{
    Visuals visuals(2560, 1440, "super awesome title");

    Signal2D heatOutput;

    int bs_id = 0;

    visuals.draw_init(args->ant.antenna_spacing, args->lambda, cows[bs_id], heatOutput);
    visuals.show_eirp(heatOutput);
    visuals.set_nodes(base_station_pos, mobile_station_pos, Dimensions(1000.f, 1000.f));
    visuals.run();
}

int main(int argc, char** argv, char** envp)
{
    //thread gui_thread(graphics_stuff);
    //gui_thread.join();

    auto program_args = argparse::parse<MyArgs>(argc, argv);
    program_args.print();

    if (program_args.isDefault())
    {
        Defaults program(args);
    }
    /* run test from the args */
    else
    {
        /* run once with specific parameters */
        program_args.interprete();
        program_args.parse_data(args);
    }

    for (auto& test : program_args.json_files)
    {
        /* setup the simulation runtime parameters */
        setup(test);
        /* run the simulation and get the SNR table */

        SimulationHelper simhelper(args, cows);
        for (auto i = 0; i < simhelper.getTimeslots(); ++i)
        {
            run(simhelper);
            //plot();

        }
    }

    for (auto& entrylist : perfmonitor->get_data())
    {
        if (program_args.showgui)
        {

        }

        for (auto& entry : entrylist.second)
        {
            cout << "SINR:" << entrylist.first
                << " " << entry << endl;

        }
    }

    close();
    return 0;
}
