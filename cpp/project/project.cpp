// project.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include "project.h"
#include "matplotlibcpp.h"
#include <cmath>
#include <queue>

#define CPRACTICE
//#define NODUPLICATES

/* in case of errors, check this */
std::string sim_error = "";
Random randomgen;


using namespace std;
namespace plt = matplotlibcpp;

vector<Cow> cows;
vector<Station> stations;
const auto& station_cords = SEATING_LOCATION;

/* plot the charts */
void plot(const vector<Placements>& input, string title = "", bool hold = false)
{

//#ifndef _DEBUG
    /* plot only */
    vector<int> locx, locy;
    for (auto& loc : input)
    {
        locx.emplace_back(loc.x);
        locy.emplace_back(loc.y);
    }

    //if (hold) plt::hold(true);

    plt::scatter(locx, locy, { {"color", "green"}, {"size", "20"} });
    //plt::
    if (title.size()) plt::title(title);
    plt::legend();
    plt::show();
//#endif
}



unordered_map<int, double> Cow::dBm2watts_lut;
vector<vector<double>> alpha_lut(SCAN_ALPHA_LUT.size(), vector<double>(SCAN_ALPHA_LUT[0].size()));

void setup()
{
    /* setup the system */
    for (int i = 0; i < COW_COUNT; ++i)
    {
        vector<Polar_Coordinates> polar_data(SEATING_LOCATION.size());
        unsigned idx = 0;
        for (auto& mstation : SEATING_LOCATION)
        {
            auto diffx = mstation.x - COW_LOCATION[i].x;
            auto diffy = mstation.y - COW_LOCATION[i].y;
            polar_data[idx++] = cart2pol(diffx, diffy);
        }

        cows.emplace_back(Cow(i, polar_data, BS_ANTENNA_ASSIGNMENTS[i]));  // initialize cow instance

    }

    for (unsigned i = 0; i < SEATING_LOCATION.size(); ++i)
    {
        stations.emplace_back(i, cows);
    }

    for (int power = MIN_POWER_dBm; power <= MAX_POWER_dBm; ++power)
    {
        Cow::dBm2watts_lut.emplace(power, dBm2watts(power));
    }

    unsigned i = 0, j = 0;
    for (auto& alpha_list : SCAN_ALPHA_LUT)
    {
        j = 0;
        for (auto& alpha : alpha_list)
        {
            alpha_lut[i][j++] = deg2rad(alpha);
        }
        ++i;
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

vector<telemetry_t> run()
{
    /*
        for power1 = -30:30
        p1 = 10^((power1-30)/10);
        for power2 = -30:30
            p2 = 10^((power2-30)/10);
            for power3 = -30:30
                p3 = 10^((power3-30)/10);

                ptable = [p1*table(1,1) / addrms([p2*table(1,2) p3*table(1,3) THER])
                    p2*table(1,2) / addrms([p1*table(1,1) p3*table(1,3) THER])
                    p3*table(1,3) / addrms([p2*table(1,2) p1*table(1,1) THER])
                    p1*table(2,1) / addrms([p2*table(2,2) p3*table(2,3) THER])
                    p2*table(2,2) / addrms([p1*table(2,1) p3*table(2,3) THER])
                    p3*table(2,3) / addrms([p2*table(2,2) p1*table(2,1) THER])
                    p1*table(3,1) / addrms([p2*table(3,2) p3*table(3,3) THER])
                    p2*table(3,2) / addrms([p1*table(3,1) p3*table(3,3) THER])
                    p3*table(3,3) / addrms([p2*table(3,2) p1*table(3,1) THER])];

                maxvec(1,idx) = power1;
                maxvec(2,idx) = power2;
                maxvec(3,idx) = power3;
                maxvec(4,idx) = max(ptable(1:3));
                maxvec(5,idx) = max(ptable(4:6));
                maxvec(6,idx) = max(ptable(7:9));
                idx = idx + 1;

%                 if ( ptable(1) >= cutofflinear || ptable(2) >= cutofflinear || ptable(3) >= cutofflinear )
%                     display('hit')
%                 else
%
%                     continue;
%                 end
%
%
%                 if ( ptable(4) >= cutofflinear || ptable(5) >= cutofflinear || ptable(6) >= cutofflinear )
%                     display('hit')
%                 else
%
%                     continue;
%                 end
%
%                 if ptable(7) >= cutofflinear || ptable(8) >= cutofflinear || ptable(9) >= cutofflinear
% %                     display(['all hit p1:' num2str(power1) ' p2:' num2str(power2) ' p3:' num2str(power3)])
% %                      (idx) = max(ptable(7:9));
%                     breakflag = 1;
%                     cowpowers = [power1 power2 power3];

%                     break;
%                 end

            end
%             if(breakflag == 1)
%                 break;
%             end
        end
%         if(breakflag == 1)
%             break;
%         end
    end
    */

    double thermal_noise = dBm2watts(thermal_n_db);
    const int slots = ceil(station_cords.size() / COW_COUNT);
    const unsigned power_range = MAX_POWER_dBm - MIN_POWER_dBm + 1;

    using slot_t = unsigned;
    using cow_t = unsigned;
    using sta_t = unsigned;
    auto& nodes = SEATING_LOCATION;


    double SINR_CUT_OFF_VALID = log2lin(24); // dBm;
    unordered_set<SimulationHelper, SimulationHelper::combo_hash> combocheck;

    /* run the simulation */
    SimulationHelper simparams(cows);
    vector<unsigned> bindings;
    PerfMon perf(simparams, bindings);

    int served = 0;
    vector<telemetry_t> sinr_list;
    unordered_map<sta_t, unordered_map<slot_t, cow_t>> station_binding;



    while (served < stations.size())
    {
            /* setup the cows first such that the simulation is set with alphas and TX powers */

                //#ifdef NODUPLICATES
             /* retrive the coefficients from base stations after antenna reconfig (alpha change) */
            simparams.setalpha();

            for (int slot = 0; slot < simparams.timeslot; ++slot)
            {
                /* only change the antenna power and recalculate the Rx signal power + get SINR */

                simparams.setpower(); // random power in the range [-30, 30]


                bindings = simparams.getbindings();

                /* check sinr for each receiving stations */
                for (unsigned bs_id = 0; bs_id < bindings.size(); ++bs_id)
                {
                    auto staid = bindings[bs_id];

                    /* setup the stations with the COW bindings to calculate RX signal and interference */

                    auto& station = stations[staid];
                    station.setRX(bs_id);

                    auto& sinr = station.getSINR();
                    if (sinr >= SINR_CUT_OFF_VALID)// && station_served[sta_id] == false) // check this settings
                    {
                        perf.update(sinr);
                        ++served;
                    }

                } //end of cow_id
            } //end of sta_id
    } //end of sta_id

    vector<bool> ms_served(nodes.size(), false);

    queue<telemetry_t> qlist;
    int slot = 0;

    while (sinr_list.size() < nodes.size())
    {
        vector<bool> bs_served(cows.size(), false);

        ++slot;
        for (auto data_list : sinr_helper) //list of same sinr
        {
            for (int i = 0; i < data_list.second.size(); ++i);
            {
                //qlist.push(data_list.second[i]);
            }

            while (qlist.size())
            {
                auto& data = qlist.front();

                if (ms_served[data.sta_id] == false || bs_served[data.cow_id] == false)
                {
                    sinr_list.emplace_back(data);

                    if (sinr_list.size() == nodes.size())
                        return sinr_list;

                    ms_served[data.sta_id] = true;
                    bs_served[data.cow_id] = true;
                    qlist.pop();
                }
                else if (ms_served[data.sta_id] == false || bs_served[data.cow_id] == true)
                {
                    qlist.push(data);
                }
            }
        }
    }

    return sinr_list;
}


int main(int argc, char** argv, char** envp)
{
    /* setup the simulation runtime parameters */
    setup();

    /* run the simulation and get the SNR table */
    run();

    plot(SEATING_LOCATION);

    return 0;
}