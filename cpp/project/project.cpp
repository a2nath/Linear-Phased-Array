// project.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

#include "project.h"
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

vector<Cow> cows;
const auto& station_cords = SEATING_LOCATION;

class Solution {
public:
    vector<int> twoSum(vector<int> nums, int target) {

        unordered_map<int, int> m;
        int index1, index2;
        int diff;

        for (int i = 0; i < nums.size(); ++i)
        {
            diff = target - nums[i];


            //if (m[nums[i]] > 0)
            //{
            //    index1 = m[nums[i]] - 1;
            //    index2 = i;
            //    break;
            //}
            if (m[nums[i]] > 0)
            {
                index1 = m[nums[i]] - 1;
                index2 = i;
                break;
            }

            m[diff] = i + 1;
        }



        return vector<int>{index1, index2};
    }
};

#include <unordered_set>

/* Setup the environment */
vector<vector<int>> NODE_SELECTION = {
    {2, 13, 11},
    {8, 15, 05},
    {1, 14, 06},
    {7, 04, 10},
    {9, 03, 12}
};

void run()
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
%                     powers = [power1 power2 power3];

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

    double matrix;                                // array matrix
    double thermal_noise = log2lin(thermal_n_db);
    const int slots = ceil(station_cords.size() / COW_COUNT);


    vector<vector<vector<unordered_map<int, vector<double>>>>> sinrtable(MAX_POWER_dBm - MIN_POWER_dBm,
        vector<vector<unordered_map<int, vector<double>>>>(MAX_POWER_dBm - MIN_POWER_dBm,
            vector<unordered_map<int, vector<double>>>(MAX_POWER_dBm - MIN_POWER_dBm)));

    vector<double> txpower_in_linear;
    for (int power = MIN_POWER_dBm; power <= MAX_POWER_dBm; ++power)
    {
        txpower_in_linear.emplace_back(dBm2watts(power));
    }

    vector<vector<double>> alpha_lut(SCAN_ALPHA_LUT.size(), vector<double>(SCAN_ALPHA_LUT[0].size()));
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

    int power_range = txpower_in_linear.size();


    //unordered_set<int> stations;
    //int i = 0;
    //for (auto& station : station_cords) { stations.insert(i++); }
    //
    //long unsigned iter = 0;
    //unordered_map<int, int> power_it;
    //int cow_it = 0;
    //
    //while (iter < power_range ^ cows.size())
    //{
    //    int cow_it = 0;
    //    while (cow_it < cows.size())
    //    {
    //        cows[cow_it].setpower(txpower_in_linear[power_it[cow_it]]);
    //        ++cow_it;
    //    }
    //
    //    ++power_it[cow_it];
    //
    //    if (power_it[cow_it] == txpower_in_linear.size())
    //    {
    //
    //    }
    //
    //    if (iter % txpower_in_linear.size() == 0)
    //    {
    //        ++cow_it;
    //        cows[]
    //
    //            cow_it = cow_it % cows.size();
    //        ++iter;
    //    }
    //}
    using slot_t = unsigned;
    using cow_t = unsigned;
    unordered_map<cow_t, int> iters;



    /* run the simulation */
    for (int slot = 0; slot < slots; ++slot)
    {
        unordered_map<cow_t, const vector<double>&> hmatrix;

        /* retrive the coefficients from base stations after antenna reconfig (alpha change) */
        for (unsigned cow_it = 0; cow_it < cows.size(); ++cow_it)
        {
            cows[cow_it].update(alpha_lut[slot][cows[cow_it].sid()]);
        }

       /* only change the antenna power and recalculate the Rx signal power + get SINR */
        double max_sta1_received_signal = 0,
            max_sta2_received_signal = 0,
            max_sta3_received_signal = 0;

        //int iter;
        //
        //unordered_map<int, vector<double>&> powers;
        //
        //for (int cow_it = 0; cow_it < cows.size(); ++cow_it)
        //{
        //    powers[cow_it] = txpower_in_linear;
        //}
        //
        //while (iter != power_range ^ cows.size())
        //{
        //    for (int power_idx = 0; power_idx < txpower_in_linear.size(); ++power_idx)
        //    {
        //        //for (int cow_it = 0; cow_it < cows.size(); ++cow_it)
        //        //{
        //
        //        cows[0].setpower(iter % power_range);// txpower_in_linear[(power_idx + cow_it) % txpower_in_linear.size()]);
        //        cows[0].setpower((iter / power_range) % power_range);// txpower_in_linear[(power_idx + cow_it) % txpower_in_linear.size()]);
        //        cows[0].setpower((iter / power_range / power_range) % power_range);// txpower_in_linear[(power_idx + cow_it) % txpower_in_linear.size()]);
        //        //}
        //
        //    }
        //}


        unordered_map<cow_t, vector<double>> sinrtable;

        for (iters[0] = 0; iters[0] < txpower_in_linear.size(); ++iters[0])
        {
            for (iters[1] = 0; iters[1] < txpower_in_linear.size(); ++iters[1])
            {
                for (iters[2] = 0; iters[2] < txpower_in_linear.size(); ++iters[2])
                {
                    double& power0 = txpower_in_linear[iters[0]];
                    double& power1 = txpower_in_linear[iters[1]];
                    double& power2 = txpower_in_linear[iters[2]];
                    vector<double> powers(cows.size());


                    //for (iters[0] = 0; iters[0] < txpower_in_linear.size(); ++iters[0])
                    //{
                    //    for (int cow_it = 0; cow_it < cows.size(); ++cow_it)
                    //    {
                    //        cows[cow_it].setpower(txpower_in_linear[cow_it]);
                    //    }
                    //}

                    auto& nodes = NODE_SELECTION[slot];


                    for (int cow_it = 0; cow_it < cows.size(); ++cow_it)
                    {
                        vector<double> rx_signals(nodes.size(), 0.00);

                        cows[cow_it].setpower(txpower_in_linear[iters[cow_it]]);

                        /* precalculate power x hmatrix coefficient */
                        for (int sta_id = 0; sta_id < nodes.size(); ++sta_id)
                        {
                            cows[cow_it].get_signal_level(nodes[sta_id], rx_signals[sta_id]);
                        }

                        /* calcalate the sinr and choose best one */
                        auto& bs_sinr_list = sinrtable[cow_it];
                        for (int sta_id = 0; sta_id < nodes.size(); ++sta_id)
                        {
                            bs_sinr_list.emplace_back(thermal_noise);
                            for (int sta_offset = sta_id; sta_offset < nodes.size(); ++sta_offset)
                            {
                                bs_sinr_list.back() += rx_signals[(sta_offset + 1) % cows.size()];
                            }
                            bs_sinr_list.back() = rx_signals[sta_id] / bs_sinr_list.back();
                        }
                    }


                    /*
                    sinr[0].emplace_back(rx_signals[0][0] /
                        (rx_signals[0][1] + rx_signals[0][2] + thermal_noise));

                    sinr[0].emplace_back(rx_signals[0][1] /
                        (rx_signals[0][0] + rx_signals[0][2] + thermal_noise));

                    sinr[0].emplace_back(rx_signals[0][2] /
                        (rx_signals[0][0] + rx_signals[0][1] + thermal_noise));

                    sinr[1].emplace_back(rx_signals[1][0] /
                        (rx_signals[1][1] + rx_signals[1][2] + thermal_noise));

                    sinr[1].emplace_back(rx_signals[1][1] /
                        (rx_signals[1][0] + rx_signals[1][2] + thermal_noise));

                    sinr[1].emplace_back(rx_signals[1][2] /
                        (rx_signals[1][0] + rx_signals[1][1] + thermal_noise));

                    sinr[2].emplace_back(rx_signals[2][0] /
                        (rx_signals[2][1] + rx_signals[2][2] + thermal_noise));

                    sinr[2].emplace_back(rx_signals[2][1] /
                        (rx_signals[2][0] + rx_signals[2][2] + thermal_noise));

                    sinr[2].emplace_back(rx_signals[2][2] /
                        (rx_signals[2][0] + rx_signals[2][1] + thermal_noise));
                        */

                }
            }
        }

    }
}

int main(int argc, char** argv, char** envp)
{
    //for (char** env = envp; *env != 0; env++)
    //{
    //    char* thisEnv = *env;
    //    printf("%s\n", thisEnv);
    //}
    //PyImport_ImportModule("numpy");
    //PyObject* numpy = PyImport_ImportModule("numpy.core._multiarray_umath");
    //if (numpy == nullptr)
    //    throw runtime_error("Numpy not found");

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

    /* plot only */
    vector<int> seatx, seaty;
    for (auto& seat : SEATING_LOCATION)
    {
        seatx.emplace_back(seat.x);
        seaty.emplace_back(seat.y);
    }


    plt::scatter(seatx, seaty);
    plt::title("Seating plan");
    plt::legend();
    plt::show();

    run();
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started:
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
