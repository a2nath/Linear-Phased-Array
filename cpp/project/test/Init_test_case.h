#pragma once
#include <vector>
#include "coordinates.h"

#define INIT_TEST

namespace init_test_case
{
    namespace Network_System
    {
        const double   frequency = 1900e6;   // center frequency
        const double   bandwidth = 20e6;     // channel bandwidth in Mz
        const double   SymbRate = 3.84e6;    // symbol rate
        const double   BlocksPerSymb = 768;  // blocks per symbol
        const double   lambda = C / frequency;
        const double   height = 200;         //antenna height in meters
        const double   ms_Grx = -2;          // dBi
        const double   system_noise = 5;     // System noise in dB;
        const unsigned base_stations = 3;    // number of base stations

        const std::vector<Placements> NODE_LOCATIONS = // seating of all mobile stations in meters
        {
            {250, 500}, {350, 500}, {450, 500}, {550, 500}, {650, 500}, {750, 500},
            {250, 600}, {350, 600}, {450, 600}, {550, 600}, {650, 600}, {750, 600},
                              {450, 325}, {500, 325}, {550, 325}
        };
    }

    namespace Antenna_Array
    {
        const int MIN_POWER_dBm = -30;
        const int MAX_POWER_dBm = +30;
        const int MIN_SCANANGLE_deg = -90;
        const int MAX_SCANANGLE_deg = +90;

        const std::vector<int> BS_ANTENNA_ASSIGNMENTS = { 5, 3, 5 };        // number of antenna allocation for each BS
        const std::vector<double> BS_PLACEMENT_THETAC = { 75, 90, 105 };    // degrees
        const std::vector<double> BS_ANTENNA_SPACING = { .3, .4, .3 };      // meters
        const std::pair<double, double> ANTENNA_DIMS = { 20e-2, 40e-2 }; // meters
        
                                                                         // dBm
        const std::vector<std::vector<double>> COW_POWER_LUT = {
            {+10.0, +65.0, +35.0},
            {-10.0, -65.0, -35.0},
            {-75.0, -40.0, +75.0},
            {+40.0, +65.0, -40.0},
            {+40.0, -65.0, -40.0}
        };
        // degrees
        const std::vector<std::vector<double>> SCAN_ALPHA_LUT = {
            {+10.0, +65.0, +35.0},
            {-10.0, -65.0, -35.0},
            {-75.0, -40.0, +75.0},
            {+40.0, +65.0, -40.0},
            {+40.0, -65.0, -40.0}
        };

        const std::vector<Placements> COW_LOCATION =
        {
            {250, 200},
            {500, 180},
            {750, 200}
        };

        /* Setup the environment */
        const std::vector<unsigned> BS_STAID_SELECTION = {

            /* BS [0] [1] [2]  */
                   2, 13, 11,
                   8, 15, 05,
                   1, 14, 06,
                   7, 04, 10,
                   9, 03, 12
        };
    }
}