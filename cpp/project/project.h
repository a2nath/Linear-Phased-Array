#pragma once
#include "common.h"
#include <cassert>

const double frequency = 1900e6;   // center frequency
const double bandwidth = 20e6;     // channel bandwidth
const double SymbRate = 3.84e6;    // symbol rate
const double BlocksPerSymb = 768;  // blocks per symbol
const double lambda = C / frequency;
const double height = 200;         //antenna height in meters
const double ms_Grx = -2;// dBi

const int MIN_POWER_dBm = -30;
const int MAX_POWER_dBm = +30;
const std::pair<double, double> ANTENNA_DIMS = { 20e-2, 40e-2 }; // meters

const double system_noise = 5; // System noise in dB;
const double thermal_n_db = -174 + round(10 * log10(bandwidth)) + system_noise;  // Thermal noise floor in dBm


const std::vector<std::vector<double>> SCAN_ALPHA_LUT = {        // degrees
    {+10.0, +65.0, +35.0},
    {-10.0, -65.0, -35.0},
    {-75.0, -40.0, +75.0},
    {+40.0, +65.0, -40.0},
    {+40.0, -65.0, -40.0}
};

const std::vector<cords> SEATING_LOCATION = // seating of all mobile stations in meters
{
    {250, 500}, {350, 500}, {450, 500}, {550, 500}, {650, 500}, {750, 500},
    {250, 600}, {350, 600}, {450, 600}, {550, 600}, {650, 600}, {750, 600},
                      {450, 325}, {500, 325}, {550, 325}
};

const std::vector<cords> COW_LOCATION =
{
    {250, 200},
    {500, 180},
    {750, 200}
};

const int COW_COUNT = 3;
const int BS_ANTENNA_ASSIGNMENTS[COW_COUNT] = { 5, 3, 5 };        // number of antenna allocation for each BS
const double BS_PLACEMENT_THETAC[COW_COUNT] = { 75, 90, 105 };    // degrees
const double BS_ANTENNA_SPACING[COW_COUNT] = { .3, .4, .3 };      // meters

//const unsigned TIMESLOTS = 5;

class AntennaSystem // think of this as multiple antennas each with the following
{
    const Dimension dims;
    std::vector<double> Gtx;
    int panel_count;
public:
    const Dimension size() const { return dims; }
    const std::vector<double>& getgain() const { return Gtx; }
    const int count() const { return panel_count; }
    void setgain(double gain) { Gtx.emplace_back(panel_count * gain); }
    AntennaSystem(double Lx, double Ly, double _count) : dims(Lx, Ly), panel_count(_count) {}
};

// Array of Antenna
class AAntenna
{
    friend class Cow;
    //static unordered_map<int, vector<Polar_Coordinates>> station_logistical_data; // seating location polar data relative to base station location
    const int antenna_id;
    const double ms_Grx_lin;
    //std::vector<double> theta_list_minus_thetaC_list; // list of angles in relation to the other stations
    AntennaSystem antennas;
    int count;                             // number of panels
    double power;                          // array power
    double theta_c;                        // array orientation
    double spacing;                        // array spacing
    double alpha;                          // array direction
    std::vector<double> phee_minus_alpha_list;
    std::vector<double> pathloss_list;
    std::vector<double> hmatrix;  // hmatrix from BS pov


public:
    /* update the antenna array from updated power and scan angle */
    void update(double alpha)//, std::vector<double>& hmatrix_output)
    {
        auto antenna_count = antennas.count();
        auto& antgain_factor = antennas.getgain();

        /* update the antenna gain Gtx */
        for (unsigned idx = 0; idx < phee_minus_alpha_list.size(); ++idx)
        {
            double phee = (phee_minus_alpha_list[idx] + alpha) / 2;

            double sin_term = antenna_count * sin(phee);

            assert(sin_term != 0);

            double gain_factor = antgain_factor[idx] * pow(sin(antenna_count * phee) / sin_term, 2);

            /* update the channel matrix */
            hmatrix[idx] = ms_Grx_lin * gain_factor / pathloss_list[idx];
        }
    }

    const std::vector<double>& getmatrix() const
    {
        return hmatrix;
    }

    AAntenna(int id, int antenna_count, const std::vector<Polar_Coordinates>& polar_sta_data) :
        antenna_id(id), ms_Grx_lin(log2lin(ms_Grx)),
        antennas(AntennaSystem(ANTENNA_DIMS.first, ANTENNA_DIMS.second, antenna_count)),
        theta_c(deg2rad(BS_PLACEMENT_THETAC[id])), spacing(BS_ANTENNA_SPACING[id]),
        hmatrix(polar_sta_data.size(), 0.0)
    {
        double phee_temp = (2 * M_PIl * spacing / lambda);
        double pl_temp_meters = 4 * M_PIl / lambda;
        auto& size = antennas.size();

        for (auto& station_data : polar_sta_data)
        {
            double theta_minus_thetaC = station_data.theta - theta_c;
            double m = (M_PIl * size.x * sin(theta_minus_thetaC)) / lambda;
            double singleant_gain = (10 * size.x * size.y / pow(lambda, 2)) * pow((1 + cos(theta_minus_thetaC)) / 2, 2);

            if (m != 0)
            {
                singleant_gain *= pow(sin(m) / m, 2);
            }

            antennas.setgain(singleant_gain);
            phee_minus_alpha_list.emplace_back(phee_temp * sin(theta_minus_thetaC));
            pathloss_list.emplace_back(pow(pl_temp_meters * station_data.hype, 2));
        }
    }
};


class Cow
{
    const int station_id;
    AAntenna antenna;
    double& power;
    //std::vector<double>* matrix;
    Coordinates location;

public:

    /* set Gtx power in linear */
    void setpower(double& antenna_power)
    {
        power = antenna_power;
    }

    /* return Gtx power in linear */
    const double& getxpower() const
    {
        return power;
    }
    const int sid() const
    {
        return station_id;
    }

    /* update the parameters of the Base Station and get channel state to each station */
    void update(double alpha)
    {
        antenna.update(alpha);
        //return antenna.getmatrix();
    }

    inline void get_signal_level(unsigned node_id, double& signal_level_lin)
    {
        signal_level_lin = antenna.hmatrix[node_id] * power;
    }

    Cow(int id, const std::vector<Polar_Coordinates>& polar_data, int antenna_count) :
        station_id(id),
        antenna(station_id, antenna_count, polar_data),
        power(antenna.power),
        location(COW_LOCATION[station_id])
    {
        /* update channel */
    }
};