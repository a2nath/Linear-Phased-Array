#pragma once
#include "common.h"
#include <unordered_set>
#include <vector>
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
const int MIN_SCANANGLE_deg = -90;
const int MAX_SCANANGLE_deg = +90;

const std::pair<double, double> ANTENNA_DIMS = { 20e-2, 40e-2 }; // meters

const double system_noise = 5; // System noise in dB;
const double thermal_n_db = -174 + round(10 * log10(bandwidth)) + system_noise;  // Thermal noise floor in dBm

using Placements = Coordinates<unsigned>;

using cow_id_distribution = std::uniform_int_distribution<int>;
extern Random randomgen;


const std::vector<std::vector<double>> SCAN_ALPHA_LUT = {        // degrees
    {+10.0, +65.0, +35.0},
    {-10.0, -65.0, -35.0},
    {-75.0, -40.0, +75.0},
    {+40.0, +65.0, -40.0},
    {+40.0, -65.0, -40.0}
};

const std::vector<Placements> SEATING_LOCATION = // seating of all mobile stations in meters
{
    {250, 500}, {350, 500}, {450, 500}, {550, 500}, {650, 500}, {750, 500},
    {250, 600}, {350, 600}, {450, 600}, {550, 600}, {650, 600}, {750, 600},
                      {450, 325}, {500, 325}, {550, 325}
};

const std::vector<Placements> COW_LOCATION =
{
    {250, 200},
    {500, 180},
    {750, 200}
};

/* Setup the environment */
std::vector<unsigned> NODE_SELECTION = {
    2, 13, 11,
    8, 15, 05,
    1, 14, 06,
    7, 04, 10,
    9, 03, 12
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
            double gain_factor_antenna_system = antgain_factor[idx]; // xN antennas already

            if (sin_term != 0)
                gain_factor_antenna_system *= pow(sin(antenna_count * phee) / sin_term, 2);

            /* update the channel matrix */
            hmatrix[idx] = ms_Grx_lin * gain_factor_antenna_system / pathloss_list[idx];
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
    Placements location;
public:
    static std::unordered_map<int, double> dBm2watts_lut;
    //static unordered_map<double, double> deg2rad_lut;

    /* set Gtx power in linear */
    const int setpower(int input_power = rand(MIN_POWER_dBm, MAX_POWER_dBm))
    {
        power = dBm2watts_lut[input_power];
        return input_power;
    }

    /* return Gtx power in linear */
    const double& getxpower() const
    {
        return power;
    }
    const int& sid() const
    {
        return station_id;
    }

    /* update the parameters of the Base Station and get channel state to each station */
    const int antenna_update(int alpha = rand(MIN_SCANANGLE_deg, MAX_SCANANGLE_deg))
    {
        antenna.update(deg2rad(alpha));
        return alpha;
        //return antenna.getmatrix();
    }

    const double& get_alpha() const
    {
        return antenna.alpha;
    }

    const Placements& get_position() const
    {
        return location;
    }

    /* set the signal level in Watts (linear): second parameter */
    inline void get_signal_level(unsigned node_id, double& signal_level_lin) const
    {
        signal_level_lin = antenna.hmatrix[node_id] * power;
    }

    Cow(unsigned id, const std::vector<Polar_Coordinates>& polar_data, int antenna_count) :
        station_id(id),
        antenna(station_id, antenna_count, polar_data),
        power(antenna.power),
        location(COW_LOCATION[station_id])
    {
        /* update channel */
    }
};

class MainCow : public Cow
{

};

class SideCow : public Cow
{

};

class Station
{
    int station_id;
    const std::vector<Cow>& cow_data; // cow that's associated to this station at any given timeslot
    double sinr;
public:

    const unsigned& sid() const
    {
        return station_id;
    }

    inline double get_rx_signal_power(unsigned cow_id)
    {
        double power;
        cow_data[cow_id].get_signal_level(station_id, power);
        return power;
    }

    double setRX(unsigned bs_id)
    {
        double signal = get_rx_signal_power(bs_id);

        double interference = 0;

        for (unsigned b = 0; b < cow_data.size(); ++b)
        {
            if (cow_data[b].sid() != bs_id)
            {
                interference += get_rx_signal_power(b);
            }
        }

        sinr = signal / (interference + dBm2watts(thermal_n_db));
    }

    /* get SINR in linear by passing signal and inteference/noise in Watts */
    const double& getSINR() const
    {
        return sinr;
    }

    Station(unsigned id, std::vector<Cow>& cows) :
        station_id(id),
        cow_data(cows),
        sinr(0)
    {
        /* update the stations */
    }

};

struct telemetry_t
{
    int power;
    double alpha;
    Placements pos;
    unsigned cow_id;
    unsigned timeslot;
    int sta_id;
    double sinr;
    telemetry_t(const Cow& cow, const unsigned& timeslot, const unsigned& sta, double& _sinr)
        :
        power(cow.getxpower()),
        alpha(cow.get_alpha()),
        pos(cow.get_position()),
        cow_id(cow.sid()),
        timeslot(timeslot),
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

struct SimulationHelper
{
    std::vector<Cow>& cows;
    unsigned cow_count;

    std::vector<int> powers;       // TX power level from base station in integer dBm
    std::vector<int> alphas;       // Antenna array directions in integer rads

    std::vector<unsigned> cow2sta; // Base stations bindings across timeslots

    unsigned timeslot;            // timeslot in the schedule
    unsigned nodes_rem;

    const std::vector<Cow>& getcows() const // this is for performance monitor
    {
        return cows;
    }

    /* gets called to get BS to MS bindings for a single timeslot. Output size: #bs */
    std::vector<unsigned> getbindings()
    {
        std::vector<unsigned> bindings;

        for (unsigned b = cow_count * timeslot; b < cow_count * timeslot + min(nodes_rem, cow_count); ++b)
        {
            bindings.emplace_back(cow2sta[b]);
        }

        ++timeslot;
        nodes_rem -= bindings.size();

        return bindings;
    }

    /* set antenna array power level in each timeslot */
    void setpower()
    {
        for (int c = 0; c < cows.size(); ++c)
        {
            double power = cows[c].setpower();
            powers[c] = power;
        }
    }

    /* set antenna array directivity before starting each simulation */
    void setalpha()
    {
        for (int c = 0; c < cows.size(); ++c)
        {
            double alpha = cows[c].antenna_update();
            alphas[c] = alpha;
        }
    }

    bool operator==(const SimulationHelper& other_combo) const
    {
        auto& thispower = this->powers;
        auto& powers = other_combo.powers;
        for (int i = 0; i < this->powers.size(); ++i)
        {
            if (thispower[i] != powers[i])
                return false;
        }

        auto& thisalphas = this->alphas;
        auto& alphas = other_combo.alphas;
        for (int i = 0; i < this->alphas.size(); ++i)
        {
            if (thisalphas[i] != alphas[i])
                return false;
        }

        return true;
    }
    struct combo_hash
    {
        std::size_t operator()(const SimulationHelper& combo) const {
            auto a = std::hash<int>()(combo.powers[0]);
            auto b = std::hash<int>()(combo.powers[1]);
            auto c = std::hash<int>()(combo.powers[2]);
            auto d = std::hash<int>()(combo.alphas[0]);
            auto e = std::hash<int>()(combo.alphas[1]);
            auto f = std::hash<int>()(combo.alphas[2]);

            return a ^ b ^ c ^ d ^ e ^ f;
        }
    };

    SimulationHelper(std::vector<Cow>& c) :
        cows(c),
        cow_count(c.size()),
        powers(c.size()),
        alphas(c.size()),
        cow2sta(NODE_SELECTION.size()),
        timeslot(0),
        nodes_rem(NODE_SELECTION.size())
    {
        /* randomize the sequence withi unique numbers */
        do {
            cow_id_distribution dist{ 1, cows.size() };
            generate(cow2sta.begin(), cow2sta.end(), [&]() { return randomgen.generate(dist); });

        } while (std::unordered_set<unsigned>(cow2sta.begin(), cow2sta.end()).size() != cow2sta.size());

        timeslot = std::ceil(NODE_SELECTION.size() / cow_count);
    }
};

/* performance logging and analysis latter */
struct PerfMon
{
    const SimulationHelper& params;
    const std::vector<unsigned>& bindings;

    /* keep track of sinr and note the configurations and bindings in decreasing order of SINR */
    std::map<double, std::vector<telemetry_t>, std::greater<double>> sinrlist;

    void update(double sinr)
    {
        auto& cows = params.getcows();
        for (unsigned c = 0; c < params.cow_count; ++c)
        {
            sinrlist[sinr].emplace_back(telemetry_t(cows[c], params.timeslot, bindings[c++], sinr));
        }
    }

    PerfMon(const SimulationHelper& parameters, std::vector<unsigned>& cowbindings) :
        params(parameters),
        bindings(cowbindings)
    {}
};