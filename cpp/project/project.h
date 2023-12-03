#pragma once
#include <unordered_set>
#include <vector>
#include <cassert>
#include "common.h"
#include "network.h"
#include "station.h"
#include "random.h"

/* tests */
#include "test/Init_test_case.h"

using namespace network_package;
using cow_id_distribution = std::uniform_int_distribution<int>;

/* simulation system parameters */
extern network_package::Input* sys;

struct SimulationHelper
{
    std::vector<Cow>& cows;
    unsigned cow_count;
    std::vector<std::vector<double>>& powers;       // TX power level from base station in integer dBm
    std::vector<std::vector<double>>& alphas;       // Antenna array directions in integer rads

    std::vector<unsigned> cow2sta; // Base stations bindings across timeslots

    unsigned timeslot;            // timeslot in the schedule
    unsigned timeslots;           // total timeslots
    unsigned nodes_rem;

    const std::vector<Cow>& getCows() const // this is for performance monitor
    {
        return cows;
    }

    /* gets called to get BS to MS bindings for a single timeslot. Output size: #bs */
    std::vector<unsigned> getBindings()
    {
        std::vector<unsigned> bindings;

        for (unsigned b = cow_count * timeslot; b < cow_count * timeslot + std::min(nodes_rem, cow_count); ++b)
        {
            bindings.emplace_back(cow2sta[b]);
        }

        nodes_rem -= bindings.size();

        return bindings;
    }

    /* set antenna array power level in each timeslot */
    void setPower()
    {
        for (int c = 0; c < cows.size(); ++c)
        {
#ifdef INIT_TEST
            cows[c].setPower(powers[timeslot][c]);
#else
            cows[c].setRandPower();
#endif
        }
    }

    /* set antenna array directivity before starting each simulation */
    void setAlpha()
    {
        for (int c = 0; c < cows.size(); ++c)
        {
#ifdef INIT_TEST
            cows[c].antennaUpdate(alphas[timeslot][c]);
#else
            cows[c].antennaRandUpdate();
#endif
        }
    }

    //bool operator==(const SimulationHelper& other_combo) const
    //{
    //    auto& thispower = this->powers;
    //    auto& powers = other_combo.powers;
    //    for (int i = 0; i < this->powers.size(); ++i)
    //    {
    //        if (thispower[i] != powers[i])
    //            return false;
    //    }
    //
    //    auto& thisalphas = this->alphas;
    //    auto& alphas = other_combo.alphas;
    //    for (int i = 0; i < this->alphas.size(); ++i)
    //    {
    //        if (thisalphas[i] != alphas[i])
    //            return false;
    //    }
    //
    //    return true;
    //}
    //struct combo_hash
    //{
    //    unsigned timeslot;
    //    std::size_t operator()(const SimulationHelper& combo) const {
    //        auto a = std::hash<int>()(combo.powers[timeslot][0]);
    //        auto b = std::hash<int>()(combo.powers[timeslot][1]);
    //        auto c = std::hash<int>()(combo.powers[timeslot][2]);
    //        auto d = std::hash<int>()(combo.alphas[timeslot][0]);
    //        auto e = std::hash<int>()(combo.alphas[timeslot][1]);
    //        auto f = std::hash<int>()(combo.alphas[timeslot][2]);
    //
    //        return a ^ b ^ c ^ d ^ e ^ f;
    //    }
    //    //combo_hash(unsigned _timeslot) : timeslot(_timeslot) {}
    //};

    SimulationHelper(std::vector<Cow>& c) :
        cows(c),
        cow_count(sys->base_stations),
        powers(sys->ant.array_power_wtts),
        alphas(sys->ant.array_scan_angle),
        timeslot(0),
        timeslots(std::ceil(sys->mobile_stations / cow_count)),
        nodes_rem(sys->mobile_stations)
    {
        cow2sta.resize(nodes_rem);

#ifdef INIT_TEST /* use this test case to cross-reference the output and thus, the functionality */
        for (unsigned i = 0; i < cow2sta.size(); ++i)
        {
            cow2sta[i] = init_test_case::Antenna_Array::BS_STAID_SELECTION[i];
        }
#else
        /* randomize the sequence with unique node IDs */
        do {
            cow_id_distribution dist{ 1, cows.size() };
            generate(cow2sta.begin(), cow2sta.end(), [&]() { return randomgen.generate(dist); });

        } while (std::unordered_set<unsigned>(cow2sta.begin(), cow2sta.end()).size() != cow2sta.size());
#endif
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
    telemetry_t() :
        power(0),
        alpha(0),
        pos(),
        cow_id(0),
        timeslot(0),
        sta_id(0),
        sinr(0) {}
};

/* performance logging and analysis latter */
struct PerfMon
{
    const SimulationHelper& params;
    const std::vector<unsigned>& bindings;
    std::vector<telemetry_t> output;

    /* keep track of sinr and note the configurations and bindings in decreasing order of SINR */
    std::map<double, std::vector<telemetry_t>, std::greater<double>> sinrlist;

    void update(double sinr)
    {
        auto& cows = params.getCows();
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


/* antenna power info info */
int get_base_station_power(std::vector<std::vector<double>>& powertable)
{
#ifdef INIT_TEST
    /* change this to reading XML/JSON input file */
    unsigned i = 0, j = 0;
    auto& LUT = init_test_case::Antenna_Array::COW_POWER_LUT;

    powertable.resize(LUT.size());
    for (auto& powerlist : LUT)
    {
        j = 0;
        powertable[i].resize(LUT[0].size());
        for (auto& power : powerlist)
        {
            powertable[i][j++] = dBm2watt(power);
        }
        ++i;
    }

#else
    /* read the JSON file [filename] */*/

#endif
    return 0;
}

/* antenna array direction info */
int get_base_station_scan_alpha(std::vector<std::vector<double>>& alphatable)
{
#ifdef INIT_TEST
    /* change this to reading XML/JSON input file */
    unsigned i = 0, j = 0;
    auto& LUT = init_test_case::Antenna_Array::SCAN_ALPHA_LUT;

    alphatable.resize(LUT.size());
    for (auto& alpha_list : LUT)
    {
        j = 0;
        alphatable[i].resize(LUT[0].size());
        for (auto& alpha : alpha_list)
        {
            alphatable[i][j++] = deg2rad(alpha);
        }
        ++i;
    }

#else
    /* read the JSON file [filename] */*/

#endif
    return 0;
}

int get_parameters(std::string filename)
{
#ifdef INIT_TEST
    auto& frequency     = init_test_case::Network_System::frequency;
    auto& bandwidth     = init_test_case::Network_System::bandwidth;
    auto& SymbRate      = init_test_case::Network_System::SymbRate;
    auto& BlocskPerSymb = init_test_case::Network_System::BlocksPerSymb;
    auto& height        = init_test_case::Network_System::height;
    auto& bs_Ptx_min    = init_test_case::Antenna_Array::MIN_POWER_dBm;
    auto& bs_Ptx_max    = init_test_case::Antenna_Array::MAX_POWER_dBm;
    auto& ms_Grx        = init_test_case::Network_System::ms_Grx;
    auto& bs_scan_min   = init_test_case::Antenna_Array::MIN_SCANANGLE_deg;
    auto& bs_scan_max   = init_test_case::Antenna_Array::MAX_SCANANGLE_deg;
    auto& system_noise  = init_test_case::Network_System::system_noise;
    auto& ms_count      = init_test_case::Network_System::mobile_stations;
    auto& bs_count      = init_test_case::Network_System::base_stations;
#else

#endif
    sys = new Input(
        frequency,
        bandwidth,
        SymbRate,
        BlocskPerSymb,
        height,
        bs_Ptx_min,
        bs_Ptx_max,
        ms_Grx,
        bs_scan_min,
        bs_scan_max,
        system_noise,
        ms_count,
        bs_count
    );

    return 0;
}

template<class U, class V>
inline int getData(const std::vector<U>& input, std::vector<V>& output, std::string filename = "")
{
#ifdef INIT_TEST
    for (auto& data : input)
    {
        output.emplace_back(data);
    }
#else
    /* read the JSON file [filename] */
#endif
    return 0;
}

/* phase array antenna count per base station */
int get_allocation_per_bs(std::vector<unsigned>& counts)
{
    return getData(init_test_case::Network_System::BS_ANTENNA_COUNTS, counts);
}

/* get base station location info */
int get_base_station_cords(std::vector<Placements>& station_position)
{
    return getData(init_test_case::Network_System::COW_LOCATION, station_position);
}

/* get mobile station location info */
int get_mobile_station_cords(std::vector<Placements>& station_position)
{
    return getData(init_test_case::Network_System::NODE_LOCATIONS, station_position);
}

/* this defines the thetaC term in the calculation or where the BS is facing */
int get_base_station_orientation(std::vector<double>& antenna_dir)
{
#ifdef INIT_TEST
    for (auto& data : init_test_case::Network_System::BS_THETAC_DEG)
    {
        antenna_dir.emplace_back(deg2rad(data));
    }
#else
    /* read the JSON file [filename] */
#endif
    return 0;
}

int get_base_station_antspace(std::vector<double>& antenna_spacing)
{
    return getData(init_test_case::Antenna_Array::ANTENNA_SPACING, antenna_spacing);
}

/* get the size of the antenna in [meters] where antennadim is Dimension<double> */
int get_base_station_anntena_dims(antennadim& dim)
{
#ifdef INIT_TEST
    /* change this to reading XML/JSON input file */
    auto& data = init_test_case::Antenna_Array::ANTENNA_DIMS;
    dim.x = data.first;
    dim.y = data.second;
#else
    /* read the JSON file [filename] */*/

#endif
    return 0;
}
