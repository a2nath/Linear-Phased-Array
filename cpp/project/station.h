#pragma once
#include "common.h"
#include "network.h"

#ifdef RAND
#include "random.h"
#endif


class Cow
{
#ifdef RAND
    Random randengine;
#endif
    const int station_id;

    double min_power_dbm;
    double max_power_dbm;
    double min_scanangle_deg;
    double max_scanangle_deg;

    AAntenna antenna;
    double& power;
    Placements& location;


public:
    static std::unordered_map<int, double> power_list_dBm2W;


    /* set Gtx power in linear */
    const int set_rand_power()
    {
        auto dist = randomgen.uni.getdist(min_power_dbm, max_power_dbm);
        int power = randomgen.generate(dist);
        return setpower(power);
    }

    /* set Gtx power in linear */
    const int setpower(int input_power)
    {
        power = power_list_dBm2W[input_power];
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
    const double antenna_rand_update()
    {
        auto dist = randomgen.uni.getdist(min_scanangle_deg, max_scanangle_deg);
        int alpha = randomgen.generate(dist);
        return antenna_update(alpha);
    }

    /* update the parameters of the Base Station and get channel state to each station */
    const double antenna_update(double alpha)
    {
        antenna.update(alpha);
        return alpha;
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


    //antennadim dim_meters, double theta, double spacing, int antenna_count,
    Cow(unsigned id,
        antennadim dim_meters,
        double theta_deg,
        double spacing,
        unsigned antenna_count,
        const std::vector<Polar_Coordinates>& polar_data,

        double minpower_dbm,
        double maxpower_dbm,
        double minscanangle_deg,
        double maxscanangle_deg
    ) :
#ifdef RAND
        randengine(engine),
#endif
        station_id(id),
        antenna(id, dim_meters, theta_deg, spacing, antenna_count, polar_data),
        power(antenna.power),
        location(base_station_pos[id]),
        min_power_dbm(minpower_dbm),
        max_power_dbm(maxpower_dbm),
        min_scanangle_deg(minscanangle_deg),
        max_scanangle_deg(maxscanangle_deg)
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
    double noise;
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

        sinr = signal / (interference + dBm2watts(noise));
    }

    /* get SINR in linear by passing signal and inteference/noise in Watts */
    const double& getSINR() const
    {
        return sinr;
    }

    Station(unsigned id, std::vector<Cow>& cows) :
        station_id(id),
        cow_data(cows),
        sinr(0),
        noise(sys->NF_dB)
    {
        /* update the stations */
    }

};