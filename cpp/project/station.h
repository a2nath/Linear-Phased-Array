#pragma once
#include "common.h"
#include "network.h"
#include "random.h"

using namespace network_package;

class Cow
{
//#ifdef RAND
//    Random randengine;
//#endif
    const unsigned station_id;

    double min_power_watt;
    double max_power_watt;
    double min_scanangle_rad;
    double max_scanangle_rad;
    AAntenna antenna;
    Placements& location;

    const std::vector<double>& powerlevels;
    unsigned power_idx;
public:
    /* set Gtx power in linear */
    void setRandPower()
    {
        //auto dist = randengine.uni.getdist(min_power_watt, max_power_watt);
        //double output = randengine.generate(dist);
        //antenna.setPower(output);
        //return output;
    }

    /* set Gtx power in linear */
    void setPower(double& input_power)
    {
        antenna.setPower(input_power);
    }

    /* return Gtx power in linear */
    const double& getxPower() const
    {
        return antenna.getPower();
    }

    const double& getAlpha() const
    {
        return antenna.getAlpha();
    }

    const Placements& getPosition() const
    {
        return location;
    }

    const int& sid() const
    {
        return station_id;
    }

    /* update the parameters of the Base Station and get channel state to each station */
    void antennaRandUpdate()
    {
        //auto dist = randengine.uni.getdist(min_scanangle_rad, max_scanangle_rad);
        //int alpha = randengine.generate(dist);
        //return antennaUpdate(alpha);
    }

    /* update the parameters of the Base Station and get channel state to each station */
    void antennaUpdate(double alpha)
    {
        antenna.update(alpha);
    }

    /* set the signal level in Watts (linear): second parameter */
    inline void getSignalLevel(unsigned node_id, double& signal_level_lin) const
    {
        signal_level_lin = antenna.coeff(node_id) * antenna.getPower();
    }


    //antennadim dim_meters, double theta, double spacing, int antenna_count,
    Cow(unsigned& id,
        const Input* parameters,
        const std::vector<double>& powerlist,
        const std::vector<Polar_Coordinates>& polar_data
        ) :
//#ifdef RAND
//        randengine(),
//#endif
        station_id(id),
        min_power_watt(parameters->bs_Ptx_min),
        max_power_watt(parameters->bs_Ptx_max),
        min_scanangle_rad(parameters->bs_scan_min),
        max_scanangle_rad(parameters->bs_scan_max),
        antenna(id, parameters, polar_data),
        location(base_station_pos[id]),
        powerlevels(powerlist),
        power_idx(0)
    {
        /* update channel */
    }
};


class Station
{
    unsigned station_id;
    const std::vector<Cow>& cow_data; // cow that's associated to this station at any given timeslot
    const double& tnf_watt;           // thermal noise floor
    double sinr;
public:

    const unsigned& sid() const
    {
        return station_id;
    }

    /* get signal power in watts */
    inline double get_rx_signal_power(unsigned cow_id)
    {
        double power;
        cow_data[cow_id].getSignalLevel(station_id, power);
        return power;
    }

    void setRX(unsigned bs_id)
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

        sinr = signal / (interference + tnf_watt);
    }

    /* get SINR in linear by passing signal and inteference/noise in dB */
    const double& getSINR() const
    {
        return lin2dB(sinr);
    }

    Station(unsigned id, std::vector<Cow>& cows, const double& noise) :
        station_id(id),
        cow_data(cows),
        tnf_watt(noise),
        sinr(0)
    {
        /* update the stations */
    }

};
