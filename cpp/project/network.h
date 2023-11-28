#pragma once
#include <cmath>
#include "coordinates.h"

#define C       3e8 //speed of light
#define M_PIl   3.141592653589793238462643383279502884L /* pi */


namespace network_package
{
    using antennadim = Dimensions<double>;

    std::vector<Placements>             mobile_station_pos;  // SEATING_LOCATION
    std::vector<Placements>             base_station_pos;    // COW_LOCATION

    unsigned                            base_station_count;  // number of base-stations
    unsigned                            mobile_station_count;// number of handsets or mobile stations

    struct Input
    {
        const double   frequency;     // center  frequency in Mz
        const double   bandwidth;     // channel bandwidth in Mz
        const double   SymbRate;      // symbol rate
        const double   BlocksPerSymb; // blocks per symbol
        const double   lambda;        // size of the waveform or wavelength
        const double   height;        // antenna height in meters
        const double   bs_Gtx;        // radiated gain in dBi
        const double   bs_Gtx_W;      // radiated gain in W
        const double   bs_Gtx_min;    // radiated gain in dBi
        const double   bs_Gtx_max;    // radiated gain in dBi
        const double   ms_Grx;        // radiated gain in dBi
        const double   ms_Grx_W;      // radiated gain in W
        const double   bs_scan_min;   // radiated gain in deg
        const double   bs_scan_max;   // radiated gain in deg
        const double   system_noise;  // System noise in dB;
        const unsigned base_stations; // number of base stations
        const double   NF_dB;

        struct Antenna
        {
            std::vector<std::vector<double>>    array_power_wtts;         // antenna array power level [W]
            std::vector<std::vector<double>>    array_scan_angle;    // NODE_SELECTION
            std::vector<unsigned>               antcount_per_base;    // antenna allocation per bs array
            std::vector<double>                 antenna_spacing;     // antenna spacing in an array

            antennadim                          antenna_dim_mtrs;    // antenna dims in [meters]
        };
        Antenna ant;

        Input(
            double   _frequency,
            double   _bandwidth,
            double   _SymbRate,
            double   _BlocksPerSymb,
            double   _lambda,
            double   _height,
            double   _bs_Gtx,
            double   _bs_Gtx_min,
            double   _bs_Gtx_max,
            double   _ms_Grx,
            double   _bs_scan_min,
            double   _bs_scan_max,
            double   _system_noise,
            unsigned _base_stations,
            double   _thermal_system_noise_dB
        ) :
            frequency(_frequency),
            bandwidth(_bandwidth),
            SymbRate(_SymbRate),
            BlocksPerSymb(),
            lambda(getLambda(frequency)),
            height(_height),
            bs_Gtx(_bs_Gtx),
            bs_Gtx_W(log2lin(_bs_Gtx)),
            bs_Gtx_min(_bs_Gtx_min),
            bs_Gtx_max(_bs_Gtx_max),
            ms_Grx(_ms_Grx),
            ms_Grx_W(log2lin(ms_Grx)),
            bs_scan_min(deg2rad(_bs_scan_min)),
            bs_scan_max(deg2rad(_bs_scan_max)),
            system_noise(_system_noise),
            base_stations(_base_stations),
            NF_dB(getThermalSystemNoise(bandwidth, system_noise))
        {
        }
    };

    inline double deg2rad(double deg)
    {
        return deg * M_PIl / 180.0;
    }
    inline double log2lin(double log)
    {
        return pow(10, log / 10);
    }
    inline double dBm2watts(double dBm)
    {
        return log2lin(dBm - 30);
    }

    inline double getLambda(double frequency)
    {
        if (frequency > 0)
            return C / frequency;

        throw std::runtime_error("Divide by zero error from passing 0 frequency in getLambda call");
    }


    /* Get system noise figure in [dBm] with thermal noise by passing B/W in [Mhz], and system noise in [dB] */
    inline double getThermalSystemNoise(const double& bandwidth, const double& system_noise)
    {
        return -174 + round(10 * log10(bandwidth)) + system_noise;
    }


    class AntennaSystem // think of this as multiple antennas each with the following
    {
        const antennadim dims;
        std::vector<double> Gtx;
        unsigned& panel_count;
    public:
        const antennadim size() const { return dims; }
        const std::vector<double>& getgain() const { return Gtx; }
        const int count() const { return panel_count; }
        void setgain(double gain) { Gtx.emplace_back(panel_count * gain); }
        AntennaSystem(double Lx, double Ly, unsigned& _count) : dims(Lx, Ly), panel_count(_count) {}
    };

    // Array of Antenna
    class AAntenna
    {
        friend class Cow;
        const int antenna_id;
        const double ms_Grx_lin;
        AntennaSystem antennas;
        int count;                             // number of panels
        double power;                          // array power
        double theta_c;                        // array orientation
        double spacing;                        // array spacing
        double alpha;                          // antenna array direction
        double lambda;                         // wavelength of the RF signal
        std::vector<double> phee_minus_alpha_list;
        std::vector<double> pathloss_list;
        std::vector<double> hmatrix;  // hmatrix from BS pov


    public:
        /* update the antenna array from updated power and scan angle */
        void update(double alpha)
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

        AAntenna(int id, antennadim dim_meters, double theta, double spacing, unsigned antenna_count, const std::vector<Polar_Coordinates>& polar_sta_data) :
            antenna_id(id),
            ms_Grx_lin(sys->ms_Grx_W),
            antennas(dim_meters.x, dim_meters.y, antenna_count),
            theta_c(theta),
            spacing(spacing),
            lambda(sys->lambda),
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
};
