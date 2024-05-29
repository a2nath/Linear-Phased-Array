#pragma once
#include <cmath>
#include <smmintrin.h>
#include <immintrin.h>
#include "coordinates.h"

namespace network_package
{
    using antennadim = Dimensions<double>;
    std::vector<Placements>             mobile_station_pos;  // SEATING_LOCATION
    std::vector<Placements>             base_station_pos;    // COW_LOCATION
    std::vector<std::vector<unsigned>>  station_ids_lut;
    std::vector<unsigned>               station_ids;


    inline double deg2rad(double deg)
    {
        return deg * M_PIl / 180.0;
    }
    inline double log2lin(double log)
    {
        return pow(10, log / 10);
    }
    inline double lin2dB(double lin)
    {
        return 10 * log10(lin);
    }
    inline double dBm2watt(double dBm)
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

    struct Input
    {
        const double   frequency;     // center  frequency in Mz
        const double   bandwidth;     // channel bandwidth in Mz
        const double   SymbRate;      // symbol rate
        const double   BlocksPerSymb; // blocks per symbol
        const double   lambda;        // size of the waveform or wavelength
        const double   height;        // antenna height in meters
        //const double   bs_Gtx;        // radiated gain in dBi
        //const double   bs_Gtx_W;      // radiated gain in W
        const double   bs_Ptx_min;    // radiated power in dBm
        const double   bs_Ptx_max;    // radiated power in dBm
        const double   ms_Grx;        // radiated gain in dBi
        const double   ms_Grx_W;      // radiated gain in W
        const double   bs_scan_min;   // radiated gain in deg
        const double   bs_scan_max;   // radiated gain in deg
        const double   system_noise;  // System noise in dB;
        const unsigned mobile_stations; // number of handsets
        const unsigned base_stations; // number of base stations
        const double   NF_watt;

        struct Antenna
        {
            std::vector<std::vector<double>>    array_power_wtts_lut;  // antenna array power level lookup table [W]
            std::vector<std::vector<double>>    array_scan_angle_lut;  // antenna scan angle in timeslots/BS lookup table [rads]
            std::vector<double>                 array_power_wtts;      // antenna array power level [W]
            std::vector<double>                 array_scan_angle;      // antenna scan angle in timeslots/BS [rads]
            std::vector<unsigned>               antcount_per_base;     // antenna allocation per bs array
            std::vector<double>                 antenna_orientation;   // antenna orientation in its placement in [rads]
            std::vector<double>                 antenna_spacing;       // antenna spacing in an array in [meters]
            antennadim                          antenna_dim_mtrs;      // antenna dims in [meters]
        };
        Antenna ant;

        Input(
            const double& _frequency,
            const double& _bandwidth,
            const double& _SymbRate,
            const double& _BlocksPerSymb,
            const double& _height,
            const double& _bs_Ptx_min,
            const double& _bs_Ptx_max,
            const double& _ms_Grx,
            const double& _bs_scan_min,
            const double& _bs_scan_max,
            const double& _system_noise,
            const unsigned& _mobile_stations,
            const unsigned& _base_stations
        ) :
            frequency(_frequency),
            bandwidth(_bandwidth),
            SymbRate(_SymbRate),
            BlocksPerSymb(_BlocksPerSymb),
            lambda(getLambda(frequency)),
            height(_height),
            //bs_Gtx(_bs_Gtx),
            //bs_Gtx_W(log2lin(_bs_Gtx)),
            bs_Ptx_min(dBm2watt(_bs_Ptx_min)),
            bs_Ptx_max(dBm2watt(_bs_Ptx_max)),
            ms_Grx(_ms_Grx),
            ms_Grx_W(log2lin(ms_Grx)),
            bs_scan_min(deg2rad(_bs_scan_min)),
            bs_scan_max(deg2rad(_bs_scan_max)),
            system_noise(_system_noise),
            mobile_stations(_mobile_stations),
            base_stations(_base_stations),
            NF_watt(log2lin(getThermalSystemNoise(bandwidth, system_noise)))
        {
        }
    };

    class AntennaSystem // think of this as multiple antennas each with the following
    {
        const antennadim& dims;
        std::vector<double> Gtx;
        const unsigned& panel_count;
    public:
        const antennadim& size() const { return dims; }
        const std::vector<double>& getGain() const { return Gtx; }
        const unsigned& count() const { return panel_count; }
        void setGain(double gain) { Gtx.emplace_back(panel_count * gain); }
        AntennaSystem(const antennadim& size, const unsigned& _count) : dims(size), panel_count(_count) {}
    };

    /* Linear Phase Array Antenna */
    class AAntenna
    {
        const double ms_Grx_watt;
        double power;                                // array power
        double alpha;                                // scan angle
        unsigned count;                              // number of panels
        double lambda;                               // wavelength of the RF signal
        double theta_c;                              // antenna array direction
        double spacing;                              // array spacing
        AntennaSystem antennas;
        std::vector<double> phee_minus_alpha_list;
        std::vector<double> pathloss_list;
        std::vector<double> hmatrix;                 // hmatrix from BS pov


    public:

        double coeff(unsigned& id) const
        {
            return hmatrix[id];
        }

        void setPanelCount(unsigned panels)
        {
            count = panels;
        }

        /* set power in watts */
        void setPower(double power_watts)
        {
            power = power_watts;
        }

        /* get power in watts */
        const double& getPower() const
        {
            return power;
        }

        /* get alpha in rads */
        const double& getAlpha() const
        {
            return alpha;
        }

        /* get the physical antenna direction */
        const double& beamdir() const
        {
            return theta_c;
        }


        /* update the antenna array from updated power and scan angle */
        void update(double new_alpha)
        {
            auto antenna_count = antennas.count();
            auto& antgain_factor = antennas.getGain();

            /* update the antenna gain Gtx */
            for (unsigned idx = 0; idx < phee_minus_alpha_list.size(); ++idx)
            {
                double phee = (phee_minus_alpha_list[idx] + new_alpha) / 2;

                double sin_term = antenna_count * sin(phee);
                double gain_factor_antenna_system = antgain_factor[idx]; // xN antennas already

                if (sin_term != 0)
                    gain_factor_antenna_system *= pow(sin(antenna_count * phee) / sin_term, 2);

                /* update the channel matrix */
                hmatrix[idx] = ms_Grx_watt * gain_factor_antenna_system / pathloss_list[idx];
            }

            alpha = new_alpha;
        }

        const std::vector<double>& getMatrix() const
        {
            return hmatrix;
        }

        AAntenna(
				unsigned& id,
				const Input* parameters,
				const std::vector<Polar_Coordinates>& polar_sta_data)
			:
			power(0),
			alpha(0),
            ms_Grx_watt(parameters->ms_Grx_W),
            count(parameters->ant.antcount_per_base[id]),
            lambda(parameters->lambda),
            theta_c(parameters->ant.antenna_orientation[id]),
            spacing(parameters->ant.antenna_spacing[id]),
            antennas(parameters->ant.antenna_dim_mtrs, parameters->ant.antcount_per_base[id]),
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

                antennas.setGain(singleant_gain);
                phee_minus_alpha_list.emplace_back(phee_temp * sin(theta_minus_thetaC));
                pathloss_list.emplace_back(pow(pl_temp_meters * station_data.hype, 2));
            }
        }
    };
};
