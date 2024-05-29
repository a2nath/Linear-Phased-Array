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
        double power;                                // array power
        double alpha;                                // scan angle
        const double ms_Grx_watt;
        unsigned count;                              // number of panels
        double lambda;                               // wavelength of the RF signal
        double spacing;                              // array spacing
        double theta_c;                              // antenna array direction
        antennadim antenna_dims;

        struct Calculations
        {
            /* Gtx gain for the whole grid */
            std::vector<double> gain_RX_grid;
            std::vector<double> phee_minus_alpha_list;
            std::vector<double> pathloss_list;
            std::vector<double> hmatrix;                 // hmatrix from BS pov

            void resize(const size_t& size)
            {
                gain_RX_grid.resize(size);
                phee_minus_alpha_list.resize(size);
                pathloss_list.resize(size);
                hmatrix.resize(size);
            }
        };

        /* Calculations needed to create/update the H-matrix coefficient table */
        Calculations simulation, graphic;

    public:
        /* accepts a RX station ID with respect to THIS station
           and returns the coefficient
        */
        double coeff(const unsigned& rx_sta) const
        {
            return simulation.hmatrix[rx_sta];
        }

        void setPanelCount(unsigned panels)
        {
            count = panels;
        }

        /* set power in watts */
        void set_power(const double& power_watts)
        {
            power = power_watts;
        }

        /* get power in watts */
        const double& get_power() const
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
        inline void update(const double& new_alpha, Calculations& calculations)
        {
            if (new_alpha != alpha)
			{
                /* update the antenna gain Gtx */
                for (unsigned idx = 0; idx < calculations.phee_minus_alpha_list.size(); ++idx)
                {
                    double phee = (calculations.phee_minus_alpha_list[idx] + new_alpha) / 2;

					double sin_term = count * sin(phee);
					double gain_factor_antenna_system = calculations.gain_RX_grid[idx]; // xN antennas already

                    if (sin_term != 0)
                    {
                        gain_factor_antenna_system *= pow(sin(count * phee) / sin_term, 2);
                    }

                    /* update the channel matrix */
                    calculations.hmatrix[idx] = ms_Grx_watt * gain_factor_antenna_system / calculations.pathloss_list[idx];
                }

                /* update the scan angle of the antenna array */
                alpha = new_alpha;
			}
        }

        /* for bare-minimum numerical calculations needed at the mobile_stations only */
        void graphics_update(const double& new_alpha)
        {
            update(new_alpha, graphic);
        }

        /* for bare-minimum numerical calculations needed at the mobile_stations only */
        void numerical_update(const double& new_alpha)
        {
            update(new_alpha, simulation);
        }

        /* re-calc the signal outs to handsets only (before calling update!) */
        inline void init(const std::vector<Polar_Coordinates>& polar_data, Calculations& calculations)
        {
            const double& pioverlambda = M_PIl / lambda;
            const double& phee_temp = 2 * spacing * pioverlambda;
            const double& pl_temp_meters = 4 * pioverlambda;
            const double& antenna_dim_factor = 10 * antenna_dims.x * antenna_dims.y / pow(lambda, 2);
            double m_factor = antenna_dims.x * pioverlambda;

            for (size_t idx = 0; idx < polar_data.size(); ++idx)
            {
                auto& cell_polar_data = polar_data[idx];

                double theta_minus_thetaC = cell_polar_data.theta - theta_c;
                double m = m_factor * sin(theta_minus_thetaC);
                double singleant_gain = antenna_dim_factor * pow((1 + cos(theta_minus_thetaC)) / 2, 2);

                if (m != 0)
                {
                    singleant_gain *= pow(sin(m) / m, 2);
                }

                calculations.phee_minus_alpha_list[idx] = phee_temp * sin(theta_minus_thetaC);
                calculations.pathloss_list[idx]         = pow(pl_temp_meters * cell_polar_data.hype, 2);
                calculations.gain_RX_grid[idx]          = singleant_gain * count;
            }
        }

        /* for GUI simulation in the whole grid */
        void graphics_init(const std::vector<Polar_Coordinates>& polar_data)
        {
            graphic.resize(polar_data.size());
            init(polar_data, graphic);
        }

        /* for bare-minimum numerical calculations needed at the mobile_stations only */
        void numerical_init(const std::vector<Polar_Coordinates>& polar_data)
        {
            simulation.resize(polar_data.size());
            init(polar_data, simulation);
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
            antenna_dims(parameters->ant.antenna_dim_mtrs)
        {
        }
    };
};
