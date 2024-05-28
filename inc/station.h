#pragma once
#include "common.h"
#include "network.h"
#include "random.h"

using namespace network_package;

class Cow
{
	const unsigned station_id;

	double min_power_watt;
	double max_power_watt;
	double min_scanangle_rad;
	double max_scanangle_rad;
	AAntenna antenna;
	Placements& location;

	const std::vector<double>& power_list;
	unsigned power_idx;
public:

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

	const unsigned& sid() const
	{
		return station_id;
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
		station_id(id),
		min_power_watt(parameters->bs_Ptx_min),
		max_power_watt(parameters->bs_Ptx_max),
		min_scanangle_rad(parameters->bs_scan_min),
		max_scanangle_rad(parameters->bs_scan_max),
		antenna(id, parameters, polar_data),
		location(base_station_pos[id]),
		power_list(powerlist),
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
				auto power = get_rx_signal_power(b);
				interference += power;
			}
		}

		sinr = signal / (interference + tnf_watt);
	}

	/* get SINR in linear by passing signal and inteference/noise in dB */
	double getSINR() const
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
