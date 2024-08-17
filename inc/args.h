#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <any>
#include <iostream>
#include <fstream>
#include "common.h"
#include "network.h"
#include <json/json.h>

#define RAPIDJSON_HAS_STDSTRING 1
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

#include "argparse.hpp"

using namespace network_package;

namespace default_t
{
	const std::string test = "tests/init_test.json";

	const double tx_power                       = 0.0;   // default power for all base-stations dBm
	const double scan_angle                     = 0.0;   // default angle for all base-stations degrees
	const double sinr_limit                     = 24.0;  // SINR limit in dB
	const unsigned timeslot                     = 0;     // current timeslot
	const unsigned timeslots                    = 1  ;   // number of timeslots to consider
	const std::vector<double> theta_c_radsdir   = { 60.0, 90.0, 120.0 }; // dBm [min max]
	const std::vector<double> power_range_dBm   = { -30.0, 30.0 };       // dBm [min max]
	const std::vector<double> scan_angle_range  = { -90.0, 90.0 };       // degress [min max]
	const int antenna_element_count             = 5;

	/* Input in deg, facing up in the XY grid */
	const double bs_direction_theta_c_deg       = 90.0;

	/* Input in hertz, default is half-wavelength */
	const double antenna_spacing(const double& frequency)
	{
		return C_SPEED * (1.0 / frequency) * (1.0/ 2.0);
	}
}

template<class T>
struct MultiData_Setup
{
	std::vector<std::vector<T>> data;

	const size_t& size() const
	{
		return data.size();
	}

	const auto& list() const
	{
		return data.front();
	}

	const auto& lut() const
	{
		return data;
	}

	/* check whether lut or not */
	bool is_lut() const
	{
		return data.size() > 1;
	}

	/* default constructor */
	MultiData_Setup() {}
	MultiData_Setup(const std::vector<T>& input) : data(1, input) {}

};

template<class T>
std::string to_string(const std::vector<T>& inputlist, std::string delim = ",")
{
	if (inputlist.size())
	{
		std::string output = std::to_string(inputlist[0]);
		output.reserve(inputlist.size() * 2);
		for (int i = 1; i < inputlist.size(); ++i)
		{
			output += delim + std::to_string(inputlist[i]);
		}
		return output;
	}

	return "";
}

/* join strings with multiple delimeters, e.g. multi_join(3, , */
template<class T>
inline std::string to_string(size_t size, T data, std::string delim = ",")
{
	std::string output_str = std::to_string(data);
	for (int i = 1; i < size; ++i)
	{
		output_str += delim + std::to_string(data);
	}
	return output_str;
}

template<class T>
inline void tokenize(const std::string& input, char token, std::vector<T>& output)
{
	int start_idx = 0;
	int index = 0;
	while ((index = input.find(token, start_idx)) != std::string::npos)
	{
		std::string string_data = input.substr(start_idx, index - start_idx);
		std::istringstream ss(string_data);

		output.resize(output.size() + 1);
		ss >> output.back();

		start_idx = index + 1;
	}

	std::string string_data = input.substr(start_idx);
	std::istringstream ss(string_data);

	output.resize(output.size() + 1);
	ss >> output.back();
}

template<class T>
inline void tokenize(const std::string& input, const std::vector<char>& tokens, std::vector<std::vector<T>>& output)
{
	int index = 0;
	int start_idx = 0;

	while ((index = input.find(tokens[0], start_idx)) != std::string::npos)
	{
		std::string substr = input.substr(start_idx, index - start_idx);
		output.resize(output.size() + 1);
		tokenize(substr, tokens[1], output.back());

		start_idx = index + 1;
	}

	std::string substr = input.substr(start_idx);
	output.resize(output.size() + 1);
	tokenize(substr, tokens[1], output.back());
}

struct Power_Values : MultiData_Setup<double>
{
	bool is_valid_and_convert(const std::vector<double>& reference_bounds)
	{
		for (auto& row : data)
		{
			for (auto& num : row)
			{
				if (!(reference_bounds[0] <= num && num <= reference_bounds[1]))
				{
					return false;
				}
			}

			row = cached::dBm2watt(row);
		}
		return true;
	}

	Power_Values() = default;
	Power_Values(const std::string& double_datalist)
	{
		tokenize(double_datalist, { ' ', ',' }, data);
	}
};

struct Scan_Values : MultiData_Setup<double>
{
	bool is_valid_and_convert(const std::vector<double>& reference_bounds)
	{
		for (auto& row : data)
		{
			for (auto& num : row)
			{
				if (!(reference_bounds[0] <= num && num <= reference_bounds[1]))
				{
					return false;
				}
			}
			row = cached::deg2rad(row);
		}
		return true;
	}

	Scan_Values() = default;
	Scan_Values(const std::string& double_datalist)
	{
		tokenize(double_datalist, { ' ', ',' }, data);
	}
};

struct Binding_Values : MultiData_Setup<unsigned>
{
	std::vector<std::vector<unsigned>> binding_data;

	bool is_valid_and_convert(const Pair<unsigned, unsigned>& bounds)
	{
		auto& num_base_stations = bounds.first;
		auto& num_clients = bounds.second;

		for (auto& ms_list : data)
		{
			if (ms_list.size() == num_base_stations && is_unique(ms_list))
			{
				for (unsigned base_station_idx = 0; base_station_idx < num_base_stations; ++base_station_idx)
				{
					auto ms_idx = ms_list[base_station_idx];

					if (!(0 <= ms_idx && ms_idx <= num_clients))
					{
						return false;
					}
				}

				binding_data.emplace_back(ms_list);
			}
			else
			{
				throw std::invalid_argument("Binding argument has invalid station binding syntax " + to_string(ms_list));
			}
		}

		return !binding_data.empty();
	}

	Binding_Values() = default;
	Binding_Values(const std::string& binding_str)
	{
		tokenize(binding_str, { ' ', ',' }, data);
	}
};

template<class T>
std::string str(std::vector<T> list)
{
	std::string output = list.front();
	for (int i = 1; i < list.size() - 1; ++i)
	{
		output += ',' + std::to_string(list[i]);
	}
	return output;
}

struct Location_Setup
{
	std::vector<Placements> data;
	const size_t& size() const
	{
		return data.size();
	}

	Location_Setup() = default;
	Location_Setup(const std::string& location_data)
	{
		size_t start_idx = 0;
		size_t index = 0;

		std::istringstream ss(location_data);
		std::string temp;

		/* argument <location_data> will be "x1,y1 x2,y2 x3,y3" */
		while (std::getline(ss, temp, ' '))
		{
			auto idx = temp.find(',');

			try
			{
				unsigned x = atol(temp.substr(0, idx).c_str());
				unsigned y = atol(temp.substr(idx + 1).c_str());

				if (!(x >= 0 && y >= 0))
					throw std::invalid_argument("argument is negative: " + location_data);

				data.emplace_back(x, y);
			}
			catch (std::invalid_argument ia)
			{
				std::cout << "Invalid argument or conversion error for location, more specifically " << ia.what() << std::endl;
				exit(-1);
			}
		}
	}
};

/* argparse for C++17 */
struct MyArgs : public argparse::Args
{

	std::string& output_dir                 = kwarg("o,output_dir", "Output directory to put results in").set_default(getcwd());
	std::string& i_json_file                = kwarg("f,file", "Input file used to run the simulation").set_default("");

	std::optional<unsigned>& timeslot       = kwarg("timeslot", "Timeslot for logging reasons").set_default(default_t::timeslot);
	double& frequency                       = kwarg("frequency", "Frequency of the signal system wide");
	double& bandwidth                       = kwarg("bandwidth", "Bandwidth of the singal");
	double& symrate                         = kwarg("symbolrate", "Symbol rate of the data stream");
	double& blockspersym                    = kwarg("blockpersymbol", "Blocks of data per symbol");
	double& antenna_height                  = kwarg("height", "Height of the antenna");
	double& gain_gtrx                       = kwarg("ms_grx", "Antenna gain of the antenna at the mobile or client station in dB");
	double& system_noise                    = kwarg("system_noise", "System noise in the transmitter and receiver");
	unsigned& base_station_count            = kwarg("base_stations", "Number of base stations in the simulation");
	unsigned& mobile_station_count          = kwarg("mobile_stations", "Number of mobile stations in the simulation");
	unsigned& timeslots                     = kwarg("timeslots", "Number of timeslots to carry out the simulation on").set_default(default_t::timeslots);
	double& sinr_limit_dB                   = kwarg("slimit", "SINR limit to consider the configuration as valid to get a good 'slimit' dB signal at the handset").set_default(default_t::sinr_limit);
	std::vector<double>& bs_theta_c         = kwarg("base_station_theta_c", "Direction antennas are facing").set_default(default_t::theta_c_radsdir);
	Location_Setup& base_stations_loc       = kwarg("base_station_location", "Location of base stattions is a list");
	Location_Setup& mobile_stations_loc     = kwarg("mobile_station_location", "Location of mobile stations is a list");
	std::vector<unsigned>& bs_antenna_count = kwarg("base_station_antenna_counts", "Number of panels in the antenna array");
	std::vector<double>& power_range_dBm    = kwarg("base_station_power_range_dBm", "Range of power that base stations will use in dBm").set_default(default_t::power_range_dBm);
	std::vector<double>& scan_angle_range   = kwarg("base_station_scan_angle_range_deg", "Scan angle of the antenna linear array at the base station in degrees").set_default(default_t::scan_angle_range);
	std::vector<double>& antenna_spacing    = kwarg("antenna_spacing", "Antenna spacing between panels").set_default(std::vector<double>());
	std::vector<double>& antenna_dims       = kwarg("antenna_dims", "Antenna dimensions in meters");

	std::optional<Power_Values>& bs_tx_power_dBm   = kwarg("antenna_txpower", "Base station transmit TX power in dBm (list or a lut)");
	std::optional<Scan_Values>& bs_scan_alpha_deg  = kwarg("scan_angle", "Base station scan angle in degrees (list or a lut)");

	Binding_Values& ms_id_selections        = kwarg("ms_selection", "Base station_ID - handset_ID binding. \
											   		Syntax is for 1 timeslot as: 2:3,3:4 or for 2+ timeslots as: 2:3,3:4 4:6,5:7");
	bool& showgui                           = flag("g", "Show gui of the simulation").set_default(false);
	bool& debug                             = flag("q,quiet", "Supress output").set_default(true);

	const std::string& get_input_filename()
	{
		return i_json_file;
	}

	void init()
	{
		if (!i_json_file.empty())
		{
			load_from_json();
		}

		/* make the directory if it does not exist */
		common::mkdir(output_dir);

		assert(antenna_height > 0);
		assert(base_station_count > 0);
		assert(mobile_station_count >= base_station_count);

		/* validate theta_C and ensure that it doesn't change between timelots (1, #bs) size */
		assert(bs_theta_c.size() == base_station_count);

		/* validate the location and antenna panel counts for each station */
		assert(base_stations_loc.data.size() == base_station_count);
		assert(mobile_stations_loc.data.size() == mobile_station_count);
		assert(bs_antenna_count.size() == base_station_count);

		/* assign default values */
		assert(power_range_dBm.size() == 2);
		assert(scan_angle_range.size() == 2);

		assert((antenna_spacing.size() > 0 && antenna_spacing.size() == base_station_count) || antenna_spacing.empty());

		if (antenna_spacing.empty())
			antenna_spacing = std::vector<double>(base_station_count, default_t::antenna_spacing(frequency));

		/* dimension of the antennas for all base stations */
		assert(antenna_dims.size() == 2);

		/* assign default values and ensure values are within range */
		if (bs_tx_power_dBm.has_value())
		{
			assert(bs_tx_power_dBm.value().data.size() == timeslots);
			assert(bs_tx_power_dBm.value().data[0].size() == base_station_count);
			assert(bs_tx_power_dBm.value().is_valid_and_convert(power_range_dBm) == true); // gets converted into watts after this
		}
		else
		{
			bs_tx_power_dBm.value().data = std::vector<std::vector<double>>(timeslots,
												std::vector<double>(base_station_count, dBm2watt(default_t::tx_power)));
		}

		/* assign default values and ensure values are within range */
		if (bs_scan_alpha_deg.has_value())
		{
			assert(bs_scan_alpha_deg.value().data.size() == timeslots);
			assert(bs_scan_alpha_deg.value().data[0].size() == base_station_count);
			assert(bs_scan_alpha_deg.value().is_valid_and_convert(scan_angle_range) == true); // gets converted into rads after this
		}
		else
		{
			bs_scan_alpha_deg.value().data = std::vector<std::vector<double>>(timeslots,
												std::vector<double>(base_station_count, deg2rad(default_t::scan_angle)));
		}

		//assert(ms_id_selections.data.size() * ms_id_selections.data[0].size() == mobile_station_count);
		assert(ms_id_selections.data.size() == timeslots);
		assert(ms_id_selections.data[0].size() == base_station_count);
		assert(ms_id_selections.is_valid_and_convert(Pair<unsigned, unsigned>{base_station_count, mobile_station_count}) == true); // gets converted into rads after this
	}

	void load_from_json()
	{
		using namespace rapidjson;

		std::ifstream ifs(i_json_file);
		if (!ifs.is_open())
		{
			std::cerr << "Cannot open JSON file: " << i_json_file << std::endl;
			exit(EXIT_FAILURE);
		}

		//IStreamWrapper isw(ifs);
		Document doc;
		doc.Parse(i_json_file.c_str());
		if (!doc.IsObject()) {
			std::cerr << "Invalid JSON format: " << i_json_file << std::endl;
			exit(EXIT_FAILURE);
		}

		auto get_double = [&doc](const char* key)
		{
			auto it = doc.FindMember(key);
			return it != doc.MemberEnd() && it->value.IsNumber() ? it->value.GetDouble() : 0.0;
		};

		auto get_string = [&doc](const char* key)
		{
			auto it = doc.FindMember(key);
			return it != doc.MemberEnd() && it->value.IsString() ? std::string(it->value.GetString()) : "";
		};

		auto get_unsigned = [&doc](const char* key)
		{
			auto it = doc.FindMember(key);
			return it != doc.MemberEnd() && it->value.IsUint() ? it->value.GetUint() : 0u;
		};

		auto get_vector_double = [&doc](const char* key)
		{
			std::vector<double> vec;
			auto it = doc.FindMember(key);
			if (doc.HasMember(key) && it->value.IsArray())
			{
				for (auto& v : it->value.GetArray()) {
					if (v.IsNumber()) {
						vec.push_back(v.GetDouble());
					}
				}
			}
			return vec;
		};

		auto get_vector_unsigned = [&doc](const char* key)
		{
			std::vector<unsigned> vec;
			auto it = doc.FindMember(key);
			if (doc.HasMember(key) && it->value.IsArray())
			{
				for (auto& v : it->value.GetArray()) {
					if (v.IsNumber()) {
						vec.push_back(v.GetDouble());
					}
				}
			}
			return vec;
		};

		output_dir = get_string("output_dir");
		timeslot = get_unsigned("timeslot");
		frequency = get_double("frequency");
		bandwidth = get_double("bandwidth");
		symrate = get_double("symbolrate");
		blockspersym = get_double("blockspersym");
		antenna_height = get_double("height");
		gain_gtrx = get_double("ms_grx");
		system_noise = get_double("system_noise");
		base_station_count = get_unsigned("base_stations");
		mobile_station_count = get_unsigned("mobile_stations");
		timeslots = get_unsigned("timeslots");
		sinr_limit_dB = get_double("slimit");
		bs_theta_c = get_vector_double("base_station_theta_c");
		base_stations_loc = Location_Setup(get_string("base_station_location"));
		mobile_stations_loc = Location_Setup(get_string("mobile_station_location"));
		bs_antenna_count = get_vector_unsigned("base_station_antenna_counts");
		power_range_dBm = get_vector_double("base_station_power_range_dBm");
		scan_angle_range = get_vector_double("base_station_scan_angle_range_deg");
		antenna_spacing = get_vector_double("antenna_spacing");
		antenna_dims = get_vector_double("antenna_dims");

		// Optional fields handling
		if (doc.HasMember("antenna_txpower"))
		{
			bs_tx_power_dBm = Power_Values(get_string("antenna_txpower"));
		}

		if (doc.HasMember("scan_angle"))
		{
			bs_scan_alpha_deg = Scan_Values(get_string("scan_angle"));
		}

		if (doc.HasMember("ms_selection"))
		{
			ms_id_selections = Binding_Values(get_string("ms_selection"));
		}
	}
};
