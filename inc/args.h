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

//#include "rapidjson/document.h"
//#include "rapidjson/writer.h"
//#include "rapidjson/stringbuffer.h"

#include "argparse.hpp"

using namespace network_package;

namespace default_t
{
	const std::string test = "tests/init_test.json";

	const double tx_power                      = 0.0;   // default power for all base-stations dBm
	const double scan_angle                    = 0.0;   // default angle for all base-stations degrees
	const double sinr_limit                    = 24.0;  // SINR limit in dB
	const unsigned timeslots                   = 1  ;   // number of timeslots to consider
	const std::vector<double> theta_c_radsdir  = { 60.0, 90.0, 120.0 }; // dBm [min max]
	const std::vector<double> power_range_dBm      = { -30.0, 30.0 };       // dBm [min max]
	const std::vector<double> scan_angle_range = { -90.0, 90.0 };       // degress [min max]
	const int antenna_element_count            = 5;

	/* Input in deg, facing up in the XY grid */
	const double bs_direction_theta_c_deg = 90.0;

	/* Input in hertz, default is half-wavelength */
	const double antenna_spacing(const double& frequency) {
		return C * (1 / frequency) * (1 / 2);
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

struct Binding_Values : MultiData_Setup<std::string>
{
	std::vector<std::vector<unsigned>> binding_data;

	bool is_valid_and_convert(const std::pair<unsigned, unsigned>& reference_bounds)
	{
		for (auto& row : data)
		{
			cached::Cache<unsigned, unsigned> cache;
			binding_data.emplace_back();

			for (auto& num : row)
			{
				auto index = num.find(':');

				if (index != std::string::npos)
				{
					auto& total_bs_stations = reference_bounds.first;
					auto& total_ms_stations = reference_bounds.second;
					unsigned base_station, ms_station;

					std::string base_station_str = num.substr(0, index);
					std::string ms_station_str = num.substr(index + 1);

					base_station = stoul(base_station_str);
					ms_station = stoul(ms_station_str);

					if (0 <= base_station && base_station < total_bs_stations
						&& 0 <= ms_station && ms_station < total_ms_stations)
					{
						if (cache.find(base_station) == cache.end())
						{
							cache.emplace(base_station, ms_station);
							binding_data.back().emplace_back(ms_station);
						}
						else
						{
							throw std::invalid_argument("Station binding " + str(base_station) + ":" + str(ms_station)
								+ " is invalid since base station is used to talk to another handset within the same timeslot");
						}
					}

					return false;
				}
				else
				{
					throw std::invalid_argument("Binding argument does not have a station binding syntax, e.g. bs:ms,bs:ms bs:ms,bs:ms");
				}
			}
		}
		return true;
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


class ProgramConfig
{
protected:
	/* simulation system parameters */

	template<typename T = double> const auto& param(const std::string& key) const
	{
		return std::any_cast<const T&>(data.at(key));
	}

	/* antenna array direction info */
	int get_station_lut(const std::vector<std::vector<unsigned>>& input, std::vector<std::vector<unsigned>>& station_ids)
	{
		station_ids = input;
		return 0;
	}

	template<class U, class V>
	inline int getData(const std::vector<U>& input, std::vector<V>& output)
	{
		for (auto& data : input)
		{
			output.emplace_back(data);
		}
		return 0;
	}

	/* get list of tx power for each base station in a timeslot */
	int get_power_list(const std::vector<double>& input, std::vector<double>& output_power)
	{
		for (unsigned i = 0; i < input.size(); ++i)
		{
			output_power.emplace_back(cached::dBm2watt(input[i]));
		}
		return 0;
	}

	/* get list of scan angle for each base station in a timeslot */
	int get_scan_angle_list(const std::vector<double>& input, std::vector<double>& output_scan_angle)
	{
		for (auto& data : input)
		{
			output_scan_angle.emplace_back(data);
		}

		return 0;
	}

	/* get a list of stations to transmit to in a timeslot */
	int get_station_list(const std::vector<unsigned>& input, std::vector<unsigned>& station_ids)
	{
		return getData(input, station_ids);
	}

	/* phase array antenna count per base station */
	int get_allocation_per_bs(const std::vector<unsigned>& input, std::vector<unsigned>& output_counts)
	{
		return getData(input, output_counts);
	}

	/* get base station location info */
	int get_base_station_cords(const std::vector<Placements>& input, std::vector<Placements>& station_position)
	{
		return getData(input, station_position);
	}

	/* get mobile station location info */
	int get_mobile_station_cords(const std::vector<Placements>& input, std::vector<Placements>& station_position)
	{
		return getData(input, station_position);
	}

	/* this defines the thetaC term in the calculation or where the BS is facing */
	int get_base_station_orientation(const std::vector<double>& input, std::vector<double>& antenna_dir)
	{
		for (auto& data : input)
		{
			antenna_dir.emplace_back(data);
		}

		return 0;
	}

	int get_base_station_antspace(const std::vector<double>& input, std::vector<double>& antenna_spacing)
	{
		return getData(input, antenna_spacing);
	}

	/* get the size of the antenna in [meters] where antennadim is Dimension<double> */
	int get_base_station_anntena_dims(const std::vector<double>& input, antennadim& dim)
	{
		auto& data = input;

		dim.x = data[0];
		dim.y = data[1];

		return 0;
	}
	std::unordered_map<std::string, std::any> data;
};

class Default;
std::filesystem::path operator+(const std::filesystem::path dir, const std::filesystem::path& name)
{
	return dir.string() + '/' + name.string();
}

/* argparse for C++17 */
struct MyArgs : public argparse::Args, public ProgramConfig {

	std::optional<std::string>& output_dir = kwarg("o,output_dir", "Output directory to put results in").set_default(getcwd());

	double& frequency                      = kwarg("frequency", "Frequency of the signal system wide");
	double& bandwidth                      = kwarg("bandwidth", "Bandwidth of the singal");
	double& symrate                        = kwarg("symbolrate", "Symbol rate of the data stream");
	double& blockspersym                   = kwarg("blockpersymbol", "Blocks of data per symbol");
	double& antenna_height                 = kwarg("height", "Height of the antenna");
	double& gain_gtrx                      = kwarg("ms_grx", "Antenna gain of the antenna at the mobile or client station in dB");
	double& system_noise                   = kwarg("system_noise", "System noise in the transmitter and receiver");
	unsigned& base_station_count           = kwarg("base_stations", "Number of base stations in the simulation");
	unsigned& mobile_station_count         = kwarg("mobile_stations", "Number of mobile stations in the simulation");
	unsigned& timeslots                    = kwarg("timeslots", "Number of timeslots to carry out the simulation on").set_default(default_t::timeslots);
	double& sinr_limit_dB                  = kwarg("slimit", "SINR limit to consider the configuration as valid to get a good 'slimit' dB signal at the handset").set_default(default_t::sinr_limit);
	std::vector<double>& bs_theta_c        = kwarg("base_station_theta_c", "Direction antennas are facing").set_default(default_t::theta_c_radsdir);
	Location_Setup& base_stations_loc      = kwarg("base_station_location", "Location of base stattions is a list");
	Location_Setup& mobile_stations_loc    = kwarg("mobile_station_location", "Location of mobile stations is a list");
	std::vector<unsigned>& bs_antenna_count = kwarg("base_station_antenna_counts", "Number of panels in the antenna array");
	std::vector<double>& power_range_dBm   = kwarg("base_station_power_range_dBm", "Range of power that base stations will use in dBm").set_default(default_t::power_range_dBm);
	std::vector<double>& scan_angle_range  = kwarg("base_station_scan_angle_range_deg", "Scan angle of the antenna linear array at the base station in degrees").set_default(default_t::scan_angle_range);
	std::vector<double>& antenna_spacing   = kwarg("antenna_spacing", "Antenna spacing between panels").set_default(std::vector<double>());
	std::vector<double>& antenna_dims      = kwarg("antenna_dims", "Antenna dimensions in meters");

	std::optional<Power_Values>& bs_tx_power_dBm   = kwarg("antenna_txpower", "Base station transmit TX power in dBm (list or a lut)");
	std::optional<Scan_Values>& bs_scan_alpha_deg = kwarg("scan_angle", "Base station scan angle in degrees (list or a lut)");

	Binding_Values& ms_id_selections       = kwarg("ms_selection", "Base station_ID - handset_ID binding. \
													Syntax is for 1 timeslot as: 2:3,3:4 or for 2+ timeslots as: 2:3,3:4 4:6,5:7");
	bool& showgui                          = flag("g", "Show gui of the simulation").set_default(false);
	bool& debug                            = flag("q,quiet", "Supress output").set_default(true);

	void init()
	{
		/* make the directory if it does not exist */
		common::mkdir(output_dir.value());

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

		if (antenna_spacing.size() != base_station_count)
			antenna_spacing = std::vector<double>(base_station_count, default_t::antenna_spacing(frequency));

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
		assert(ms_id_selections.is_valid_and_convert(std::pair<unsigned, unsigned>{base_station_count, mobile_station_count}) == true); // gets converted into rads after this
	}
};


template <typename T> std::string stringify(T x) {
	std::stringstream ss;
	ss << x;
	return ss.str();
}

class JasonHelper : public ProgramConfig
{
	/* helper function to break the cvs delimeters */
	template<class T = double> const auto get_list(std::string input)
	{
		size_t idx2 = input.find(',');
		std::vector<T> out;

		if (idx2 < std::string::npos)
		{
			size_t idx = 0;
			while (idx2 < std::string::npos)
			{
				out.emplace_back(stod(input.substr(idx, idx2 - idx)));
				idx = idx2 + 1;
				idx2 = input.find(',', idx2 + 1);
			}
			out.emplace_back(stod(input.substr(idx, idx2 - idx)));
		}
		return out.size() ? out : std::vector<T>();
	}

	template<class T = double> auto get_table(std::string input)
	{
		size_t idx2 = input.find(',');
		std::vector<std::vector<T>> out;

		if (idx2 < std::string::npos)
		{
			size_t idx = 0;
			while (idx2 < std::string::npos)
			{
				out.emplace_back(stod(input.substr(idx, idx2 - idx)));
				idx = idx2 + 1;
				idx2 = input.find(',', idx2 + 1);
			}
			out.emplace_back(stod(input.substr(idx, idx2 - idx)));
		}
		return out.size() ? out : std::vector<std::vector<T>>();
	}

	template<class T> const auto get_pair(Json::Value::const_iterator input)
	{
		std::vector<Placements> output;

		for (Json::Value::const_iterator ite = (*input).begin(); ite != (*input).end();)
		{
			output.emplace_back();

			output.back().x = stod((*ite).asString());
			++ite;
			output.back().y = stod((*ite).asString());
			++ite;
		}

		return output;
	}

	// RapidJSON
	//struct MyHandler {
	//    const char* type;
	//    std::string data;
	//
	//    MyHandler() : type(), data() {}
	//
	//    bool Null() { type = "Null"; data.clear(); return true; }
	//    bool Bool(bool b) { type = "Bool:"; data = b ? "true" : "false"; return true; }
	//    bool Int(int i) { type = "Int:"; data = stringify(i); return true; }
	//    bool Uint(unsigned u) { type = "Uint:"; data = stringify(u); return true; }
	//    bool Int64(int64_t i) { type = "Int64:"; data = stringify(i); return true; }
	//    bool Uint64(uint64_t u) { type = "Uint64:"; data = stringify(u); return true; }
	//    bool Double(double d) { type = "Double:"; data = stringify(d); return true; }
	//    bool RawNumber(const char* str, SizeType length, bool) { type = "Number:"; data = std::string(str, length); return true; }
	//    bool String(const char* str, SizeType length, bool) { type = "String:"; data = std::string(str, length); return true; }
	//    bool StartObject() { type = "StartObject"; data.clear(); return true; }
	//    bool Key(const char* str, SizeType length, bool) { type = "Key:"; data = std::string(str, length); return true; }
	//    bool EndObject(SizeType memberCount) { type = "EndObject:"; data = stringify(memberCount); return true; }
	//    bool StartArray() { type = "StartArray"; data.clear(); return true; }
	//    bool EndArray(SizeType elementCount) { type = "EndArray:"; data = stringify(elementCount); return true; }
	//private:
	//    MyHandler(const MyHandler& noCopyConstruction);
	//    MyHandler& operator=(const MyHandler& noAssignment);
	//};
public:
	JasonHelper(const std::string& file)
	{
		Json::Value root;
		std::ifstream ifs;

		// Set exceptions to be thrown on failure
		ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);

		std::string filecontents = "";
		try
		{
			ifs.open(file);
			// Read from the file
		}
		catch (const std::ifstream::failure& e)
		{
			std::cerr << "Exception opening/reading file: " << e.what() << '\n';
			std::cerr << "Error code: " << strerror(errno) << '\n';
			exit(errno);
		}

		/* parse the json file and fill in args for the program */
		if (ifs.good())
		{
			Json::CharReaderBuilder builder;
			builder["collectComments"] = false;
			JSONCPP_STRING errs;
			if (!parseFromStream(builder, ifs, &root, &errs))
			{
				std::cout << errs << std::endl;
				throw std::runtime_error("Parser problem: " + str(EXIT_FAILURE));
			}
			//Document d;
			//d.Parse(json);
			//
			//MyHandler handler;
			//StringStream ss(filecontents.c_str());
			//
			//Reader reader;
			//reader.IterativeParseInit();
			//while (!reader.IterativeParseComplete()) {
			//    reader.IterativeParseNext<kParseDefaultFlags>(ss, handler);
			//    cout << handler.type << handler.data << endl;
			//}

			//Json::CharReaderBuilder builder;
			//builder["collectComments"] = false;
			//JSONCPP_STRING errs;
			//if (!parseFromStream(builder, ifs, &root, &errs))
			//{
			//    std::cerr << "JSON parser problem: " << errs << std::endl;
			//    exit(EXIT_FAILURE);
			//}
			//
			/* gather the data */
			for (Json::Value::const_iterator outer = root.begin(); outer != root.end(); outer++)
			{
				if (outer.key().compare("base_station_antenna_counts") == 0)
				{
					for (Json::Value::const_iterator inner = (*outer).begin(); inner != (*outer).end(); inner++)
					{
						data.emplace(inner.key().asString(), get_list<unsigned>((*inner).asString()));
					}
				}
				else if (outer.key().compare("base_station_power_range_dBm") == 0
					|| outer.key().compare("base_station_scan_angle_range_deg") == 0)
				{
					for (Json::Value::const_iterator inner = (*outer).begin(); inner != (*outer).end(); inner++)
					{
						data.emplace(inner.key().asString(), get_list((*inner).asString()));
					}
				}
				else if (outer.key().compare("base_station_theta_c") == 0
					|| outer.key().compare("antenna_spacing") == 0
					|| outer.key().compare("antenna_dims") == 0)
				{
					for (Json::Value::const_iterator inner = (*outer).begin(); inner != (*outer).end(); inner++)
					{
						data.emplace(inner.key().asString(), get_list((*inner).asString()));
					}
				}
				else if (outer.key().compare("base_station_location") == 0
					|| outer.key().compare("mobile_station_location") == 0)
				{
					for (Json::Value::const_iterator inner = (*outer).begin(); inner != (*outer).end(); inner++)
					{
						data.emplace(inner.key().asString(), get_pair<Placements>((*inner).begin()));
					}

				}
				else if (outer.key().compare("base_station_power_dBm_lut") == 0
					|| outer.key().compare("base_station_scan_alpha_deg_lut") == 0)
				{
					for (Json::Value::const_iterator inner = (*outer).begin(); inner != (*outer).end(); inner++)
					{
						for (Json::Value::const_iterator last_ite = (*inner).begin(); last_ite != (*inner).end(); last_ite++)
						{
							data.emplace(last_ite.key().asString(), get_table<double>((*last_ite).asString()));
						}
					}
				}
				else if (outer.key().compare("base_to_mobile_station_id_selection_lut") == 0)
				{
					for (Json::Value::const_iterator inner = (*outer).begin(); inner != (*outer).end(); inner++)
					{
						for (Json::Value::const_iterator last_ite = (*inner).begin(); last_ite != (*inner).end(); last_ite++)
						{
							data.emplace(last_ite.key().asString(), get_table<unsigned>((*last_ite).asString()));
						}
					}
				}
				else
				{
					data.emplace(outer.key().asString(), get_list<double>((*outer).asString()));
				}
			}
		}
		else
		{
			std::cerr << "JasonHelper: json file not read" << std::endl;
		}
	}
};
