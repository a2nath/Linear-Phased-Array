#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <any>
#include <iostream>
#include <fstream>
#include <optional>
#include <json/json.h>
#include "common.h"
#include "json_config.h"
#ifdef __CUDACC__
#include "network.cuh"
#else
#include "network.h"
#endif

#include "argparse.hpp"

using namespace network_package;

template<class T>
struct MultiData_Setup
{
	std::vector<std::vector<T>> data;

	const size_t size() const
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
	std::istringstream stream(input);
	std::string item;

	while (std::getline(stream, item, token))
	{
		std::istringstream value_stream(item);
		T value{};
		value_stream >> std::ws;

		if ((std::is_unsigned<T>::value && value_stream.peek() == '-') || !(value_stream >> value))
		{
			throw std::invalid_argument("Invalid numeric value: " + item);
		}

		value_stream >> std::ws;
		if (!value_stream.eof())
		{
			throw std::invalid_argument("Invalid numeric value: " + item);
		}

		output.push_back(value);
	}

	if (input.empty() || input.back() == token)
	{
		throw std::invalid_argument("Missing numeric value in: " + input);
	}
}

/*
-----------------------
Inputs such as
3,4,5 6,7,8 9,2,1
means
timstot1: TX1 --> 3, TX2 --> 4, TX3 --> 5
timstot2: TX1 --> 6, TX2 --> 7, TX3 --> 8
timstot3: TX1 --> 9, TX2 --> 2, TX3 --> 1
-----------------------
*/
template<class T>
inline void tokenize(const std::string& input, const std::vector<char>& tokens, std::vector<std::vector<T>>& output)
{
	std::istringstream stream(input);
	std::string row;

	while (stream >> row) // " " space delim is for between time slots, "," comma is for different TXs
	{
		output.emplace_back();
		tokenize(row, tokens[1], output.back());
	}
}

/* 2D table of double to store and process power (from dBm to Watt) */
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

			row = rf_math::dBm2watt(row);
		}
		return true;
	}

	Power_Values() = default;
	Power_Values(const std::string& double_datalist)
	{
		tokenize(double_datalist, { ' ', ',' }, data);
	}
};

/* 2D table of double to store and process scan angles (from deg to rad) */
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
			row = rf_math::deg2rad(row);
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
	/* Validates the size of inputs, checks for bounds, and converts the values to 0-indexing */
	bool is_valid_and_convert(unsigned num_base_stations, unsigned num_clients)
	{
		for (auto& rx_tx_mapping_list : data)
		{	// num_base_stations = cols
			if (rx_tx_mapping_list.size() == num_base_stations && is_unique(rx_tx_mapping_list))
			{
				for (unsigned base_station_idx = 0; base_station_idx < num_base_stations; ++base_station_idx)
				{
					unsigned& ms_idx = rx_tx_mapping_list[base_station_idx];

					if (!(0 < ms_idx && ms_idx <= num_clients))
					{
						return false;
					}

					ms_idx--; // need to use 0-indexing, CRITICAL.
				}
			}
			else
			{
				throw std::invalid_argument("Binding argument has invalid station binding syntax " + to_string(rx_tx_mapping_list));
			}
		}

		return !data.empty();
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
	const size_t size() const
	{
		return data.size();
	}

	Location_Setup() = default;
	Location_Setup(const std::string& location_data)
	{
		std::istringstream ss(location_data);
		std::string temp;

		/* argument <location_data> will be "x1,y1 x2,y2 x3,y3" */
		while (ss >> temp)
		{
			unsigned_v coordinates;
			tokenize(temp, ',', coordinates);
			if (coordinates.size() != 2)
			{
				throw std::invalid_argument("Location must be an x,y pair: " + temp);
			}

			data.emplace_back(coordinates[0], coordinates[1]);
		}
	}
};


namespace defaults_t
{
	const std::string test = "tests/init_test.json";
	/* Input in hertz, default is half-wavelength */

	double antenna_spacing(const double& frequency)
	{
		return C_SPEED * (1.0 / frequency) * (1.0 / 2.0);
	}

	//const std::string output_dir = getcwd();
	//const std::string& i_json_file = ""
	constexpr unsigned timeslot = 0;
	constexpr unsigned timeslots = 1;
	constexpr unsigned base_station_count = 3;
	constexpr unsigned mobile_station_count = 3;
	constexpr double frequency = 1.9e9;
	constexpr double bandwidth = 20e6;
	constexpr double symrate = 3.84e6;
	constexpr double blockspersym = 768;
	constexpr double bs_antenna_height = 200;
	constexpr double mobile_gain_dB = -2;
	constexpr double system_noise_dB = 5;
	constexpr double ms_sinr_limit_dB = 24;
	const unsigned_v field_size = { 1000, 1000 };  // size of the render area not including controls etc.
	const double_v bs_theta_c_rads_dir = { 75.0, 90.0, 105.0 }; // Facing directions in degrees. East = 0 degrees
	const Location_Setup tx_loc("250,199 500,180 750,199");
	// Source clients 7, 4, 10, compacted into client IDs 1, 2, 3.
	const Location_Setup rx_loc("250,600 550,500 550,600");
	const unsigned_v bs_antenna_count = { 5, 3, 5 };
	const Location_Setup tx_location_range("0,0 1000,200");
	const Location_Setup rx_location_range("0,200 1000,1000");
	const double_v bs_theta_c_range_deg = { 0.0, 360.0 };
	const unsigned_v bs_antenna_count_limit = { 6, 6, 6 };
	const double_v bs_power_range_dBm = { -30.0, 30.0 };
	const double_v bs_scan_angle_range_deg = { -90.0, 90.0 };
	const double_v bs_antenna_spacing_m = { 0.3, 0.4, 0.3 };
	const double_v bs_antenna_dims = { 0.2, 0.4 };
	const Power_Values bs_tx_power_dBm = Power_Values("30.0,29.0,28.0");
	const Scan_Values bs_scan_alpha_deg = Scan_Values("40.0,65.0,-40.0");
	const Binding_Values ms_id_selections = Binding_Values("1,2,3");
	constexpr bool gui = false;
	constexpr bool plot = true;
}

/* argparse for C++17 */
struct MyArgs : public argparse::Args
{
	std::string& output_dir            = kwarg("o,output_dir", "Output directory to put results in").set_default(getcwd());
	std::string& i_json_file           = kwarg("f,file", "Input file used to run the simulation").set_default("");

	unsigned& timeslot                 = kwarg("timeslot", "Timeslot for logging reasons").set_default(defaults_t::timeslot);
	unsigned& timeslots                = kwarg("timeslots", "Number of timeslots to carry out the simulation on").set_default(defaults_t::timeslots);
	unsigned& base_station_count       = kwarg("base_stations", "Number of base stations in the simulation").set_default(defaults_t::base_station_count);
	unsigned& mobile_station_count     = kwarg("mobile_stations", "Number of mobile stations in the simulation").set_default(defaults_t::mobile_station_count);
	double& frequency                  = kwarg("frequency", "Frequency of the signal system wide").set_default(defaults_t::frequency);
	double& bandwidth                  = kwarg("bandwidth", "Bandwidth of the signal").set_default(defaults_t::bandwidth);
	double& symrate                    = kwarg("symbolrate", "Symbol rate of the data stream").set_default(defaults_t::symrate);
	double& blockspersym               = kwarg("blockpersymbol", "Blocks of data per symbol").set_default(defaults_t::blockspersym);
	double& bs_antenna_height          = kwarg("height", "Height of the antenna").set_default(defaults_t::bs_antenna_height);
	double& gain_gtrx                  = kwarg("ms_grx", "Antenna gain of the antenna at the mobile or client station in dB").set_default(defaults_t::mobile_gain_dB);
	double& system_noise               = kwarg("system_noise", "System noise in the transmitter and receiver").set_default(defaults_t::system_noise_dB);
	double& ms_sinr_limit_dB           = kwarg("slimit", "SINR limit to consider the configuration as valid to get a good 'slimit' dB signal at the handset").set_default(defaults_t::ms_sinr_limit_dB);
	unsigned_v& field_size             = kwarg("area", "Size of the field to render GUI. Must be compiled with -DGRAPHICS").set_default(defaults_t::field_size);
	double_v& bs_theta_c               = kwarg("base_station_theta_c", "Direction antennas are facing").set_default(defaults_t::bs_theta_c_rads_dir);
	Location_Setup& tx_loc             = kwarg("base_station_location", "Locations of base stations").set_default(defaults_t::tx_loc);
	Location_Setup& rx_loc             = kwarg("mobile_station_location", "Locations of mobile stations").set_default(defaults_t::rx_loc);
	unsigned_v& bs_antenna_count       = kwarg("base_station_antenna_counts", "Number of panels in the antenna array").set_default(defaults_t::bs_antenna_count);

	// TX and RX ranges are independent rectangles in absolute map coordinates.
	Location_Setup& tx_location_range  = kwarg("base_station_location_range_meters", "TX rectangle: xmin,ymin xmax,ymax").set_default(defaults_t::tx_location_range);
	Location_Setup& rx_location_range  = kwarg("mobile_station_location_range_meters", "RX rectangle: xmin,ymin xmax,ymax").set_default(defaults_t::rx_location_range);
	double_v& bs_theta_c_range_deg     = kwarg("base_station_theta_c_range_deg", "Allowed antenna facing directions in degrees").set_default(defaults_t::bs_theta_c_range_deg);
	unsigned_v& bs_antenna_count_limit = kwarg("base_station_antenna_counts_limit", "Maximum panel count per transmitter").set_default(defaults_t::bs_antenna_count_limit);
	double_v& bs_power_range_dBm       = kwarg("base_station_power_range_dBm", "Range of power that base stations will use in dBm").set_default(defaults_t::bs_power_range_dBm);
	double_v& bs_scan_angle_range_deg  = kwarg("base_station_scan_angle_range_deg", "Scan angle of the antenna linear array at the base station in degrees").set_default(defaults_t::bs_scan_angle_range_deg);
	double_v& bs_antenna_spacing_m     = kwarg("antenna_spacing", "Antenna spacing between panels").set_default(defaults_t::bs_antenna_spacing_m);
	double_v& bs_antenna_dims          = kwarg("antenna_dims", "Antenna dimensions in meters").set_default(defaults_t::bs_antenna_dims);
	Power_Values& bs_tx_power_dBm      = kwarg("antenna_txpower", "Base station transmit TX power in dBm (list or a lut)").set_default(defaults_t::bs_tx_power_dBm);
	Scan_Values& bs_scan_alpha_deg     = kwarg("scan_angle", "Base station scan angle in degrees (list or a lut)").set_default(defaults_t::bs_scan_alpha_deg);
	Binding_Values& ms_id_selections   = kwarg("ms_selection", "1-based client IDs, one per transmitter column; rows separated by spaces: 2,3,4 5,6,7").set_default(defaults_t::ms_id_selections);
	bool& gui                          = flag("gui", "Real-time render after printing simulation results").set_default(false);
	bool& plot                         = flag("plot", "Capture and download the heatmaps").set_default(true);

	void parse(int argc, const char* const* argv, const bool& raise_on_error = false)
	{
		// With JSON, only the input filename and output directory are CLI options.
		// Filtering before argparse also ignores malformed CLI configuration overrides.
		bool has_json_file = false;
		std::vector<const char*> file_arguments{argv[0]};

		for (int index = 1; index < argc; ++index)
		{
			const std::string option = argv[index];
			const std::string key = option.substr(0, option.find('='));
			const bool is_file = key == "--file" || key == "-f";
			const bool is_output = key == "--output_dir" || key == "-o";
			if (!is_file && !is_output)
			{
				continue;
			}

			has_json_file = has_json_file || is_file;
			file_arguments.push_back(argv[index]);
			if (option.find('=') == std::string::npos)
			{
				if (index + 1 >= argc || std::string(argv[index + 1]).rfind('-', 0) == 0)
				{
					throw std::invalid_argument("Argument missing a value: " + key);
				}
				file_arguments.push_back(argv[++index]);
			}
		}

		if (has_json_file)
		{
			argparse::Args::parse(static_cast<int>(file_arguments.size()), file_arguments.data(), raise_on_error);
		}
		else
		{
			argparse::Args::parse(argc, argv, raise_on_error);
		}
	}

	void print() const
	{
		if (i_json_file.empty())
		{
			argparse::Args::print();
		}
		else
		{
			// argparse records CLI text, which does not describe a JSON configuration.
			std::cout << "Input JSON (before unit conversion): " << i_json_file << '\n';
			std::ifstream input(i_json_file);
			std::cout << input.rdbuf() << '\n';
		}

		std::cout << "Computed timeslots: " << timeslots << '\n';
		std::cout << "Output directory: " << output_dir << std::endl;
	}

	const std::string& get_json_filename()
	{
		return i_json_file;
	}

	const Power_Values& tx_powerlist() const
	{
		return bs_tx_power_dBm;
	}

	const Scan_Values& tx_alphalist() const
	{
		return bs_scan_alpha_deg;
	}

	void init()
	{
		std::cout << "init args" << std::endl;

		if (!i_json_file.empty())
		{
			load_from_json();
		}

		// These checks must also run in Release builds; assert() is not input validation.
		auto require = [](bool valid, const std::string& message)
		{
			if (!valid)
			{
				throw std::invalid_argument("Invalid configuration: " + message);
			}
		};

		auto positive = [](double value)
		{
			return std::isfinite(value) && value > 0;
		};

		require(positive(frequency), "frequency must be positive and finite");
		require(positive(bandwidth), "bandwidth must be positive and finite");
		require(positive(symrate), "symbolrate must be positive and finite");
		require(positive(blockspersym), "blockpersymbol must be positive and finite");
		require(positive(bs_antenna_height), "antenna_height must be positive and finite");
		require(std::isfinite(gain_gtrx), "ms_grx_dB must be finite");
		require(std::isfinite(system_noise), "system_noise_dB must be finite");
		require(std::isfinite(ms_sinr_limit_dB), "sinr_limit_dB must be finite");
		require(base_station_count > 0, "base_stations must be positive");
		require(mobile_station_count >= base_station_count, "number of mobile stations must be at least that of base stations");

		// One client per TX per slot. Division before addition avoids unsigned overflow.
		timeslots = mobile_station_count / base_station_count;
		if (mobile_station_count % base_station_count != 0)
		{
			++timeslots;
		}

		require(field_size.size() == 2 && field_size[0] > 0 && field_size[1] > 0,
			"area must contain positive X and Y dimensions");
		require(tx_loc.size() == base_station_count, "base_station_location count must match base_stations");
		require(rx_loc.size() == mobile_station_count, "mobile_station_location count must match mobile_stations");

		// Each range gives two absolute corners. Neither rectangle is repositioned.
		require(tx_location_range.size() == 2, "base_station_location_range_meters needs two corners");
		require(rx_location_range.size() == 2, "mobile_station_location_range_meters needs two corners");
		const auto& tx_min = tx_location_range.data[0];
		const auto& tx_max = tx_location_range.data[1];
		const auto& rx_min = rx_location_range.data[0];
		const auto& rx_max = rx_location_range.data[1];

		require(tx_min.x < tx_max.x && tx_min.y < tx_max.y,
			"base_station_location_range_meters corners must increase on both axes");
		require(rx_min.x < rx_max.x && rx_min.y < rx_max.y,
			"mobile_station_location_range_meters corners must increase on both axes");
		require(tx_max.x <= field_size[0] && tx_max.y <= field_size[1],
			"base_station_location_range_meters must fit inside area");
		require(rx_max.x <= field_size[0] && rx_max.y <= field_size[1],
			"mobile_station_location_range_meters must fit inside area");

		auto is_inside = [](const Placements& point, const Location_Setup& range)
		{
			const auto& minimum = range.data[0];
			const auto& maximum = range.data[1];
			return minimum.x <= point.x && point.x <= maximum.x &&
				minimum.y <= point.y && point.y <= maximum.y;
		};

		for (size_t index = 0; index < rx_loc.size(); ++index)
		{
			require(is_inside(rx_loc.data[index], rx_location_range),
				"mobile_station_location client " + rx_loc.data[index].str() +
				" must be inside the RX rectangle");
		}

		for (size_t index = 0; index < tx_loc.size(); ++index)
		{
			const auto& point = tx_loc.data[index];
			require(is_inside(point, tx_location_range),
				"base_station_location transmitter " + std::to_string(index + 1) +
				" must be inside the TX rectangle");
			// Negating the inclusive check also excludes every RX edge.
			require(!is_inside(point, rx_location_range),
				"base_station_location transmitter " + std::to_string(index + 1) +
				" must be outside the RX rectangle, including its edges");
		}

		// Validate every per-transmitter value, not just the first table column.
		auto validate_bounds = [&](const double_v& bounds, const std::string& name)
		{
			require(bounds.size() == 2, name + " must contain [minimum, maximum]");
			require(std::isfinite(bounds[0]) && std::isfinite(bounds[1]) && bounds[0] <= bounds[1],
				name + " must have finite, ordered bounds");
		};

		validate_bounds(bs_theta_c_range_deg, "base_station_theta_c_range_deg");
		validate_bounds(bs_power_range_dBm, "base_station_power_range_dBm");
		validate_bounds(bs_scan_angle_range_deg, "base_station_scan_angle_range_deg");
		require(bs_theta_c.size() == base_station_count, "base_station_theta_c needs one value per transmitter");
		require(bs_antenna_count.size() == base_station_count, "base_station_antenna_counts needs one value per transmitter");
		require(bs_antenna_count_limit.size() == base_station_count,
			"base_station_antenna_counts_limit needs one value per transmitter");

		// Half-wavelength spacing is a CLI default only. JSON must supply its spacing list.
		if (i_json_file.empty() && bs_antenna_spacing_m.empty())
		{
			bs_antenna_spacing_m.assign(base_station_count, defaults_t::antenna_spacing(frequency));
		}
		require(bs_antenna_spacing_m.size() == base_station_count, "antenna_spacing needs one value per transmitter");
		require(bs_antenna_dims.size() == 2 && positive(bs_antenna_dims[0]) && positive(bs_antenna_dims[1]),
			"antenna_dims must contain two positive, finite dimensions");

		for (unsigned column = 0; column < base_station_count; ++column)
		{
			require(std::isfinite(bs_theta_c[column]) && bs_theta_c_range_deg[0] <= bs_theta_c[column] &&
				bs_theta_c[column] <= bs_theta_c_range_deg[1], "base_station_theta_c is outside its declared bounds");
			require(bs_antenna_count[column] > 0 && bs_antenna_count[column] <= bs_antenna_count_limit[column],
				"base_station_antenna_counts must be positive and within the per-transmitter limit");
			require(positive(bs_antenna_spacing_m[column]), "antenna_spacing values must be positive and finite");
		}

		// Complete omitted CLI schedules in input units. JSON schedules are never filled in.
		//if (i_json_file.empty())
		//{
		//	if (bs_tx_power_dBm.data.empty())
		//	{
		//		bs_tx_power_dBm.data.assign(timeslots, double_v(base_station_count, defaults_t::bs_tx_power_dBm));
		//	}
		//	if (bs_scan_alpha_deg.data.empty())
		//	{
		//		bs_scan_alpha_deg.data.assign(timeslots, double_v(base_station_count, defaults_t::scan_angle));
		//	}
		//	if (ms_id_selections.data.empty())
		//	{
		//		ms_id_selections.data.assign(timeslots, unsigned_v(base_station_count));
		//		unsigned next_client = 1;
		//		for (auto& row : ms_id_selections.data)
		//		{
		//			for (auto& client : row)
		//			{
		//				client = next_client;
		//				next_client = next_client == mobile_station_count ? 1 : next_client + 1;
		//			}
		//		}
		//	}
		//}

		auto validate_table = [&](const auto& table, const std::string& name)
		{
			require(table.size() == timeslots, name + " must have " + std::to_string(timeslots) + " timeslot rows");
			for (size_t row = 0; row < table.size(); ++row)
			{
				require(table[row].size() == base_station_count,
					name + " row " + std::to_string(row + 1) + " must have one column per transmitter");
			}
		};

		validate_table(bs_tx_power_dBm.data, "base_station_power_dBm_lut");
		validate_table(bs_scan_alpha_deg.data, "base_station_scan_alpha_deg_lut");
		validate_table(ms_id_selections.data, "base_to_mobile_station_id_selection_lut");

		for (unsigned row = 0; row < timeslots; ++row)
		{
			for (unsigned column = 0; column < base_station_count; ++column)
			{
				const double power = bs_tx_power_dBm.data[row][column];
				const double angle = bs_scan_alpha_deg.data[row][column];
				require(std::isfinite(power) && bs_power_range_dBm[0] <= power && power <= bs_power_range_dBm[1],
					"base_station_power_dBm_lut has a value outside its declared bounds");
				require(std::isfinite(angle) && bs_scan_angle_range_deg[0] <= angle && angle <= bs_scan_angle_range_deg[1],
					"base_station_scan_alpha_deg_lut has a value outside its declared bounds");

				if (row > 0)
				{
					const double change = power - bs_tx_power_dBm.data[row - 1][column];
					require(std::abs(change) < 1e-9 || std::abs(std::abs(change) - 1.0) < 1e-9,
						"base_station_power_dBm_lut changes must be -1, 0, or +1 dB per timeslot");
				}
			}
		}

		std::unordered_set<unsigned> served_clients;
		for (const auto& row : ms_id_selections.data)
		{
			const std::unordered_set<unsigned> unique_clients(row.begin(), row.end());
			require(unique_clients.size() == row.size(), "ms_selection must contain distinct clients within each timeslot");
			for (unsigned client : row)
			{
				require(client >= 1 && client <= mobile_station_count, "ms_selection client IDs must be in 1..mobile_stations");
				served_clients.insert(client);
			}
		}
		require(served_clients.size() == mobile_station_count, "ms_selection must cover every client");

		// All validation succeeded. Convert input units and 1-based bindings exactly here.
		const bool valid_power = bs_tx_power_dBm.is_valid_and_convert(bs_power_range_dBm);
		const bool valid_scan = bs_scan_alpha_deg.is_valid_and_convert(bs_scan_angle_range_deg);
		const bool valid_bindings = ms_id_selections.is_valid_and_convert(base_station_count, mobile_station_count);
		require(valid_power && valid_scan && valid_bindings, "schedule conversion failed");

		std::cout << "args built" << std::endl;
	}

	void load_from_json()
	{
		JsonConfig json(i_json_file);

		// Scalars: every physical parameter must come from the JSON file.
		json.read({"frequency"}, frequency);
		json.read({"bandwidth"}, bandwidth);
		json.read({"symbolrate"}, symrate);
		json.read({"blockpersymbol", "blockspersym"}, blockspersym);
		json.read({"antenna_height", "height"}, bs_antenna_height);
		json.read({"ms_grx_dB", "ms_grx"}, gain_gtrx);
		json.read({"system_noise_dB", "system_noise"}, system_noise);
		json.read({"base_stations"}, base_station_count);
		json.read({"mobile_stations"}, mobile_station_count);
		json.read({"sinr_limit_dB", "slimit"}, ms_sinr_limit_dB);

		// Geometry and per-transmitter limits. Old scenario spellings remain supported.
		json.read({"area"}, field_size);
		json.read({"base_station_location"}, tx_loc.data);
		json.read({"mobile_station_location"}, rx_loc.data);
		json.read({"base_station_location_range_meters", "base_station_field_range_meters"}, tx_location_range.data);
		json.read({"mobile_station_location_range_meters", "mobile_station_field_range_meters"}, rx_location_range.data);
		json.read({"base_station_theta_c"}, bs_theta_c);
		json.read({"base_station_theta_c_range_deg", "base_station_theta_c_lim_deg"}, bs_theta_c_range_deg);
		json.read({"base_station_antenna_counts"}, bs_antenna_count);
		json.read({"base_station_antenna_counts_limit", "antenna_count_limit"}, bs_antenna_count_limit);

		// RF bounds, antenna dimensions, and schedules are required, including every table row.
		json.read({"base_station_power_range_dBm"}, bs_power_range_dBm);
		json.read({"base_station_scan_angle_range_deg"}, bs_scan_angle_range_deg);
		json.read({"antenna_spacing"}, bs_antenna_spacing_m);
		json.read({"antenna_dims"}, bs_antenna_dims);
		json.read({"base_station_power_dBm_lut", "antenna_txpower"}, bs_tx_power_dBm.data);
		json.read({"base_station_scan_alpha_deg_lut", "base_station_scan_angle", "scan_angle"}, bs_scan_alpha_deg.data);
		json.read({"base_to_mobile_station_id_selection_lut", "ms_selection"}, ms_id_selections.data);

		// Run controls may be omitted. JSON never sets the output directory.
		timeslot = defaults_t::timeslot;
		json.read({"gui"}, gui, false);
		json.read({"plot"}, plot, false);

		// timeslots is derived in init(), including for files containing a legacy value.
	}
};
