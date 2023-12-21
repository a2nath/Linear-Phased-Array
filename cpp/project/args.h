#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <any>
#include <iostream>
#include <fstream>
#include "common.h"
#include "network.h"
#include "json/json.h"
#include "../inc/argparse/argparse.hpp"

using namespace network_package;

class ProgramConfig
{
protected:
    /* simulation system parameters */

    template<typename T = double> const auto& param(const std::string& key) const
    {
        return std::any_cast<const T&>(data.at(key));
    }

    /* antenna power info info */
    int get_power_lut(const std::vector<std::vector<double>>& input, std::vector<std::vector<double>>& powertable)
    {
        powertable.resize(input.size(), std::vector<double>(input[0].size()));

        for (unsigned i = 0; i < input.size(); ++i)
        {
            for (unsigned j = 0; j < input[i].size(); ++j)
            {
                powertable[i][j] = dBm2watt(input[i][j]);
            }
        }

        return 0;
    }

    /* antenna array direction info */
    int get_scan_angle_lut(const std::vector<std::vector<double>>& input, std::vector<std::vector<double>>& alphatable)
    {


        alphatable.resize(input.size(), std::vector<double>(input[0].size()));

        for (unsigned i = 0; i < input.size(); ++i)
        {
            for (unsigned j = 0; j < input[i].size(); ++j)
            {
                alphatable[i][j] = deg2rad(input[i][j]);
            }
        }

        return 0;
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
            output_power.emplace_back(dBm2watt(input[i]));
        }
        return 0;
    }

    /* get list of scan angle for each base station in a timeslot */
    int get_scan_angle_list(const std::vector<double>& input, std::vector<double>& output_scan_angle)
    {
        for (auto& data : input)
        {
            output_scan_angle.emplace_back(deg2rad(data));
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
            antenna_dir.emplace_back(deg2rad(data));
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
public:
    virtual void parse_data(Input*& args)
    {
        auto& power_range = param<std::vector<double>>("base_station_power_range_dBm");
        auto& scana_range = param<std::vector<double>>("base_station_scan_angle_range_deg");

        args = new Input(
            param("frequency"),
            param("bandwidth"),
            param("SymbRate"),
            param("BlocksPerSymb"),
            param("height"),
            power_range[0],
            power_range[1],
            param("ms_Grx"),
            scana_range[0],
            scana_range[1],
            param("system_noise"),
            param<unsigned>("mobile_stations"),
            param<unsigned>("base_stations")
        );

        get_base_station_cords(param<std::vector<Placements>>("base_station_location"), base_station_pos);
        get_mobile_station_cords(param<std::vector<Placements>>("mobile_station_location"), mobile_station_pos);

        //if (data.find("mobile_station_location") != data.end())
        //    get_power_list(param<std::vector<double>>("base_station_power_dBm"), args->ant.array_power_wtts);

        if (data.find("base_station_power_dBm_lut") != data.end())
            get_power_lut(param<std::vector<std::vector<double>>>("base_station_power_dBm_lut"), args->ant.array_power_wtts_lut);

        //if (data.find("base_station_scan_alpha_deg") != data.end())
        //    get_scan_angle_list(param<std::vector<double>>("base_station_scan_alpha_deg"), args->ant.array_scan_angle);

        if (data.find("base_station_scan_alpha_deg_lut") != data.end())
            get_scan_angle_lut(param<std::vector<std::vector<double>>>("base_station_scan_alpha_deg_lut"), args->ant.array_scan_angle_lut);

        //if (data.find("base_to_mobile_station_id_selection") != data.end())
        //    get_station_list(param<std::vector<unsigned>>("base_to_mobile_station_id_selection"), station_ids);

        if (data.find("base_to_mobile_station_id_selection_lut") != data.end())
            get_station_lut(param<std::vector<std::vector<unsigned>>>("base_to_mobile_station_id_selection_lut"), station_ids_lut);

        //get_power_list(param<std::vector<double>>("base_station_power_dBm"), args->ant.array_power_wtts);
        get_power_lut(param<std::vector<std::vector<double>>>("base_station_power_dBm_lut"), args->ant.array_power_wtts_lut);
        //get_scan_angle_list(param<std::vector<double>>("base_station_scan_alpha_deg"), args->ant.array_scan_angle);
        get_scan_angle_lut(param<std::vector<std::vector<double>>>("base_station_scan_alpha_deg_lut"), args->ant.array_scan_angle_lut);
        //get_station_list(param<std::vector<unsigned>>("base_to_mobile_station_id_selection"), station_ids);
        get_station_lut(param<std::vector<std::vector<unsigned>>>("base_to_mobile_station_id_selection_lut"), station_ids_lut);
        get_allocation_per_bs(param<std::vector<unsigned>>("base_station_antenna_counts"), args->ant.antcount_per_base);
        get_base_station_orientation(param<std::vector<double>>("base_station_theta_c"), args->ant.antenna_orientation);
        get_base_station_antspace(param<std::vector<double>>("antenna_spacing"), args->ant.antenna_spacing);
        get_base_station_anntena_dims(param<std::vector<double>>("antenna_dims"), args->ant.antenna_dim_mtrs);
    }
};

class Default;
std::filesystem::path operator+(const std::filesystem::path dir, const std::filesystem::path& name)
{
    return dir.string() + '/' + name.string();
}

/* argparse for C++17 */
struct MyArgs : public argparse::Args, public ProgramConfig {

    std::vector<std::string> json_files;

    //std::optional<bool>& defaults                             = kwarg("d,default", "Use defaults from the program without any input file");
    std::optional<std::string>& input_dir                     = kwarg("i,input_dir", "Input directory of the test files");
    std::optional<std::string>& output_dir                    = kwarg("o,output_dir", "Output directory to put results in");
    std::optional<std::string>& json_file                     = kwarg("f,test,test_file", "Input parameter file to run in .json format");
    std::optional<double>& frequency                          = kwarg("frequency", "Frequency of the signal system wide");
    std::optional<double>& bandwidth                          = kwarg("bandwidth", "Bandwidth of the singal");
    std::optional<double>& symrate                            = kwarg("symrate", "Symbol rate of the data stream");
    std::optional<double>& blockspersym                       = kwarg("blockps", "Blocks of data per symbol");
    std::optional<double>& antenna_height                     = kwarg("height", "Height of the antenna");
    std::optional<double>& gain_gtrx                          = kwarg("ms_grx", "Antenna gain of the antenna at the mobile or client station in dB");
    std::optional<double>& system_noise                       = kwarg("system_noise", "System noise in the transmitter and receiver");
    std::optional<unsigned>& mobile_stations_count            = kwarg("mobile_stations", "Number of mobile stations in the simulation");
    std::optional<unsigned>& base_stations_count              = kwarg("base_stations", "Number of base stations in the simulation");
    std::optional<std::vector<double>>& bs_theta_c            = kwarg("theta_c", "Directions antennas are placed at");
    std::optional<std::vector<unsigned>>& base_stations_loc   = kwarg("base_stations_loc", "Location of base stattions, get a list");
    std::optional<std::vector<unsigned>>& mobile_stations_loc = kwarg("mobile_stations_loc", "Location of mobile stations, get a list");
    std::optional<std::vector<unsigned>>& bs_antenna_counts   = kwarg("antenna_counts", "Location of mobile stations, get a list");
    std::optional<std::vector<double>>& power_range           = kwarg("power_range", "Range of power that base stations will use in dBm");
    std::optional<std::vector<double>>& scan_angle_range      = kwarg("scan_angle_range", "Scan angle of the antenna linear array at the base station in degrees").set_default("-90,90");;
    std::optional<std::vector<double>>& antenna_spacing       = kwarg("antenna_spacing", "Antenna spacing between panels");
    std::optional<std::vector<double>>& antenna_dims          = kwarg("antenna_dims", "Antenna dimensions in meters");
    std::optional<std::vector<double>>& bs_tx_power_dBm       = kwarg("antenna_txpower", "Base station transmit TX power in dBm");
    std::optional<std::vector<std::vector<double>>>& bs_tx_power_dBm_lut    = kwarg("antenna_txpower_lut", "Base station transmit TX power in dBm lookup table");
    std::optional<std::vector<double>>& bs_scan_alpha_deg                   = kwarg("scan_angle", "Base station scan angle");
    std::optional<std::vector<std::vector<double>>>& bs_scan_alpha_deg_lut  = kwarg("scan_angle_lut", "Base station scan angle lookup table");
    std::optional<std::vector<unsigned>>& ms_selection                      = kwarg("ms_selection", "Mobile station id selection");
    std::optional<std::vector<std::vector<unsigned>>>& ms_selection_lut     = kwarg("ms_selection_lut", "Mobile station id selection lookup table");


    /* get a list of config files in json format */
    int get_tests(const std::string& directory)
    {
        std::string path(directory);
        std::string ext(".json");
        for (auto& p : std::filesystem::recursive_directory_iterator(path))
        {
            if (p.path().extension().compare(ext) == 0) // if that extension
                json_files.emplace_back(p.path().string());
        }
        return json_files.size();
    }

    bool isFile(std::string& path)
    {
        bool answer = false;

        const std::string& dir = input_dir != std::nullopt && isDir(input_dir.value()) ? input_dir.value() : getcwd();
        std::filesystem::path filename(path);

        if (filename.parent_path().string().size() == 0) // just a filename so append input dir
        {
            answer = std::filesystem::is_regular_file(dir + filename.filename());
            path = dir + filename.filename().string();   // filename has input dir as well
        }
        else
        {
            answer = std::filesystem::is_regular_file(path);
        }

        json_files.emplace_back(path);

        return answer;
    }

    bool isDir(std::filesystem::path path)
    {

        bool answer = std::filesystem::is_directory(path);

        if (!answer)
            std::cout << "WARNING: " << "input directory is not valid" << std::endl;

        return answer;
    }

    /* check if default values to be used from the hardcoded file */
    bool isDefault()
    {
        bool answer = !(
            /* check if the path is correct and whether there are valid tests to add, if file is not separately specified */
               (input_dir != std::nullopt && json_file == std::nullopt && isDir(input_dir.value()) && get_tests(input_dir.value()) > 0)

            /* check if the ouput dir is valid */
            || (output_dir != std::nullopt && isDir(output_dir.value()))

            /* check if the file is valid and whether it can be added to the list with or without the input dir */
            || (json_file != std::nullopt && isFile(json_file.value()))

            /* check the rest of the parameters if they are defined */
            || frequency != std::nullopt
            || bandwidth != std::nullopt
            || symrate != std::nullopt
            || blockspersym != std::nullopt
            || antenna_height != std::nullopt
            || gain_gtrx != std::nullopt
            || system_noise != std::nullopt
            || mobile_stations_count != std::nullopt
            || base_stations_count != std::nullopt
            || bs_theta_c != std::nullopt
            || base_stations_loc != std::nullopt
            || mobile_stations_loc != std::nullopt
            || bs_antenna_counts != std::nullopt
            || power_range != std::nullopt
            || scan_angle_range != std::nullopt
            || antenna_spacing != std::nullopt
            || antenna_dims != std::nullopt
            || bs_tx_power_dBm != std::nullopt
            || bs_tx_power_dBm_lut != std::nullopt
            || bs_scan_alpha_deg != std::nullopt
            || bs_scan_alpha_deg_lut != std::nullopt
            || ms_selection != std::nullopt
            || ms_selection_lut != std::nullopt
            );

        /* add at least one so that test can start */
        if (json_files.empty())
            json_files.emplace_back("");

        return answer;
    }

    void interprete()
    {
        data.emplace("frequency", frequency.value());
        data.emplace("bandwidth", bandwidth.value());
        data.emplace("SymbRate", symrate.value());
        data.emplace("BlocksPerSymb", blockspersym.value());
        data.emplace("height", antenna_height.value());
        data.emplace("ms_Grx", gain_gtrx.value());
        data.emplace("system_noise", system_noise.value());
        data.emplace("mobile_stations", mobile_stations_count.value());
        data.emplace("base_stations", base_stations_count.value());

        if (bs_antenna_counts.value().size() != base_stations_count.value())
            throw std::runtime_error("Antenna allocation info is incorrect");

        data.emplace("base_station_antenna_counts", bs_antenna_counts.value());

        if (power_range.value().size() != 2 || scan_angle_range.value().size() != 2)
            throw std::runtime_error("Power range or scan angle range is not correct");

        data.emplace("base_station_power_range_dBm", power_range.value());
        data.emplace("base_station_scan_angle_range_deg", scan_angle_range.value());


        if (bs_theta_c.value().size() != base_stations_count.value())
            throw std::runtime_error("Theta_C info is incorrect");

        data.emplace("base_station_theta_c", bs_theta_c.value());

        if (antenna_spacing.value().size() != base_stations_count.value())
            throw std::runtime_error("Antenna spacing is not correct");

        data.emplace("antenna_spacing", antenna_spacing.value());

        if (antenna_dims.value().size() != 2)
            throw std::runtime_error("Antenna dimension is not correct");

        data.emplace("antenna_dims", antenna_dims.value());

        if (base_stations_loc.value().size() % 2 > 0 || mobile_stations_loc.value().size() % 2 > 0)
            throw std::runtime_error("Station location data is not correct");

        std::vector<Placements> bs_locs, ms_locs;
        for (unsigned i = 0; i < base_stations_loc.value().size(); i += 2)
        {
            bs_locs.emplace_back(base_stations_loc.value()[i], base_stations_loc.value()[i + 1]);
        }
        data.emplace("base_station_location", bs_locs);

        for (unsigned i = 0; i < mobile_stations_loc.value().size(); i += 2)
        {
            ms_locs.emplace_back(mobile_stations_loc.value()[i], mobile_stations_loc.value()[i + 1]);
        }
        data.emplace("mobile_station_location", ms_locs);

        if (bs_tx_power_dBm != std::nullopt)
            data.emplace("base_station_power_dBm", bs_tx_power_dBm.value());

        if (bs_tx_power_dBm_lut != std::nullopt)
            data.emplace("base_station_power_dBm_lut", bs_tx_power_dBm_lut.value());

        if (bs_scan_alpha_deg != std::nullopt)
            data.emplace("base_station_scan_alpha_deg", bs_scan_alpha_deg.value());

        if (bs_scan_alpha_deg_lut != std::nullopt)
            data.emplace("base_station_scan_alpha_deg_lut", bs_scan_alpha_deg_lut.value());

        if (ms_selection != std::nullopt)
            data.emplace("base_to_mobile_station_id_selection", ms_selection.value());

        if (ms_selection_lut != std::nullopt)
            data.emplace("base_to_mobile_station_id_selection_lut", ms_selection_lut.value());
    }

    void parse_data(Input*& args)
    {
        auto& power_range = param<std::vector<double>>("base_station_power_range_dBm");
        auto& scana_range = param<std::vector<double>>("base_station_scan_angle_range_deg");

        args = new Input(
            param("frequency"),
            param("bandwidth"),
            param("SymbRate"),
            param("BlocksPerSymb"),
            param("height"),
            power_range[0],
            power_range[1],
            param("ms_Grx"),
            scana_range[0],
            scana_range[1],
            param("system_noise"),
            param<unsigned>("mobile_stations"),
            param<unsigned>("base_stations")
        );

        get_base_station_cords(param<std::vector<Placements>>("base_station_location"), base_station_pos);
        get_mobile_station_cords(param<std::vector<Placements>>("mobile_station_location"), mobile_station_pos);

        if (data.find("base_station_power_dBm") != data.end())
            get_power_list(param<std::vector<double>>("base_station_power_dBm"), args->ant.array_power_wtts);

        if (data.find("base_station_power_dBm_lut") != data.end())
            get_power_lut(param<std::vector<std::vector<double>>>("base_station_power_dBm_lut"), args->ant.array_power_wtts_lut);

        if (data.find("base_station_scan_alpha_deg") != data.end())
            get_scan_angle_list(param<std::vector<double>>("base_station_scan_alpha_deg"), args->ant.array_scan_angle);

        if (data.find("base_station_scan_alpha_deg_lut") != data.end())
            get_scan_angle_lut(param<std::vector<std::vector<double>>>("base_station_scan_alpha_deg_lut"), args->ant.array_scan_angle_lut);

        if (data.find("base_to_mobile_station_id_selection") != data.end())
            get_station_list(param<std::vector<unsigned>>("base_to_mobile_station_id_selection"), station_ids);

        if (data.find("base_to_mobile_station_id_selection_lut") != data.end())
            get_station_lut(param<std::vector<std::vector<unsigned>>>("base_to_mobile_station_id_selection_lut"), station_ids_lut);

        get_allocation_per_bs(param<std::vector<unsigned>>("base_station_antenna_counts"), args->ant.antcount_per_base);
        get_base_station_orientation(param<std::vector<double>>("base_station_theta_c"), args->ant.antenna_orientation);
        get_base_station_antspace(param<std::vector<double>>("antenna_spacing"), args->ant.antenna_spacing);
        get_base_station_anntena_dims(param<std::vector<double>>("antenna_dims"), args->ant.antenna_dim_mtrs);
    }
};

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
    };

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
    };

public:
    JasonHelper(const std::string& file)
    {
        Json::Value root;
        std::ifstream ifs;
        ifs.open(file);

        /* parse the json file and fill in args for the program */
        if (ifs.good())
        {
            Json::CharReaderBuilder builder;
            builder["collectComments"] = false;
            JSONCPP_STRING errs;
            if (!parseFromStream(builder, ifs, &root, &errs)) {
                std::cout << errs << std::endl;
                throw std::runtime_error("Parser problem: " + str(EXIT_FAILURE));
            }

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
            sim_error = "JasonHelper: json file not read";
        }
    }
};

class Defaults : public ProgramConfig
{
    const double   frequency = 1900e6;   // center frequency
    const double   bandwidth = 20e6;     // channel bandwidth in Mz
    const double   SymbRate = 3.84e6;    // symbol rate
    const double   BlocksPerSymb = 768;  // blocks per symbol
    const double   lambda = C / frequency;
    const double   height = 200;         //antenna height in meters
    const double   ms_Grx = -2;          // dBi
    const double   system_noise = 5;     // System noise in dB;
    const unsigned mobile_stations = 15; // number of handdset
    const unsigned base_stations = 3;    // number of base stations

    const std::vector<double> BS_THETAC_DEG = { 75, 90, 105 };    // degrees

    const std::vector<Placements> COW_LOCATION =    // in meters
    {
        {250, 200}, {500, 180}, {750, 200}
    };

    const std::vector<Placements> NODE_LOCATIONS =  // in meters
    {
        {250, 500}, {350, 500}, {450, 500}, {550, 500}, {650, 500}, {750, 500},
        {250, 600}, {350, 600}, {450, 600}, {550, 600}, {650, 600}, {750, 600},
                          {450, 325}, {500, 325}, {550, 325}
    };

    const std::vector<int> BS_ANTENNA_COUNTS = { 5, 3, 5 };        // number of antenna allocation for each BS


    const int MIN_POWER_dBm = -30;       // dBm
    const int MAX_POWER_dBm = +30;       // dBm
    const int MIN_SCANANGLE_deg = -90;   // degrees
    const int MAX_SCANANGLE_deg = +90;   // degrees

    const std::vector<double> ANTENNA_SPACING = { .3, .4, .3 }; // meters
    const std::vector<double> ANTENNA_DIMS = { .2, .4 };     // meters

    const std::vector<std::vector<double>> COW_POWER_LUT = {
        {+30.0, +30.0, +30.0},
        {+30.0, +30.0, +30.0},
        {+30.0, +30.0, +30.0},
        {-30.0, +30.0, -30.0},
        {-30.0, +30.0, -30.0}
    };

    // degrees
    const std::vector<std::vector<double>> SCAN_ALPHA_LUT = {
        {+10.0, +65.0, +35.0},
        {-10.0, -65.0, -35.0},
        {-75.0, -40.0, +75.0},
        {+40.0, +65.0, -40.0},
        {+40.0, -65.0, -40.0}
    };

    /* Setup the environment */
    const std::vector<unsigned> BS_STAID_SELECTION = {

        /* BS [0] [1] [2]  */
               2, 13, 11,
               8, 15, 05,
               1, 14, 06,
               7, 04, 10,
               9, 03, 12
    };

    /* antenna power info info */
    int get_power_lut(std::vector<std::vector<double>>& powertable)
    {
        auto& list = COW_POWER_LUT;

        powertable.resize(list.size(), std::vector<double>(list[0].size()));

        for (unsigned i = 0; i < list.size(); ++i)
        {
            for (unsigned j = 0; j < list[i].size(); ++j)
            {
                powertable[i][j] = dBm2watt(list[i][j]);
            }
        }

        return 0;
    }

    /* antenna array direction info */
    int get_scan_angle_lut(std::vector<std::vector<double>>& alphatable)
    {
        auto& list = SCAN_ALPHA_LUT;

        alphatable.resize(list.size(), std::vector<double>(list[0].size()));

        for (unsigned i = 0; i < list.size(); ++i)
        {
            for (unsigned j = 0; j < list[i].size(); ++j)
            {
                alphatable[i][j] = deg2rad(list[i][j]);
            }
        }

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

    /* phase array antenna count per base station */

    int get_allocation_per_bs(std::vector<unsigned>& output_counts)
    {
        return getData(BS_ANTENNA_COUNTS, output_counts);

    }

    /* get base station location info */

    int get_base_station_cords(std::vector<Placements>& station_position)
    {
        return getData(COW_LOCATION, station_position);

    }

    /* get mobile station location info */
    int get_mobile_station_cords(std::vector<Placements>& station_position)
    {
        return getData(NODE_LOCATIONS, station_position);
    }

    /* this defines the thetaC term in the calculation or where the BS is facing */
    int get_base_station_orientation(std::vector<double>& antenna_dir)
    {
        for (auto& data : BS_THETAC_DEG)

        {
            antenna_dir.emplace_back(deg2rad(data));
        }

        return 0;
    }

    int get_base_station_antspace(std::vector<double>& antenna_spacing)
    {
        return getData(ANTENNA_SPACING, antenna_spacing);
    }

    /* get the size of the antenna in [meters] where antennadim is Dimension<double> */
    int get_base_station_anntena_dims(antennadim& dim)
    {
        auto& data = ANTENNA_DIMS;

        dim.x = data[0];
        dim.y = data[1];

        return 0;
    }
public:
    Defaults(Input*& args)
    {
        args = new Input(
            frequency,
            bandwidth,
            SymbRate,
            BlocksPerSymb,
            height,
            MIN_POWER_dBm,
            MAX_POWER_dBm,
            ms_Grx,
            MIN_SCANANGLE_deg,
            MAX_SCANANGLE_deg,
            system_noise,
            mobile_stations,
            base_stations
        );
        get_mobile_station_cords(mobile_station_pos);
        get_base_station_cords(base_station_pos);
        get_power_lut(args->ant.array_power_wtts_lut);
        get_scan_angle_lut(args->ant.array_scan_angle_lut);
        get_allocation_per_bs(args->ant.antcount_per_base);
        get_base_station_orientation(args->ant.antenna_orientation);
        get_base_station_antspace(args->ant.antenna_spacing);
        get_base_station_anntena_dims(args->ant.antenna_dim_mtrs);
    }
};