#pragma once
#include <string>
#include <map>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <exception>
#include <fstream>

#define C       3e8 //speed of light
#define M_PIl   3.141592653589793238462643383279502884L /* pi */


#ifdef _WIN32
#include <process.h>
/* get process id and attach to process */
#define GetCurrentProcessId getpid
#else
#include <unistd.h>
#endif

extern std::string sim_error;

template<class T> std::string str(const T& input) { return std::to_string(input); }
template<class U, class V> using Pair = std::pair<U, V>;

template<class T>
inline T rand(T min, T max)
{
	return min + rand() % (max - min + 1);
}

inline std::string getcwd()
{
	return std::filesystem::current_path().string();
}

namespace common
{
	bool mkdir(const std::string& path)
	{
		if (!(std::filesystem::exists(path) && std::filesystem::is_directory(path)))
		{
			try
			{
				std::filesystem::create_directories(path);
					return true; // Directory created successfully
			}
			catch (const std::filesystem::filesystem_error& e)
			{
				std::cerr << "Error creating directory " << path << ": " << e.what() << std::endl;
					return false; // Error creating directory
			}
		}

		return true;
	}
}

/* Inverse map with limited size. Starts discarding lower entries with full */
template<class U, class V>
struct revmap_t
{
	std::map<U, V, std::greater<V>> container;
	size_t size;
	void emplace(const U& key, const V& value) override
	{
		//if ()
		container.emplace(key, value);

		if (container.size() == size)
		{
			container.erase(container.rend());
		}
	}
	const V& operator[](const U& key)
	{
		return container[key];
	}
	revmap_t(size_t _size) : size(size) {}
};



class Logger
{
	using path = std::filesystem::path;
	std::ofstream stream;
public:
	inline void write(const std::string& str)
	{
		stream << str;
	}

	void writeline(const std::string& str)
	{
		write(str);
		stream << std::endl;
	}

	void close()
	{
		stream.close();
	}

	Logger(const std::string& filename)
	{
		stream.open(filename);
	}
};

namespace cached
{
	template<class U, class V> using Cache = std::unordered_map<U, V>;
	template<class U> using HashSet        = std::unordered_set<U>;
	Cache<double, double> cache_sin, cache_dbm2w, cache_db2lin;

	inline double deg2rad(const double& deg)
	{
		return deg * M_PIl / 180.0;
	}

	inline double log2lin(double log)
	{
		auto ite = cache_db2lin.find(log);
		if (ite == cache_db2lin.end())
		{
			double answer = pow(10, log / 10);
			cache_db2lin.emplace(log, answer);
			return answer;
		}

		return ite->second;
	}

	inline double dBm2watt(double dBm)
	{
		auto ite = cache_dbm2w.find(dBm);
		if (ite == cache_dbm2w.end())
		{
			double answer = log2lin(dBm - 30);
			cache_dbm2w.emplace(dBm, answer);
			return answer;
		}

		return ite->second;
	}

	inline std::vector<double> deg2rad(const std::vector<double>& log_values)
	{
		std::vector<double> output(log_values.size());
		size_t idx = 0;
		for (auto& item : log_values)
		{
			output[idx++] = deg2rad(item);
		}

		return output;
	}

	inline std::vector<double> dBm2watt(const std::vector<double>& dbm_values)
	{
		std::vector<double> output(dbm_values.size());
		size_t idx = 0;
		for (auto& item : dbm_values)
		{
			output[idx++] = dBm2watt(item);
		}
		return output;
	}

	inline double sin(const double& rads)
	{
		auto ite = cache_sin.find(rads);
		if (ite != cache_sin.end())
		{
			return ite->second;
		}
		else
		{
			double val = std::sin(rads);
			cache_sin.emplace(rads, val);
			return val;
		}
	}

	inline double cos(const double& rads)
	{
		return cached::sin((M_PIl / 2.0) - rads);
	}

}
