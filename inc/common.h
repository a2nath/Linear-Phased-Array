#pragma once
#include <string>
#include <map>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <exception>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstring>
#include <queue>
#include <iostream>
#include <thread>
#include <future>
#include <condition_variable>


#include "spdlog/spdlog.h"


#define GRAPHICS

#define C_SPEED 3e8L /* speed of light */
#define M_PIl   3.141592653589793238462643383279502884L /* pi */


#ifdef _WIN32
#define localtime(current_time, local_time) localtime_s(local_time, current_time);
#else
#define localtime(current_time, local_time) localtime_r(&current_time, &local_time);
#include <unistd.h>
#endif

extern std::string sim_error;

template<typename T>
std::string str(const T& input, int precision = std::numeric_limits<T>::digits10 + 1)
{
	// For floating-point types, apply precision
	if constexpr (std::is_floating_point<T>::value)
	{
		std::ostringstream ss;
		ss << std::setprecision(precision) << input;
		return ss.str();
	}

	// For non-floating-point types, no special formatting
	return std::to_string(input);
}

template<class U, class V = U> using Pair = std::pair<U, V>;
using unsigned_v = std::vector<unsigned>;
using double_v = std::vector<double>;
using float_v = std::vector<float>;

inline std::string timestamp() {
	// Get the current time as a time_point
	auto now = std::chrono::system_clock::now();

	// Convert to time_t to be able to format the time
	std::time_t currentTime = std::chrono::system_clock::to_time_t(now);

	// Create a tm structure for the local time
	std::tm localTime;

	localtime(&currentTime, &localTime);

	// Use a stringstream to format the timestamp
	std::ostringstream tstamp;
	tstamp << std::put_time(&localTime, "%Y%m%d_%H%M%S");

	return tstamp.str();
}
struct Bench
{
	std::chrono::steady_clock::time_point begin;
	std::chrono::steady_clock::time_point end;
	std::chrono::steady_clock::time_point mark_time;
	long long delta;

	inline auto get()
	{
		return delta;
	}

	inline void bench_start()
	{
		begin = std::chrono::steady_clock::now();
		mark_time = begin;
	}

	inline void mark()
	{
		delta = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - mark_time).count();
		mark_time = std::chrono::steady_clock::now();
	}

	inline void bench_stop()
	{
		end = std::chrono::steady_clock::now();
		delta = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	}
};



template<class T>
inline T rand(T min, T max)
{
	return min + rand() % (max - min + 1);
}

inline std::string getcwd()
{
	return std::filesystem::current_path().string();
}

template<class T, class Compare = std::less<T>>
bool is_unique(std::vector<T>& vec, Compare comp = Compare())
{
	std::sort(vec.begin(), vec.end(), comp);
	return std::adjacent_find(vec.begin(), vec.end()) == vec.end();
}

namespace common
{
	inline bool mkdir(const std::string& path)
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
	template<class T>
	inline void write(const T& contents)
	{
		if (stream.is_open())
		{
			stream << contents;
		}
	}

	template<class T>
	void writeline(const T& contents)
	{
		if (stream.is_open())
		{
			write(contents);
			stream << std::endl;
		}
	}

	void close()
	{
		if (stream.is_open())
		{
			stream.close();
		}
	}

	void copy_settings(path& inputfile, path& output_dirname)
	{
		/* place the input configfile into the output directory to quickly check the inputs */
		std::filesystem::copy_file(inputfile, path(output_dirname / inputfile.filename()).string());

	}

	Logger(path inputfile, path output_dirname, path output_filename = "")
	{
		if (std::filesystem::is_regular_file(inputfile))
		{
			copy_settings(inputfile, output_dirname);
			if (output_filename.string().empty())
			{
				output_filename = inputfile.stem().string() + "." + timestamp() + ".out";
			}
		}
		else
		{
			std::cout << "No reference input file found that could be associated to the sim outputs" << std::endl;
			if (output_filename.string().empty())
			{
				output_filename = "output." + timestamp() + ".out";
			}
		}

		/* open the output file stream during sim */
		stream.open(path(output_dirname / output_filename).string());
	}
};

namespace cached
{
	template<class U, class V> using Cache = std::unordered_map<U, V>;
	template<class U> using HashSet        = std::unordered_set<U>;
	extern Cache<float, float> cache_sin, cache_dBm2w, cache_w2dBm, cache_db2lin, cache_lin2dB, cache_pow, cache_deg2rad;

	inline double deg2rad(const double& deg)
	{
		auto ite = cache_deg2rad.find(deg);
		if (ite == cache_deg2rad.end())
		{
			double answer = deg * M_PIl / 180.0;
			cache_deg2rad.emplace(deg, answer);
			return answer;
		}

		return ite->second;
	}

	inline double pow_2(const double& base)
	{
		auto ite = cache_pow.find(base);
		if (ite == cache_pow.end())
		{
			double answer = pow(base, 2);
			cache_pow.emplace(base, answer);
			return answer;
		}

		return ite->second;
	}

	inline double lin2dB(const double& lin)
	{
		auto ite = cache_lin2dB.find(lin);
		if (ite == cache_lin2dB.end())
		{
			double answer = 10 * log10(lin);
			cache_lin2dB.emplace(lin, answer);
			return answer;
		}

		return ite->second;
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

	inline double watt2dBm(const float& watt)
	{
		auto ite = cache_w2dBm.find(watt);
		if (ite == cache_w2dBm.end())
		{
			double answer = lin2dB(watt) + 30;
			cache_w2dBm.emplace(watt, answer);
			return answer;
		}

		return ite->second;
	}

	inline double dBm2watt(double dBm)
	{
		auto ite = cache_dBm2w.find(dBm);
		if (ite == cache_dBm2w.end())
		{
			double answer = log2lin(dBm - 30);
			cache_dBm2w.emplace(dBm, answer);
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



template<class Type>
struct Dimensions
{
	Type x, y;

	Dimensions& operator=(const Dimensions& ref)
	{
		x = ref.x;
		y = ref.y;
		return *this;
	}
	Dimensions& operator=(const Type& ref)
	{
		x = ref;
		y = ref;
		return *this;
	}

	const std::string str() const
	{
		return "{ x:" + std::to_string(x) + ",y:" + std::to_string(y) + " }";
	}

	bool operator==(const Dimensions& ref)
	{
		return x == ref.x && y == ref.y;
	}

	bool operator!=(const Dimensions& ref)
	{
		return !operator==(ref);
	}

	friend bool operator==(const Dimensions& d1, const Dimensions& d2)
	{
		return d1.x == d2.x && d1.y == d2.y;
	}

	friend bool operator!=(const Dimensions& d1, const Dimensions& d2)
	{
		return d1.x != d2.x || d1.y != d2.y;
	}

	bool is_zero() const
	{
		return x == 0 || y == 0;
	}

	Dimensions(const Type& i, const Type& j) : x(i), y(j) {}
	Dimensions() : x(0), y(0) {}
};

template<class Type>
struct Coordinates : public Dimensions<Type>
{
	Coordinates& operator-(const Coordinates& a)
	{
		Coordinates temp = *this;

		if ((long double)temp.x - (long double)a.x < 0
			|| (long double)temp.y - (long double)a.y < 0)
		{
			throw std::runtime_error("Tried to subtract too much from coordinates data structure leading to less than zero");
		}
		else
		{
			temp.x -= a.x;
			temp.y -= a.y;
		}

		return temp;
	}

	Coordinates(Type _x, Type _y) : Dimensions<Type>(_x, _y) {}
	Coordinates() {}
};
using Placements = Coordinates<unsigned>;
using placement_v = std::vector<Placements>;

struct Polar_Coordinates : private Coordinates<double>
{
	double theta;
	double hype;
	Polar_Coordinates& operator=(const Polar_Coordinates& input)
	{
		theta = input.theta;
		hype = input.hype;
		return *this;
	}
	Polar_Coordinates(double i, double j) : Coordinates(i, j), theta(x), hype(y) {}
	Polar_Coordinates() : Coordinates(), theta(x), hype(y) {}
};

struct PolarArray
{
	Polar_Coordinates* data_ptr;
	size_t array_size;

	void set(const size_t& size)
	{
		if (array_size)
		{
			data_ptr = (Polar_Coordinates*)realloc(data_ptr, size * sizeof(Polar_Coordinates));
		}
		else
		{
			data_ptr = (Polar_Coordinates*)malloc(size * sizeof(Polar_Coordinates));
		}

		if (data_ptr == NULL)
		{
			spdlog::error("Memory allocation failed");
			throw std::runtime_error("Ask for less memory");
		}

		array_size = size;
	}

	PolarArray() : array_size(0), data_ptr(nullptr) {}
};

inline Polar_Coordinates cart2pol(const double& x, const double& y)
{
	double theta = atan2(y, x);
	double hype = hypot(x, y);
	return Polar_Coordinates(theta, hype);
}

template<class T>
inline Polar_Coordinates cart2pol(const Coordinates<T>& c)
{
	return cart2pol(c.x, c.y);
}

inline Coordinates<double> pol2cart(const double& theta, const double& hype)
{
	auto x = cos(theta) * hype;
	auto y = sin(theta) * hype;
	return Coordinates<double>(x, y);
}

inline Coordinates<double> pol2cart(const Polar_Coordinates& c)
{
	return pol2cart(c.theta, c.hype);
}


/* antenna states */
using antennadim = Dimensions<float>;
struct Settings
{
	float     power; // -30 to +30 dBm, default value is 0
	float     alpha; // -90 to +90, default value is 91 or 1.58825
	unsigned   panel_count;
	float     lambda;
	float     spacing;
	float     theta_c;
	antennadim antenna_dims;

	friend bool operator==(const Settings& settings1, const Settings& settings2)
	{
		return settings1.power == settings2.power &&
			settings1.alpha == settings2.alpha &&
			settings1.panel_count == settings2.panel_count &&
			settings1.lambda == settings2.lambda &&
			settings1.spacing == settings2.spacing &&
			settings1.theta_c == settings2.theta_c &&
			settings1.antenna_dims == settings2.antenna_dims;
	}

	friend bool operator!=(const Settings& settings1, const Settings& settings2)
	{
		return !operator==(settings1, settings2);
	}


	const std::string str() const
	{
		return "{ power:" + std::to_string(power)
			+ ", alpha:" + std::to_string(alpha)
			+ ", panel_count:" + std::to_string(panel_count)
			+ ", lambda:" + std::to_string(lambda)
			+ ", spacing:" + std::to_string(spacing)
			+ ", spacing:" + std::to_string(spacing)
			+ ", theta_c:" + std::to_string(theta_c)
			+ ", antenna_dims:" + antenna_dims.str() + " }";
	}


	Settings(const double& ipower,
		const double& ialpha,
		const unsigned& ipanel_count,
		const double& ilambda,
		const double& ispacing,
		const double& itheta_c,
		const antennadim& iantenna_dims)
		:
		power(ipower),
		alpha(ialpha),
		panel_count(ipanel_count),
		lambda(ilambda),
		spacing(ispacing),
		theta_c(itheta_c),
		antenna_dims(iantenna_dims)
	{
	}

	Settings() :
		power(std::numeric_limits<double>::max()),
		alpha(std::numeric_limits<double>::max()),
		panel_count(std::numeric_limits<unsigned>::max()),
		lambda(std::numeric_limits<double>::max()),
		spacing(std::numeric_limits<double>::max()),
		theta_c(std::numeric_limits<double>::max()),
		antenna_dims(std::numeric_limits<float>::max(), std::numeric_limits<float>::max())
	{
	}
};


namespace graphics
{
	extern std::mutex compute_sim_mutex;          // Mutex to protect the shared queue
	extern std::mutex render_mutex;
	extern std::mutex graphics_data_mutex;
	extern std::condition_variable consig;      // Condition variable to signal the consumer thread
	extern Dimensions<unsigned> render_space;

	struct State
	{
		int tx_idx;
		Settings settings;
		Placements location;

		State& operator=(const State& b)
		{
			settings = b.settings;
			location = b.location;
			return *this;
		}

		friend bool operator==(const State& a, const State& b)
		{
			return true;// a.settings == b.settings && a.location == b.location;
		}

		friend bool operator!=(const State& a, const State& b)
		{
			return a.settings != b.settings || a.location != b.location;
		}

		State(const int& id = -1) : tx_idx(id)
		{
			settings.alpha = -1;
			settings.lambda = -1; // TODO
			settings.panel_count = 0; // TODO
			settings.power = -1;
			settings.spacing = -1; // TODO
			settings.theta_c = -1;
		}
	};

	struct DataSync
	{
		/* signal the tx number */
		int compute_tx_id, render_tx_id, def_render_tx_id;
		bool resize_event;
		std::queue<State*> mainq;
		std::vector<std::queue<State*>> pending;

		inline void emplace_state(State& state)
		{
			std::scoped_lock<std::mutex> lock(graphics::graphics_data_mutex);  // Lock the mutex
			mainq.emplace(&state);

			compute_tx_id = 1;
		}

		inline void event_resize(const unsigned& width, const unsigned& height, const int& def_render_id)
		{
			def_render_tx_id = def_render_id;
			if (width != render_space.x || height != render_space.y)
			{
				resize_event = true;
				render_space = { width, height };
			}
		}

		inline void event_render(const int& idx)
		{
			std::scoped_lock<std::mutex> lock(graphics::render_mutex);  // Lock the mutex
			render_tx_id = idx;
		}

		bool got_updates(const int& def_render_id)
		{
			def_render_tx_id = def_render_id;
			return compute_tx_id > 0 || render_tx_id >= 0 || resize_event;
		}

		DataSync(int num_tx) : compute_tx_id(0), render_tx_id(-1), def_render_tx_id(-1), \
			resize_event(false), pending(num_tx)
		{
		}
	};
}
using state_v = std::vector<graphics::State>;
