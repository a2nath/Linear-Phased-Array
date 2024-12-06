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

template<class U, class V> using Pair = std::pair<U, V>;
using unsigned_v = std::vector<unsigned>;
using double_v = std::vector<double>;

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
	extern Cache<double, double> cache_sin, cache_dbm2w, cache_db2lin, cache_pow;

	inline double deg2rad(const double& deg)
	{
		return deg * M_PIl / 180.0;
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



template<class Width>
struct Dimensions
{
	Width x, y;

	Dimensions& operator=(const Dimensions& ref)
	{
		x = ref.x;
		y = ref.y;
		return *this;
	}
	Dimensions& operator=(const Width& ref)
	{
		x = ref;
		y = ref;
		return *this;
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
		return !operator==(d1, d2);
	}


	Dimensions(const Width& i, const Width& j) : x(i), y(j) {}
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
	double& theta;
	double& hype;
	Polar_Coordinates& operator=(const Polar_Coordinates& input)
	{
		theta = input.theta;
		hype = input.hype;
		return *this;
	}
	Polar_Coordinates(double i, double j) : Coordinates(i, j), theta(x), hype(y) {}
	Polar_Coordinates() : Coordinates(), theta(x), hype(y) {}
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

#include <thread>
#include <future>
#include <condition_variable>

namespace graphics
{
	extern std::mutex queue_mutex;          // Mutex to protect the shared queue
	extern std::mutex finished_mutex;
	extern std::condition_variable consig;      // Condition variable to signal the consumer thread

	struct State
	{
		int tx_idx;
		double ant_power, ant_dir, ant_angle;
		long row, col;
		Placements location;

		State(const int& id, const long& irow, const long& icol, const double& power, const double& dir, const double& scan, const long& x, const long& y) :
			tx_idx(id), row(irow), col(icol), ant_power(power), ant_dir(dir), ant_angle(scan), location({ (unsigned)x, (unsigned)y }) {}
		State() : tx_idx(-1), row(-1), col(-1), ant_power(-1), ant_dir(-1), ant_angle(-1), location() {}
	};

	struct DataSync
	{
		/* signal the tx number */
		int finished;
		long size;
		std::queue<State> mainq;
		std::vector<std::queue<State*>> pending;


		auto& front()
		{
			std::lock_guard<std::mutex> lock(queue_mutex);  // Lock the mutex

			auto state = mainq.front();
			State new_state;

			/* is this in the pending references */
			while (pending[state.tx_idx].front() != &state)
			{
				mainq.pop();

				if (mainq.size())
				{
					state = mainq.front();
				}
				else
				{
					state = new_state;
				}
			}

			/* the newest valid request for that tx */
			return state;
		}

		void pop()
		{
			std::lock_guard<std::mutex> lock(queue_mutex);  // Lock the mutex

			if (mainq.size())
			{
				auto& state = mainq.front();
				pending[state.tx_idx].pop();
				mainq.pop();

				std::unique_lock<std::mutex> lock(graphics::finished_mutex);
				finished = state.tx_idx;

				size--;

				consig.notify_one();
			}
		}

		void emplace(const int& idx, const long& row, const long& col, const double& power, const double& dir, const double& scan, const long& x, const long& y)
		{
			std::lock_guard<std::mutex> lock(queue_mutex);  // Lock the mutex

			/* queue the latest request from tx */
			mainq.emplace(idx, row, col, power, dir, scan, x, y);

			pending[idx].emplace(&mainq.back());


			/* just process the newest request from other tx */
			for (int i = 0; i < pending.size(); ++i)
			{
				if (i != idx && pending[i].size() > 1)
				{
					std::queue<State*> empty_queue({ pending[i].back() });
					std::swap(pending[i], empty_queue);
				}
			}

			size++;

			consig.notify_one();  // Notify the consumer thread that there's new data
		}

		DataSync(int num_tx) : size(0), finished(-1), pending(num_tx)
		{

		}
	};
}
