#pragma once
#include <random>
#include <cmath>
#include <string>
#include <map>

#define C       3e8 //speed of light
#define M_PIl   3.141592653589793238462643383279502884L /* pi */

extern std::string sim_error;

inline double deg2rad(double deg)
{
    return deg * M_PIl / 180.0;
}
inline double log2lin(double log)
{
    return pow(10, log / 10);
}
inline double dBm2watts(double dBm)
{
    return log2lin(dBm - 30);
}
template<class T>
inline T rand(T min, T max)
{
    return min + rand() % (max - min + 1);
}
struct Dimension
{
    double x, y;
    Dimension(double i, double j) : x(i), y(j) {}
};

template<class T>
struct Coordinates
{
    T x, y;
    Coordinates& operator=(unsigned a)
    {
        x = a;
        y = a;
    }

    Coordinates operator-(Coordinates& a)
    {
        Coordinates temp = *this;

        if ((long double)temp.x - (long double)a.x < 0
            || (long double)temp.y - (long double)a.y < 0)
        {
            sim_error = "Tried to subtract too much from coordinates data structure leading to less than zero";
            temp = 0;
        }
        else
        {
            temp.x -= a.x;
            temp.y -= a.y;
        }

        if (sim_error)
        {
            throw std::runtime_error(sim_error + ", errno " + strerror(errno));
        }

        return temp;
    }

    Coordinates(T _x, T _y) : x(_x), y(_y) {}
    Coordinates() : x(0), y(0) {}
};

struct Polar_Coordinates : private Coordinates<double>
{
    double& theta;
    double& hype;
    Polar_Coordinates& operator=(const Polar_Coordinates& input)
    {
        theta = input.theta;
        hype  = input.hype;
        return *this;
    }
    Polar_Coordinates(double i, double j) : Coordinates(i, j), theta(x), hype(y) {}
    Polar_Coordinates() : Coordinates(), theta(x), hype(y) {}
};

inline Polar_Coordinates cart2pol(double x, double y)
{
    double theta = atan2(y, x);
    double hype = hypot(x, y);
    return Polar_Coordinates(theta, hype);
}

template<class T>
inline Polar_Coordinates cart2pol(Coordinates<T>& c)
{
    return cart2pol(c.x, c.y);
}

inline Coordinates<double> pol2cart(double theta, double hype)
{
    auto x = cos(theta) * hype;
    auto y = sin(theta) * hype;
    return Coordinates<double>(x, y);
}

inline Coordinates<double> pol2cart(Polar_Coordinates& c)
{
    return pol2cart(c.theta, c.hype);
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
        return container[U];
    }
    revmap_t(size_t _size) : size(size) {}
};

/* Random number generator module */
class Random
{

    // First create an instance of an engine.
    std::random_device rnd_device;
    long long unsigned seed;

    // Specify the engine and distribution.
    std::mt19937 mersenne_engine;
public:

    /* Get a random number by passing one distribution obj */
    template<class T>
    inline auto generate(T& dist)
    {
        return dist(mersenne_engine);
    }

    Random() : mersenne_engine(rnd_device)
    {

    }


};