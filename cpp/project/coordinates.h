#pragma once
#include <cmath>

template<class Width>
struct Dimensions
{
    Width x, y;

    Dimensions& operator=(const Dimensions& ref)
    {
        x = ref.x;
        y = ref.y;
    }
    Dimensions& operator=(Width& a)
    {
        x = a;
        y = a;
    }
    bool operator==(const Dimensions& ref)
    {
        return x == ref.x && y == ref.y;
    }
    Dimensions(Width& i, Width& j) : x(i), y(j) {}
    Dimensions() : x(0), y(0) {}
};

template<class T>
struct Coordinates : public Dimensions<T>
{
    Coordinates& operator-(Coordinates& a)
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

    Coordinates(T _x, T _y) : Dimensions(_x, _y) {}
    Coordinates() : Dimensions() {}
};
using Placements = Coordinates<unsigned>;

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