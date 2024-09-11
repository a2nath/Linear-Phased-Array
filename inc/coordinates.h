#pragma once
#include <cmath>
#include <vector>
#include <cstring>

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
