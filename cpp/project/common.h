#pragma once
#include <cmath>

#define C       3e8 //speed of light
#define M_PIl   3.141592653589793238462643383279502884L /* pi */

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
struct Dimension
{
    double x, y;
    Dimension(double i, double j) : x(i), y(j) {}
};
struct Coordinates
{
    int x, y;
    //Coordinates(double i, double j) :x(i), y(j) {}
    //Coordinates() :x(0), y(0) {}
    Coordinates operator-(Coordinates& a)
    {
        Coordinates temp = *this;
        temp.x -= a.x;
        temp.y -= a.y;
        return temp;
    }
};
using cords = Coordinates;

struct Polar_Coordinates
{
    double theta;
    double hype;
    Polar_Coordinates(double i, double j) :theta(i), hype(j) {}
    Polar_Coordinates() : theta(0), hype(0) {}
};

inline Polar_Coordinates cart2pol(double x, double y)
{
    double theta = atan2(y, x);
    double hype = hypot(x, y);
    return Polar_Coordinates(theta, hype);
}

inline Polar_Coordinates cart2pol(Coordinates& c)
{
    return cart2pol(c.x, c.y);
}
