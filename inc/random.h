#pragma once
#include <random>

#define RAND

/* Random number generator module */
class Random
{
    /* First create an instance of an engine */
    std::random_device rnd_device;
    unsigned seed;

    /* Specify the engine and distribution */
    std::mt19937 mersenne_engine;

    /* uniform distribution */
    struct Uniform
    {
        std::uniform_int_distribution<int> getdist(int min, int max)
        {
            return std::uniform_int_distribution<int>{min, max};
        }

        std::uniform_real_distribution<double> getdist(double min, double max)
        {
            return std::uniform_real_distribution<double>{min, max};
        }
    };
public:

    Uniform uni;

    /* Get a random number by passing one distribution obj */
    template<class T>
    inline auto generate(T& dist)
    {
        return dist(mersenne_engine);
    }

    const unsigned& getSeed() const
    {
        return seed;
    }

    Random(unsigned _seed) : seed(_seed), mersenne_engine(rnd_device())
    {
        mersenne_engine.seed(seed);
    }

    Random() : seed(rand()), mersenne_engine(rnd_device())
    {
        mersenne_engine.seed(seed);
    }
};