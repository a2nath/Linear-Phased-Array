#pragma once
#include <string>
#include <map>
#include <filesystem>

extern std::string sim_error;

template<class T> std::string str(const T& input) { return std::to_string(input); }

template<class T>
inline T rand(T min, T max)
{
    return min + rand() % (max - min + 1);
}

inline std::string getcwd()
{
    return std::filesystem::current_path().string();
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
