#pragma once

#include <fstream>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.h"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"
#include "rapidjson/istreamwrapper.h"

// Read JSON types here. Simulation rules and unit conversions belong in MyArgs::init().
class JsonConfig
{
    std::string filename;
    rapidjson::Document document;

    static void decode(const rapidjson::Value& value, double& output, const std::string& key)
    {
        if (!value.IsNumber())
        {
            throw std::invalid_argument(key + " must be a number");
        }
        output = value.GetDouble();
    }

    static void decode(const rapidjson::Value& value, unsigned& output, const std::string& key)
    {
        if (!value.IsUint())
        {
            throw std::invalid_argument(key + " must be a nonnegative integer");
        }
        output = value.GetUint();
    }

    static void decode(const rapidjson::Value& value, bool& output, const std::string& key)
    {
        if (!value.IsBool())
        {
            throw std::invalid_argument(key + " must be true or false");
        }
        output = value.GetBool();
    }

    static void decode(const rapidjson::Value& value, Placements& output, const std::string& key)
    {
        if (!value.IsArray() || value.Size() != 2)
        {
            throw std::invalid_argument(key + " must be an [x, y] pair");
        }
        decode(value[0], output.x, key + "[0]");
        decode(value[1], output.y, key + "[1]");
    }

    template<class T>
    static void decode(const rapidjson::Value& value, std::vector<T>& output, const std::string& key)
    {
        if (!value.IsArray())
        {
            throw std::invalid_argument(key + " must be an array");
        }

        output.resize(value.Size());
        for (rapidjson::SizeType index = 0; index < value.Size(); ++index)
        {
            decode(value[index], output[index], key + "[" + std::to_string(index) + "]");
        }
    }

public:

    explicit JsonConfig(const std::string& ifilename) : filename(ifilename)
    {
        std::ifstream stream(filename);
        if (!stream.is_open())
        {
            throw std::runtime_error("Cannot open JSON file: " + filename);
        }

        rapidjson::IStreamWrapper input(stream);
        document.ParseStream(input);

        if (document.HasParseError())
        {
            throw std::runtime_error(filename + ": Invalid JSON format at byte " +
                std::to_string(document.GetErrorOffset()) + ": " +
                rapidjson::GetParseError_En(document.GetParseError()));
        }

        if (!document.IsObject())
        {
            throw std::runtime_error(filename + ": expected a JSON object");
        }
    }

    // The first name is preferred. Remaining names support existing scenario files.
    template<class T>
    void read(std::initializer_list<const char*> names, T& destination, bool required = true) const
    {
        const rapidjson::Value* value = nullptr;
        std::string key;

        for (const char* name : names)
        {
            const auto member = document.FindMember(name);
            if (member == document.MemberEnd())
            {
                continue;
            }

            if (value != nullptr)
            {
                throw std::invalid_argument(filename + ": supply only one spelling of " + *names.begin());
            }

            value = &member->value;
            key = filename + ": " + name;
        }

        if (value != nullptr)
        {
            decode(*value, destination, key);
        }
        else if (required)
        {
            throw std::invalid_argument(filename + ": missing JSON parameter: " + *names.begin());
        }
    }
};
