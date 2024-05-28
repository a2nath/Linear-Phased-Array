#pragma once
#include <SFML/Graphics.hpp>
#include <string>
#include <iostream>
#include <vector>
#include <chrono>
#include <unordered_map>
#include "station.h"
#include "coordinates.h"

using namespace std;
using namespace sf;

constexpr unsigned TPS = 60; //ticks per seconds
const sf::Time     timePerUpdate = sf::seconds(1.0f / float(TPS));

template<class T>
using Pattern2D = vector<vector<T>>;
using Signal2D = Pattern2D<double>;

struct station_details_t
{
    double ant_spacing;
    unsigned ant_panel_count;
    double ant_tx_power;
    double ant_scan_angle;
    Placements position;
};

int compare(const station_details_t& sta1, const station_details_t& sta2)
{
    return (sta1.ant_spacing == sta2.ant_spacing
        && sta1.ant_panel_count == sta2.ant_panel_count
        && sta1.ant_tx_power == sta2.ant_tx_power
        && sta1.ant_scan_angle == sta2.ant_scan_angle
        && sta1.position.x == sta2.position.x
        && sta1.position.y == sta2.position.y) == true ? 0 : 1;
}


namespace draw
{
    template<class T>
    Vector2u shape(const std::vector<std::vector<T>>& grid, bool gettotal = false)
    {
        if (gettotal)
            return { unsigned(grid.size() * grid[0].size()), (unsigned)0 };
        else
            return { (unsigned)grid.size(), (unsigned)grid[0].size() };
    }
}


struct Node
{
    enum class Type
    {
        handset,
        base
    } type;
    CircleShape shape;
    Vector2f station_position;

    static float scale_factor()
    {
        return 5.f;
    }

    CircleShape& operator()()
    {
        return shape;
    }



    Node(
        Type station_type,
        float radius,
        size_t pointCount = 30)
        :
        type(station_type),
        shape(radius, pointCount)
    {
    }
};

struct Tile
{
    RectangleShape shape;
    const Vector2f og_size;

    void reset()
    {
        reshape({ 1, 1 });
    }

    void reshape(const Vector2f& scale_factor)
    {
        if (1 < scale_factor.x && 1 < scale_factor.y)
        {
            shape.setScale(scale_factor);
        }
    }

    static float scale_factor()
    {
        return 100.f;
    }
    Tile(Vector2f size) : og_size(size), shape(size)
    {

    }
};

/* Returns timestamp based on whether it is for a filename. Default false */
inline void getTimeStamp(char* buffer, bool for_filename = false)
{
    auto now = std::chrono::system_clock::now();
    auto now_c = std::chrono::system_clock::to_time_t(now);

    if (for_filename)
        std::strftime(buffer, sizeof(buffer), "%m%d%y_%H%M%S", std::localtime(&now_c));
    else
        std::strftime(buffer, sizeof(buffer), "%m%d%y %H%M%S", std::localtime(&now_c));
}

class Visuals
{
    Vector2u WINDOW_SIZE;
    Vector2f TILE_SIZE;
    RenderWindow window;
    vector<vector<Tile>> tile_grid;
    vector<Node> nodes;

    /* track the state of the sta to trigger the update function */
    std::unordered_map<unsigned, station_details_t> state_data;

    // Function to convert HSV to RGB color
    void hsvToRgb(float h, float s, float v, float& r, float& g, float& b) {
        int i = int(h * 6);
        float f = h * 6 - i;
        float p = v * (1 - s);
        float q = v * (1 - f * s);
        float t = v * (1 - (1 - f) * s);

        switch (i % 6) {
        case 0: r = v, g = t, b = p; break;
        case 1: r = q, g = v, b = p; break;
        case 2: r = p, g = v, b = t; break;
        case 3: r = p, g = q, b = v; break;
        case 4: r = t, g = p, b = v; break;
        case 5: r = v, g = p, b = q; break;
        }
    }

    // Function to map exponential values to HSV color
    void mapToColor(float value, float& h, float& s, float& v) {
        // Normalize value between 0 and 1
        float normalizedValue = std::min(std::max(value, 0.0f), 1.0f);

        // Map value to hue (blue to green to yellow)
        h = 0.66f - normalizedValue * 0.66f;

        // Set saturation and value to maximum
        s = 1.0f;
        v = 1.0f;
    }

    /* change the state of the heat map based on the bs_id */
    //void updateHeatMap(unsigned bs_id)
    //{
    //    draw_init(
    //}

public:
    /* can be signal from one sta or overall snr map */
    void show_eirp(Signal2D& signal_strength)
    {
        TILE_SIZE = { min(WINDOW_SIZE.x, WINDOW_SIZE.y) / Tile::scale_factor(), min(WINDOW_SIZE.x, WINDOW_SIZE.y) / Tile::scale_factor() };
        Tile tile({ TILE_SIZE.x, TILE_SIZE.y });
        tile.shape.setFillColor({ 50, 50, 50, 180 });
        tile.shape.setOutlineThickness(1.0f);
        tile.shape.setOutlineColor({ 50, 50, 50, 250 });
        tile_grid.resize(WINDOW_SIZE.y / TILE_SIZE.y, vector<Tile>(WINDOW_SIZE.x / TILE_SIZE.x, tile));

        for (auto& row : signal_strength)
        {
            for (auto& cell : row)
            {
                float h, s, v;
                mapToColor(cell, h, s, v);

                // Convert HSV to RGB
                float r, g, b;
                hsvToRgb(h, s, v, r, g, b);
            }
        }

        for (auto y = 0; y < tile_grid.size(); ++y)
        {
            for (auto x = 0; x < tile_grid[0].size(); ++x)
            {
                auto& tile = tile_grid[y][x];
                tile.shape.setPosition(x * TILE_SIZE.x, y * TILE_SIZE.y);
            }
        }
    }


    void set_nodes(vector<Placements>& base, vector<Placements>& clients, Dimensions<float> field_size)
    {
        nodes.resize(clients.size(), { Node::Type::handset, TILE_SIZE.x / Node::scale_factor() });
        nodes.resize(clients.size() + base.size(), { Node::Type::base, TILE_SIZE.x / Node::scale_factor() });

        unsigned idx = 0;

        for (auto& station : base)
        {
            auto& node = nodes[idx++];
            node.station_position = { float(station.x / field_size.x), float(station.y / field_size.y) };
            cout << "x " << node.station_position.x << " y " << node.station_position.y << endl;
            node.shape.setPosition(node.station_position.x * WINDOW_SIZE.x, node.station_position.y * WINDOW_SIZE.y);
            node.shape.setFillColor({ 255,255,255 });
            node.shape.setOutlineColor(sf::Color::Transparent);
            node.shape.setOutlineThickness(0.f);
        }

        for (auto& station : clients)
        {
            auto& node = nodes[idx++];
            node.station_position = { float(station.x / field_size.x), float(station.y / field_size.y) };
            cout << "x " << node.station_position.x << " y " << node.station_position.y << endl;
            node.shape.setPosition(node.station_position.x * WINDOW_SIZE.x, node.station_position.y * WINDOW_SIZE.y);
            node.shape.setFillColor({ 150,150,255 });
            node.shape.setOutlineColor(sf::Color::Transparent);
            node.shape.setOutlineThickness(0.f);
        }
    }

    unsigned run()
    {
        window.setPosition(sf::Vector2i{ window.getPosition().x, 0 });

        bool is_focused = false;
        auto view = window.getDefaultView();
        auto speed = 100.f;

        sf::Clock clock;
        sf::Time timeSinceLastUpdate = sf::Time::Zero;
        sf::Time FrameTime = sf::seconds(1.f / 60.f);

        /* background tile when there is no signal */
        Tile bg_tile({ TILE_SIZE.x, TILE_SIZE.y });
        bg_tile.shape.setFillColor({ 50, 50, 50, 180 });
        bg_tile.shape.setOutlineThickness(1.0f);
        bg_tile.shape.setOutlineColor({ 50, 50, 50, 250 });


        while (window.isOpen())
        {
            sf::Time dt = clock.restart();
            timeSinceLastUpdate += dt;

            sf::Event event{};
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed || event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
                    window.close();

                switch (event.type)
                {
                case sf::Event::KeyPressed:
                {
                    switch (event.key.code)
                    {
                    case sf::Keyboard::Enter:
                    {
                        std::cout << "Enter Pressed\n";
                        break;
                    }
                    case sf::Keyboard::Space:
                    {
                        std::cout << "Space Pressed\n";
                        break;
                    }
                    case sf::Keyboard::S:
                    {
                        if (event.key.control)
                        {
                            sf::Texture texture;
                            const auto& size = window.getSize();
                            texture.create(size.x, size.y);
                            texture.update(window);
                            sf::Image screenshot = texture.copyToImage();

                            //char buffer[80]; // set a timestamp to avoid overwriting
                            //auto name = "screenshot_" + string(buffer) + ".png";

                            string name;
                            name.reserve(95);
                            name = "screenshot_";
                            char* buffer = name.data();
                            buffer += name.size();

                            getTimeStamp(buffer, true);
                            name += ".png";

                            // Save the screenshot to a file
                            if (screenshot.saveToFile(name))
                            {
                                std::cout << "Screenshot saved " << name << std::endl;
                            }
                            else
                            {
                                std::cerr << "Failed to save screenshot" << std::endl;
                            }
                        }
                        break;
                    }
                    }
                    break;
                }
                case sf::Event::Resized:
                {
                    window.clear({ Color(100,100,100) }); // called every fame
                    // update the view to the new size of the window
                    ;
                    //cout << "old viewport " << viewport.left << "," << viewport.top << endl;
                    view.setSize(static_cast<float>(event.size.width), static_cast<float>(event.size.height));

                    //cout << " new viewport " << viewport.left << "," << viewport.top << endl;
                    //sf::FloatRect visibleArea(0, 0, event.size.width, event.size.height);

                    //cout << "Window size {" << WINDOW_SIZE.x << "," << WINDOW_SIZE.y << "}, ";
                    //Vector2f scale_size = { float(event.size.width / WINDOW_SIZE.x), float(event.size.height / WINDOW_SIZE.y) };
                    WINDOW_SIZE = { event.size.width, event.size.height };
                    FloatRect visibleArea(0, 0, event.size.width, event.size.height);
                    window.setView(sf::View(visibleArea));
                    view = window.getDefaultView();
                    //cout << "new size {" << WINDOW_SIZE.x << "," << WINDOW_SIZE.y << "}, ";

                    //cout << "tile size {" << TILE_SIZE.x << "," << TILE_SIZE.y << "}, ";

                    cout << "grid size {" << tile_grid.size() << "," << tile_grid[0].size() << "}, ";
                    tile_grid.clear();
                    tile_grid.resize(WINDOW_SIZE.y / TILE_SIZE.y, vector<Tile>(WINDOW_SIZE.x / TILE_SIZE.x, bg_tile));
                    cout << "new size {" << tile_grid.size() << "," << tile_grid[0].size() << "}" << endl;

                    for (auto y = 0; y < tile_grid.size(); ++y)
                    {
                        //cout << "row " << y << ": ";
                        for (auto x = 0; x < tile_grid[0].size(); ++x)
                        {
                            auto& tile = tile_grid[y][x];
                            tile.shape.setPosition(x * TILE_SIZE.x, y * TILE_SIZE.y);
                            //cout << x * TILE_SIZE.x << "," << y * TILE_SIZE.y << ", ";
                            //tile.reshape(scale_size);
                        }
                        //cout << endl;
                    }

                    for (auto& node : nodes)
                    {
                        node.shape.setPosition(node.station_position.x * WINDOW_SIZE.x, node.station_position.y * WINDOW_SIZE.y);
                    }

                    break;
                }
                case sf::Event::LostFocus:
                {
                    is_focused = false;
                    break;
                }
                case sf::Event::GainedFocus:
                {
                    is_focused = true;
                    break;
                }
                }
            }

            if (is_focused)
            {
                auto mousePos = sf::Mouse::getPosition(window);
                auto mouseWorldPos = window.mapPixelToCoords(mousePos, view);

                //window.setTitle("Mouse Position: (" + std::to_string(int(mouseWorldPos.x / 64.f)) + ", " +
                //    std::to_string(int(mouseWorldPos.y / 64.f)) + ")");

                if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                    if (mouseWorldPos.x >= 0 && mouseWorldPos.y >= 0 &&
                        mouseWorldPos.x < WINDOW_SIZE.x && mouseWorldPos.y < WINDOW_SIZE.y)
                    {
                        auto position = mousePos;
                        position.x /= TILE_SIZE.x;
                        position.y /= TILE_SIZE.y;
                        tile_grid[position.y][position.x].shape.setFillColor(sf::Color::Red);
                    }
                }
                //! ** INPUT SECTION **
                //if (sf::Keyboard::isKeyPressed(sf::Keyboard::S)) {
                //    shape.move(0.f, -speed * timePerUpdate.asSeconds());
                //}
                //if (sf::Keyboard::isKeyPressed(sf::Keyboard::D)) {
                //    shape.move(-speed * timePerUpdate.asSeconds(), 0.f);
                //}
                //if (sf::Keyboard::isKeyPressed(sf::Keyboard::W)) {
                //    shape.move(0.f, speed * timePerUpdate.asSeconds());
                //}
                //if (sf::Keyboard::isKeyPressed(sf::Keyboard::A)) {
                //    shape.move(speed * timePerUpdate.asSeconds(), 0.f);
                //}

                while (timeSinceLastUpdate > FrameTime) {
                    timeSinceLastUpdate -= FrameTime;

                    while (timeSinceLastUpdate > FrameTime) {
                        timeSinceLastUpdate -= FrameTime;
                        //! ** UPDATE SECTION**
                    }
                }
            }

            //view.setCenter(shape.getPosition());
            //window.setView(view);

            //! ** DRAW SECTION **
            window.clear(); // called every fame

            for (auto& row : tile_grid)
                for (auto& tile : row)
                    window.draw(tile.shape);

            for (auto& node : nodes)
                window.draw(node.shape);

            window.setView(view);
            window.display();
        }

        return EXIT_SUCCESS;
    }

    static void test()
    {
        unsigned window_width = 200;
        unsigned window_height = 200;

        sf::RenderWindow window(sf::VideoMode(window_width, window_height), "SFML works!");


        sf::CircleShape shape(100.f);
        shape.setFillColor(sf::Color::Green);
        shape.setPosition(10.f, 50.f);

        while (window.isOpen())
        {
            sf::Event event;
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                    window.close();
            }
            window.clear();
            window.draw(shape);
            window.display();
        }
    }

    void draw_init(const std::vector<double>& spacings, const double& lambda, Cow& cow, Signal2D& signal_data)
    {
        auto& current_state = state_data[cow.sid()];

        station_details_t next;
        cow.details(next.ant_spacing, next.ant_panel_count, next.ant_tx_power, next.ant_scan_angle, next.position);

        if (compare(current_state, next)) // if changed, then update
        {
            std::unordered_map<unsigned, std::vector<double>> constants_table; // stores constants calculation

            auto shape = draw::shape(tile_grid, false); // get the total elements and use param(1)
            std::vector<double> signalHeat(shape.x * shape.y);

            cow.setHighResMatrix(spacings[cow.sid()], lambda, constants_table, shape.x, shape.y);
            cow.getHeatMapData(signalHeat);

            signal_data.resize(shape.x, std::vector<double>(shape.y));

            size_t i = 0;
            size_t offset = 0;
            size_t size = signalHeat.size();
            size_t index = 0;

            while (signalHeat.begin() + offset < signalHeat.end())
            {
                signal_data[index++].assign(signalHeat.begin() + offset, signalHeat.begin() + offset + size);
                offset += size;
            }

            current_state = next;
        }
        //signal heat map on grid ready to present
    }

    Visuals(
        unsigned width,
        unsigned height,
        string title
    ) :
    WINDOW_SIZE(width, height),
    window(VideoMode(WINDOW_SIZE.x, WINDOW_SIZE.y), title)
    {
    }
};
