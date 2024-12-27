#pragma once
#include <iostream>
#include <csignal>
#include "common.h"
#include <deque>
#include <tuple>
#include <future>
//#include <limits>


#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>

#define CONTROLS

#ifdef _DEBUG

//#undef CONTROLS

#endif

#ifdef CONTROLS
#include "imgui.h"
#include "imgui-SFML.h"
#endif

namespace graphics
{
    using namespace std;
    using vertex_v = std::vector<sf::VertexArray>;

    struct txvertex
    {
        //sf::VertexArray transmitter;
        int id;
        vertex_v indicators;
        sf::Vector2i size;
        sf::Vector2i location;
        sf::Color ogcolor;
        sf::RectangleShape transmitter;

        void setPosition(const float& new_x, const float& new_y)
        {
            auto& curloc = transmitter.getPosition();
            sf::Vector2f offset = { new_x - curloc.x, new_y - curloc.y };

            transmitter.setPosition(new_x, new_y);
            for (auto& indicator : indicators)
            {
                const auto& vertices = indicator.getVertexCount();

                for (int v = 0; v < vertices; ++v)
                {
                    indicator[v].position += { offset.x, offset.y };
                }
            }
        }

        void setPosition()
        {
            setPosition(location.x, location.y);
        }

        txvertex(int iid, const Placements& ilocation, sf::Color color)
            :
            id(iid),
            size(10, 10),
            transmitter(sf::Vector2f(size)),
            location(ilocation.x - size.x / 2, ilocation.y - size.y / 2),
            ogcolor(color)
        {
            transmitter.setPosition(sf::Vector2f(location));
            transmitter.setFillColor(ogcolor);
        }
    };

    struct HeatGrid
    {
        const float padding;
        const float offset_height, offset_width;
        const unsigned& data_height;
        const unsigned& data_width;
        unsigned pixel_height, pixel_width;

        float init_pxl_range[2], prev_pxl_range[2], curr_pxl_range[2];
        float init_thresholds[3], prev_thresholds[3], curr_thresholds[3];

        sf::Vector2u bounds_lower, bounds_upper, window_size;

        const placement_v& rx_locations;
        const std::vector<State>& txstates;

        sf::VertexArray grid;
        std::vector<sf::CircleShape> rx_cicles;

        std::vector<txvertex> txdata;
        sf::Font font;

        bool debug_mode;

        //sf::Color monocolorgrid(const double& raw, unsigned& tx_id)
        //{
        //    double color_value = (max_pxl == min_pxl) ? 0.0 : ((raw - min_pxl) / (max_pxl - min_pxl));
        //
        //    // Interpolated Color : Transition from red -> yellow -> green -> blue
        //    if (color_value <= thresholds[0]) {
        //        return sf::Color(tx_id & 0b11 ? 255 : 0, static_cast<sf::Uint8>(color_value * 4 * 255), 255); // Blue to Cyan
        //    }
        //    else if (color_value <= thresholds[1]) {
        //        return sf::Color(0, 255, static_cast<sf::Uint8>((1 - (color_value - thresholds[0]) * 4) * 255)); // Cyan to Green
        //    }
        //    else if (color_value <= thresholds[2]) {
        //        return sf::Color(static_cast<sf::Uint8>((color_value - thresholds[1]) * 4 * 255), 255, 0); // Green to Yellow
        //    }
        //    else {
        //        return sf::Color(255, static_cast<sf::Uint8>((1 - (color_value - thresholds[2]) * 4) * 255), 0); // Yellow to Red
        //    }
        //}


        sf::Color colorgrid(const size_t& tx_id, const double& minval, const double& maxval)
        {
            double color_value = (maxval == minval) ? 0.0 : (double(tx_id - minval) / double(maxval - minval));

            // Interpolated Color : Transition from red -> yellow -> green -> blue
            if (color_value <= curr_thresholds[0])
            {
                return sf::Color(0, static_cast<sf::Uint8>(color_value * 4 * 255), 255); // Blue to Cyan
            }
            else if (color_value <= curr_thresholds[1])
            {
                return sf::Color(0, 255, static_cast<sf::Uint8>((1 - (color_value - curr_thresholds[0]) * 4) * 255)); // Cyan to Green
            }
            else if (color_value <= curr_thresholds[2])
            {
                return sf::Color(static_cast<sf::Uint8>((color_value - curr_thresholds[1]) * 4 * 255), 255, 0); // Green to Yellow
            }
            else
            {
                return sf::Color(255, static_cast<sf::Uint8>((1 - (color_value - curr_thresholds[2]) * 4) * 255), 0); // Yellow to Red
            }
        }

        sf::Color colorgrid(const float& raw, const float& minval, const float& maxval)
        {
            float color_value = (maxval == minval) ? 0.0 : ((raw - minval) / (maxval - minval));

            // Interpolated Color : Transition from red -> yellow -> green -> blue
            sf::Uint8 r, g, b;
            if (color_value > curr_thresholds[2])  // Yellow to Red
            {
                r = 255; g = static_cast<sf::Uint8>((1.0 - color_value) / (1.0 - curr_thresholds[2]) * 255); b = 0;
            }
            else if (color_value > curr_thresholds[1]) // Green to Yellow
            {
                r = static_cast<sf::Uint8>((color_value - curr_thresholds[1]) / (curr_thresholds[2] - curr_thresholds[1]) * 255); g = 255; b = 0;
            }
            else if (color_value > curr_thresholds[0]) // Cyan to Green
            {
                r = 0; g = 255; b = static_cast<sf::Uint8>((curr_thresholds[1] - color_value) / (curr_thresholds[1] - curr_thresholds[0]) * 255);
            }
            else // Blue to Cyan
            {
                r = 0; g = static_cast<sf::Uint8>(color_value / curr_thresholds[0] * 255); b = 255;
            }

            return sf::Color(r, g, b);
        }

        void draw_legend(sf::RenderWindow& window) {
            // Legend dimensions
            const float legendWidth = 55.f;
            const float legendHeight = 400.f;
            //const float bar_height = legendHeight / 50.f;
            const sf::Vector2f legendPos(data_width + offset_width + 30.f, offset_height); // Place on the right side

            // Create the legend gradient
            float range_normalized = (curr_pxl_range[1] - curr_pxl_range[0]) / static_cast<float>(legendHeight - 1);
            float height_offset = legendPos.y + legendHeight;

            sf::RectangleShape gradientRect(sf::Vector2f(legendWidth, 1.f));
            for (int i = 0; i < legendHeight; ++i)
            {

                float signalValue = curr_pxl_range[0] + i * range_normalized;
                //sf::Color color = ;
                gradientRect.setFillColor(colorgrid(signalValue, curr_pxl_range[0], curr_pxl_range[1]));

                //std::cout << "[" << i << "] signal:" << signalValue << " colorvalue:" << ((max_pxl == min_pxl) ? 0.0 : ((signalValue - min_pxl) / (max_pxl - min_pxl))) \
                //   << " color:" << to_string(color) << " is at height:" << height_offset - i << " maxval:" << max_pxl << " minval:" << min_pxl << std::endl;

                gradientRect.setPosition(legendPos.x, height_offset - i);
                window.draw(gradientRect);
            }

            // Labels for the legend

            sf::Text labelMin, labelMax, labelMid;

            labelMax.setFont(font);
            labelMax.setString(str(static_cast<int>(curr_pxl_range[1])) + " dB" + (debug_mode ? "m" : ""));
            labelMax.setCharacterSize(15);
            labelMax.setFillColor(sf::Color::White);
            labelMax.setPosition(legendPos.x + legendWidth + 5.f, legendPos.y - 10.f);
            window.draw(labelMax);

            labelMid.setFont(font);
            labelMid.setCharacterSize(15);
            labelMid.setFillColor(sf::Color::White);
            int steps = legendHeight / 100.0;
            float range = curr_pxl_range[0] + curr_pxl_range[1];

            for (int i = 1; i < steps; ++i)
            {
                labelMid.setString(str(static_cast<int>(i * range / steps)) + " dB" + (debug_mode ? "m" : ""));
                labelMid.setPosition(legendPos.x + legendWidth + 5.f, legendPos.y + i * 100);
                window.draw(labelMid);
            }

            labelMin.setFont(font);
            labelMin.setString(str(static_cast<int>(curr_pxl_range[0])) + " dB" + (debug_mode ? "m" : ""));
            labelMin.setCharacterSize(15);
            labelMin.setFillColor(sf::Color::White);
            labelMin.setPosition(legendPos.x + legendWidth + 5.f, legendPos.y + legendHeight - 10.f);
            window.draw(labelMin);
        }

        void update_panning(sf::Vector2i& moved_offset)
        {
            bounds_lower = { (unsigned)max((float)moved_offset.x, (float)0.0), (unsigned)max((float)moved_offset.y, (float)0.0) };
            bounds_upper = { (unsigned)min((float)window_size.x, (float)window_size.x + moved_offset.x), (unsigned)min((float)window_size.y, (float)window_size.y + moved_offset.y) };
            moved_offset = { 0, 0 };
        }

        /* update the heat colors from raw calculations */
        void update_heat(const double_v& tx_raw_sigdata)
        {
            long long index = 0;
            for (long long row_idx = data_height - 1; row_idx >= 0; --row_idx)
            {
                auto row_offset = row_idx * data_width;
                for (long long cols_idx = 0; cols_idx < data_width; ++cols_idx)
                {
                    size_t v_index = (row_offset + cols_idx) << 2;

                    sf::Color v_color = colorgrid(tx_raw_sigdata[index++], curr_pxl_range[0], curr_pxl_range[1]);
                    grid[v_index + 0].color = v_color;
                    grid[v_index + 1].color = v_color;
                    grid[v_index + 2].color = v_color;
                    grid[v_index + 3].color = v_color;
                }
            }
        }

        inline void resize(const Dimensions<unsigned>& shape)
        {
            render_space = shape;
        }

        void tx_lines(const float& dir_radians, txvertex& txvertex_obj)
        {
            txvertex_obj.indicators.clear();

            auto& position = txvertex_obj.location;
            auto& idx = txvertex_obj.id;

            sf::VertexArray line11(sf::Lines, 2), line12(sf::Lines, 2);

            // Calculate the endpoint of the first line based on direction
            float length = 30.0f / 2;  // Halfway across the rectangle

            // Line 1's position: calculated from the antenna direction
            float line1StartX = position.x + 5.0f; // Middle of the rectangle;
            float line1StartY = position.y + 5.0f;

            // Endpoint using the direction to calculate the x and y offset
            float line11EndX = line1StartX - length * cos(dir_radians);
            float line11EndY = line1StartY - length * sin(dir_radians);

            float line12EndX = line1StartX + length * cos(dir_radians);
            float line12EndY = line1StartY + length * sin(dir_radians);

            line11[0].position = sf::Vector2f(line1StartX, line1StartY);
            line11[1].position = sf::Vector2f(line11EndX, line11EndY);

            line11[0].color = sf::Color::White;
            line11[1].color = sf::Color::White;

            line12[0].position = sf::Vector2f(line1StartX, line1StartY);
            line12[1].position = sf::Vector2f(line12EndX, line12EndY);

            line12[0].color = sf::Color::White;
            line12[1].color = sf::Color::White;

            txvertex_obj.indicators.emplace_back(line11);
            txvertex_obj.indicators.emplace_back(line12);

            /* draw the scan angle of the linear phase array */
            sf::VertexArray arrow(sf::Lines, 6);

            // Calculate the beam angle relative to the first line (add scan_angle to direction)
            float beam_angle_radians = txstates[idx].settings.theta_c + txstates[idx].settings.alpha;

            float beam_length = 25.0f;  // Length of the signal beam line

            // Line 1's position: calculated from the antenna direction
            float line2StartX = position.x + 5.0f; // Middle of the rectangle
            float line2StartY = position.y + 5.0f;

            // Calculate the endpoint of the second line based on beam angle
            float line2EndX = line2StartX + beam_length * cos(beam_angle_radians);
            float line2EndY = line2StartY - beam_length * sin(beam_angle_radians);

            arrow[0].position = sf::Vector2f(line2StartX, line2StartY);  // Starts from the end of the first line
            arrow[1].position = sf::Vector2f(line2EndX, line2EndY);

            arrow[0].color = sf::Color::Blue;
            arrow[1].color = sf::Color::Blue;

            // Arrowhead parameters
            float arrowhead_angle = M_PIl / 6;  // 30 degrees for the arrowhead angle
            float arrowhead_length = 5.0f;    // Length of the arrowhead sides

            // Calculate the points for the arrowhead
            float left_head_angle = beam_angle_radians + arrowhead_angle;   // Angle for the left arrowhead
            float right_head_angle = beam_angle_radians - arrowhead_angle;  // Angle for the right arrowhead

            float left_head_x = line2EndX - arrowhead_length * cos(left_head_angle);
            float left_head_y = line2EndY + arrowhead_length * sin(left_head_angle);

            float right_head_x = line2EndX - arrowhead_length * cos(right_head_angle);
            float right_head_y = line2EndY + arrowhead_length * sin(right_head_angle);

            // Draw the left side of the arrowhead
            arrow[2].position = sf::Vector2f(line2EndX, line2EndY);  // From the beam endpoint
            arrow[3].position = sf::Vector2f(left_head_x, left_head_y);  // To the left side of the arrowhead
            arrow[2].color = sf::Color::Blue;
            arrow[3].color = sf::Color::Blue;

            // Draw the right side of the arrowhead
            arrow[4].position = sf::Vector2f(line2EndX, line2EndY);  // From the beam endpoint
            arrow[5].position = sf::Vector2f(right_head_x, right_head_y);  // To the right side of the arrowhead
            arrow[4].color = sf::Color::Blue;
            arrow[5].color = sf::Color::Blue;

            txvertex_obj.indicators.emplace_back(arrow);
        }

        void init(const sf::Vector2u& ibounds_lower, const sf::Vector2u& ibounds_upper)
        {
            bounds_lower = ibounds_lower;
            bounds_upper = ibounds_upper;

            size_t index = 0;
            for (long long row_idx = data_height - 1; row_idx >= 0; --row_idx)
            {
                auto row_offset = row_idx * data_width;
                for (long long cols_idx = 0; cols_idx < data_width; ++cols_idx)
                {
                    size_t v_index = (row_offset + cols_idx) << 2;

                    grid[v_index + 0].position = sf::Vector2f(offset_width + cols_idx * pixel_width, offset_height + row_idx * pixel_height);
                    grid[v_index + 1].position = sf::Vector2f(offset_width + (cols_idx + 1) * pixel_width, offset_height + row_idx * pixel_height);
                    grid[v_index + 2].position = sf::Vector2f(offset_width + (cols_idx + 1) * pixel_width, offset_height + (row_idx + 1) * pixel_height);
                    grid[v_index + 3].position = sf::Vector2f(offset_width + cols_idx * pixel_width, offset_height + (row_idx + 1) * pixel_height);
                }
            }

            sf::CircleShape sta(8.0);
            sta.setFillColor(sf::Color(90, 90, 90));
            sta.setOutlineColor(sf::Color::Black);
            sta.setOutlineThickness(2.0f);

            for (auto& loc : rx_locations)
            {
                sta.setPosition(loc.x - 8.0 + offset_width, data_height - 8.0 - loc.y + offset_height);
                rx_cicles.emplace_back(sta);
            }

            for (int i = 0; i < txstates.size(); ++i)
            {
                auto& loc = txstates[i].location;
                float dir_radians = M_PIl / 2 - txstates[i].settings.theta_c;

                txdata.emplace_back(i, Placements{ loc.x + (unsigned)offset_width, unsigned(data_height + offset_height) - loc.y }, sf::Color(90, 90, 90));

                /* draw the placement direction indicators for each tower */
                tx_lines(dir_radians, txdata.back());
            }
        }

        Placements bounded(const int& x, const int& y)
        {
            Placements p;
            p.x = min((int)data_width, x);
            p.x = max(0, x);

            p.y = min((int)data_height, y);
            p.y = max(0, y);

            return p;
        }

        void reset()
        {
            std::copy(std::begin(curr_thresholds), std::end(curr_thresholds), std::begin(prev_thresholds));
            std::copy(std::begin(init_thresholds), std::end(init_thresholds), std::begin(curr_thresholds));

            std::copy(std::begin(curr_pxl_range), std::end(curr_pxl_range), std::begin(prev_pxl_range));
            std::copy(std::begin(init_pxl_range), std::end(init_pxl_range), std::begin(curr_pxl_range));
        }

        /* undo the action performed */
        void undo()
        {
            std::swap(prev_thresholds, curr_thresholds);
            std::swap(prev_pxl_range, curr_pxl_range);
        }

        HeatGrid(
            const unsigned& icols,
            const unsigned& irows,
            const float& imin,
            const float& imax,
            const sf::Vector2u& iwindow_size,
            const placement_v& irx_locations,
            const vector<State>& curr_state)
            :
            padding(30.0f),
            offset_height(padding),
            offset_width(padding),
            data_width(icols),
            data_height(irows),
            pixel_width(iwindow_size.x / icols),
            pixel_height(iwindow_size.y / irows),
            bounds_lower({0, 0}),
            bounds_upper(iwindow_size),
            rx_locations(irx_locations),
            txstates(curr_state),
            grid(sf::Quads, irows * icols * 4),
            debug_mode(false)
        {
            init_thresholds[0] = 0.25; // Cyan to Green
            init_thresholds[1] = 0.50; // Green to Yellow
            init_thresholds[2] = 0.81; // Yellow to Red

            init_pxl_range[0] = imin;
            init_pxl_range[1] = imax;

            std::copy(std::begin(init_thresholds), std::end(init_thresholds), std::begin(prev_thresholds));
            std::copy(std::begin(init_thresholds), std::end(init_thresholds), std::begin(curr_thresholds));

            std::copy(std::begin(init_pxl_range), std::end(init_pxl_range), std::begin(prev_pxl_range));
            std::copy(std::begin(init_pxl_range), std::end(init_pxl_range), std::begin(curr_pxl_range));

            init(bounds_lower, bounds_upper);

            if (!font.loadFromFile("font/OpenSans-Light.ttf"))
            {
                std::cout << "current dir:" << system("cd") << endl;
                spdlog::warn("Font file did not load for SFML library");
            }

            resize({ data_width, data_height });
        }
    };

    // Function to handle Ctrl + C
    void sig_handler(int signal)
    {
        if (signal == SIGINT)
        {
            std::cout << "Ctrl + C pressed, exiting...\n";
            std::exit(0);
        }
    }

    /* write to a file: loggerfile, string, debug-print */
    void printerr(Logger& logger, const string& str, bool debug = true)
    {
        logger.write(str);
        if (debug)
            cerr << str << endl;
    }

    void print(Logger& logger, const string& str, bool debug = true)
    {
        logger.write(str);
        if (debug)
            cout << str << endl;
    }

    inline void zoom_in(sf::RenderWindow& window, sf::View& view, float& zoom, const float& adjustent)
    {
        zoom /= adjustent;
        view.setSize(window.getDefaultView().getSize());
        view.zoom(zoom);
        window.setView(view);
    }

    inline void zoom_out(sf::RenderWindow& window, sf::View& view, float& zoom, const float& adjustent)
    {
        zoom *= adjustent;
        view.setSize(window.getDefaultView().getSize());
        view.zoom(zoom);
        window.setView(view);
    }

    inline void pan_window(sf::RenderWindow& window, sf::View& view, sf::Vector2f& curr_pos, const sf::Vector2i& offset)
    {
        sf::Vector2f delta = { (float)offset.x, (float)offset.y };
        view.move(delta);  // Move the view by the delta
        window.setView(view);
        curr_pos += delta;
    }

    /*    G U I    */
    int render(
        Logger& logger,
        const placement_v& init_rx_locations,
        const state_v& tx_states,
        const std::vector<double_v>& raw_cow_data,
        const std::vector<double_v>& mrg_cow_data,
        const unsigned& grid_rows,
        const unsigned& grid_cols,
        const float& min_color_span,
        const float& max_color_span,
        DataSync& sync,
        bool& is_rendering)
    {
        std::signal(SIGINT, sig_handler);

        size_t render_width = grid_cols + 800;
        size_t render_height = grid_rows + 500;

        sf::RenderWindow window(sf::VideoMode(render_width, render_height), "SFML Grid Plot");
        window.setVerticalSyncEnabled(true);

        auto window_size = window.getSize();

        sf::Vector2f curr_position;
        sf::Vector2i panning_view;
        sf::Vector2i moved_offset;

        sf::Clock delta_clock;
        sf::View view = window.getDefaultView();

        /* init variables */
        bool state_changed = true;
        bool panning = false;
        float zoomLevel = 1.0f;
        float mouse_delta_thresh = 0.01f;
        float zoom_change_factor = 1.1f;
        long pan_adj_factor = 10;
        int render_tx_id = 0;
        size_t tx_count = tx_states.size();

        /* set the states, current (working variable), previous (for undo), init (for reset) */
        auto& init = tx_states;
        state_v curr, prev;
        txvertex* tx_dragging = nullptr;

        /* heat data contains vertices too */
        HeatGrid griddata(grid_cols, grid_rows, min_color_span, max_color_span, window_size, init_rx_locations, init);

        std::vector<std::string> tx_header, tx_x_slider, tx_x_inp, tx_y_slider, tx_y_inp, tx_power_slider, tx_power_inp, tx_dir_slider, tx_dir_inp, tx_scan_slider, tx_scan_inp;
        std::vector<Coordinates<int>> grid_tx_offsets;
        std::vector<float> power_dBm(init.size());

        for (auto i = 0; i < init.size(); ++i)
        {

            grid_tx_offsets.emplace_back(griddata.txdata[i].location.x - init[i].location.x, \
                                        griddata.txdata[i].location.y - init[i].location.y);

            auto sidx = str(i);
            tx_header.emplace_back("TX " + sidx);
            tx_x_slider.emplace_back("X##slider" + sidx);
            tx_x_inp.emplace_back("X##input" + sidx);
            tx_y_slider.emplace_back("Y##slider" + sidx);
            tx_y_inp.emplace_back("Y##input" + sidx);
            tx_power_slider.emplace_back("Power##slider" + sidx);
            tx_power_inp.emplace_back("Power##input" + sidx);
            tx_dir_slider.emplace_back("Direction##slider" + sidx);
            tx_dir_inp.emplace_back("Direction##input" + sidx);
            tx_scan_slider.emplace_back("Scan Angle##slider" + sidx);
            tx_scan_inp.emplace_back("Scan Angle##input" + sidx);

            power_dBm[i] = cached::watt2dBm(init[i].settings.power);
        }

        curr = init;
        prev = init;

        /* init the heatmap to display heat from TX id */
        const std::vector<std::vector<double>>* ptr_live_data = &mrg_cow_data;
        sync.render_tx_id = render_tx_id;

#ifdef CONTROLS
        ImGui::SFML::Init(window);
#endif
        /* only update the heat when INIT or moving MOVING tx on the map */
        std::thread heat_checker([&]()
            {
                while (is_rendering)
                {
                    std::unique_lock<std::mutex> lock(graphics::graphics_data_mutex);  // Lock the mutex
                    graphics::consig.wait(lock, [&]()
                        {
                            return sync.render_tx_id >= 0 || !is_rendering;
                        }
                    );

                    if (!is_rendering) break;

                    griddata.update_heat((*ptr_live_data)[sync.render_tx_id]);
                    sync.render_tx_id = -1;
                }
            }
        );

        // Main loop
        while (window.isOpen() && is_rendering)
        {
            sf::Event event;

            while (window.pollEvent(event))
            {
                // OpenGL rendering here
                //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
#ifdef CONTROLS
                ImGui::SFML::ProcessEvent(event);
#endif
                switch (event.type)
                {
                case sf::Event::Closed:
                {
                    window.close();
                    break;
                }
                case sf::Event::Resized:
                {
                    auto new_size = window.getSize();
                    griddata.resize({ new_size.x, new_size.y });
                    sync.emplace_resize(new_size.x, new_size.y, render_tx_id);
                    break;
                }
                case sf::Event::MouseWheelScrolled:
                {
                    if (event.mouseWheelScroll.delta - mouse_delta_thresh > 0 && !ImGui::IsWindowHovered())
                    {
                        zoom_in(window, view, zoomLevel, zoom_change_factor);
                    }
                    else if (event.mouseWheelScroll.delta - mouse_delta_thresh < 0 && !ImGui::IsWindowHovered())
                    {
                        zoom_out(window, view, zoomLevel, zoom_change_factor);
                    }
                    break;
                }
                case sf::Event::MouseButtonPressed:
                {
                    // Mouse press: check if the click was inside the object
                    switch (event.mouseButton.button)// == sf::Mouse::Left)
                    {
                    case sf::Mouse::Left:
                    {
                        auto& txvertex = griddata.txdata;
                        for (auto& tx : txvertex)
                        {
                            if (tx.transmitter.getGlobalBounds().contains(event.mouseButton.x, event.mouseButton.y))
                            {   // found it
                                tx_dragging = &tx;

                                auto& id = tx_dragging->id;
                                prev[id].location = { (unsigned)tx_dragging->location.x - grid_tx_offsets[id].x, \
                                                    (unsigned)tx_dragging->location.y - grid_tx_offsets[id].y };

                                sync.render_tx_id = tx_dragging->id;
                                break;
                            }
                        }

                        break;
                    }
                    case sf::Mouse::Right:
                    {
                        panning_view = sf::Mouse::getPosition(window);
                        panning = true;
                        break;
                    }
                    default:
                        break;
                    }

                    break;
                }
                case sf::Event::MouseButtonReleased:
                {
                    // Mouse press: check if the click was inside the object
                    switch (event.mouseButton.button)
                    {
                    case sf::Mouse::Left:
                    {
                        if (tx_dragging)
                        {
                            /* mouse release causes heat update */
                            sync.render_tx_id = tx_dragging->id;
                            tx_dragging = nullptr;
                        }
                        break;
                    }
                    case sf::Mouse::Right:
                    {
                        griddata.update_panning(moved_offset);

                        panning = false;
                        break;
                    }
                    default:
                        break;
                    } // end switch

                    break;
                }
                case sf::Event::MouseMoved:
                {
                    if (tx_dragging)
                    {
                        auto new_loc = griddata.bounded(event.mouseMove.x, event.mouseMove.y);
                        auto& id = tx_dragging->id;

                        curr[id].location = { new_loc.x - grid_tx_offsets[id].x, new_loc.y - grid_tx_offsets[id].y };
                        tx_dragging->setPosition(event.mouseMove.x, event.mouseMove.y);

                        sync.emplace_state(curr[id]);
                    }

                    if (panning)
                    {
                        auto new_view = sf::Mouse::getPosition(window);
                        moved_offset = new_view - panning_view;

                        if (abs(panning_view.x - mouse_delta_thresh) > 0 || abs(panning_view.y - mouse_delta_thresh) > 0)
                        {
                            pan_window(window, view, curr_position, moved_offset);
                            panning_view = new_view;
                        }
                    }

                    break;
                }
                case sf::Event::KeyPressed:
                {
                    switch (event.key.scancode)
                    {
                    case sf::Keyboard::Scan::Z:
                    {
                        if (event.key.control)
                        {
                            for (int i = 0; i < tx_count; ++i)
                            {
                                if (prev[i] != curr[i])
                                {
                                    auto& tx = curr[i];

                                    griddata.txdata[tx.tx_idx].setPosition(\
                                        tx.location.x + grid_tx_offsets[tx.tx_idx].x, \
                                        tx.location.y + grid_tx_offsets[tx.tx_idx].y);

                                    sync.emplace_state(tx);
                                }

                                curr[i] = prev[i];
                            }

                            if (sync.mainq.empty())
                            {
                                sync.render_tx_id = render_tx_id;
                            }

                            /* undo the color thresholds */
                            griddata.undo();
                        }
                        break;
                    }
                    case sf::Keyboard::Scan::D:
                    {
                        if (event.key.control)
                        {
                            if (griddata.debug_mode == false)
                            {
                                ptr_live_data = &raw_cow_data;
                                griddata.debug_mode = true;
                            }
                            else
                            {
                                ptr_live_data = &mrg_cow_data;
                                griddata.debug_mode = false;
                            }

                            sync.render_tx_id = render_tx_id;
                        }
                        break;
                    }
                    case sf::Keyboard::Scan::R:
                    {
                        zoomLevel = 1.0f;  // Reset zoom level
                        view = window.getDefaultView();

                        if (event.key.shift) // Ctrl + R - erases all changes
                        {
                            for (int i = 0; i < curr.size(); ++i)
                            {
                                if (init[i] != curr[i])
                                {
                                    prev[i] = curr[i];
                                    curr[i] = init[i];
                                    auto& tx = curr[i];

                                    griddata.txdata[tx.tx_idx].setPosition(\
                                        tx.location.x + grid_tx_offsets[tx.tx_idx].x, \
                                        tx.location.y + grid_tx_offsets[tx.tx_idx].y);

                                    sync.emplace_state(tx);
                                }
                            }
                        }
                        else
                        {
                            sync.render_tx_id = render_tx_id;
                        }

                        griddata.reset();

                        window.setView(view);
                        break;
                    }
                    case sf::Keyboard::Scan::Tab:
                    {
                        render_tx_id = (render_tx_id + 1) % tx_count;
                        sync.render_tx_id = render_tx_id;
                        break;
                    }
                    case sf::Keyboard::Scan::Left:
                    {
                        moved_offset = { 0 + pan_adj_factor, 0 };
                        pan_window(window, view, curr_position, moved_offset);
                        griddata.update_panning(moved_offset);
                        break;
                    }
                    case sf::Keyboard::Scan::Up:
                    {
                        moved_offset = { 0, 0 + pan_adj_factor };
                        pan_window(window, view, curr_position, moved_offset);
                        griddata.update_panning(moved_offset);
                        break;
                    }
                    case sf::Keyboard::Scan::Right:
                    {
                        moved_offset = { -pan_adj_factor, 0 };
                        pan_window(window, view, curr_position, moved_offset);
                        griddata.update_panning(moved_offset);
                        break;
                    }
                    case sf::Keyboard::Scan::Down:
                    {
                        moved_offset = { 0, 0 - pan_adj_factor };
                        pan_window(window, view, curr_position, moved_offset);
                        griddata.update_panning(moved_offset);
                        break;
                    }
                    }
                }
                default:
                    break;
                }

            } // end while pollEvent

#ifdef CONTROLS
            // ImGui logic for GUI buttons
            ImGui::SFML::Update(window, delta_clock.restart());
            ImGui::Begin("Control Panel");

            if (ImGui::CollapsingHeader("Map Control", ImGuiTreeNodeFlags_DefaultOpen))
            {
                if (ImGui::Button("Zoom In"))
                {
                    zoom_in(window, view, zoomLevel, zoom_change_factor);
                }
                else if (ImGui::Button("Zoom Out"))
                {
                    zoom_out(window, view, zoomLevel, zoom_change_factor);
                }
                else if (ImGui::Button("Reset View"))
                {
                    zoomLevel = 1.0f;  // Reset zoom level
                    window.setView(window.getDefaultView());
                }
                else if (ImGui::Button("Left"))
                {
                    moved_offset = { 0 + pan_adj_factor, 0 };
                    pan_window(window, view, curr_position, moved_offset);
                    griddata.update_panning(moved_offset);
                }
                else if (ImGui::Button("Up"))
                {
                    moved_offset = { 0, 0 + pan_adj_factor };
                    pan_window(window, view, curr_position, moved_offset);
                    griddata.update_panning(moved_offset);

                }
                else if (ImGui::Button("Right"))
                {
                    moved_offset = { -pan_adj_factor, 0 };
                    pan_window(window, view, curr_position, moved_offset);
                    griddata.update_panning(moved_offset);

                }
                else if (ImGui::Button("Down"))
                {
                    moved_offset = { 0, 0 - pan_adj_factor };
                    pan_window(window, view, curr_position, moved_offset);
                    griddata.update_panning(moved_offset);
                }
                else if (ImGui::Button("Debug"))
                {
                    if (griddata.debug_mode == false)
                    {
                        ptr_live_data = &raw_cow_data;
                        griddata.debug_mode = true;
                    }
                    else
                    {
                        ptr_live_data = &mrg_cow_data;
                        griddata.debug_mode = false;
                    }

                    sync.render_tx_id = render_tx_id;
                }
            }

            if (ImGui::CollapsingHeader("Heatmap Control", ImGuiTreeNodeFlags_DefaultOpen))
            {
                /* Signal Threshold Sliders */
                ImGui::Text("Signal Threshold Control");

                if (ImGui::SliderFloat("Low##slider", &griddata.curr_thresholds[0], 0.0f, 1.0f, "%.2f") ||
                    ImGui::InputFloat("Low##input", &griddata.curr_thresholds[0], 0.1f, 1.0f, "%.2f"))
                    sync.render_tx_id = render_tx_id;

                if (ImGui::SliderFloat("Mid##slider", &griddata.curr_thresholds[1], 0.0f, 1.0f, "%.2f") ||
                    ImGui::InputFloat("Mid##input", &griddata.curr_thresholds[1], 0.1f, 1.0f, "%.2f"))
                    sync.render_tx_id = render_tx_id;

                if (ImGui::SliderFloat("High##slider", &griddata.curr_thresholds[2], 0.0f, 1.0f, "%.2f") ||
                    ImGui::InputFloat("High##input", &griddata.curr_thresholds[2], 0.1f, 1.0f, "%.2f"))
                    sync.render_tx_id = render_tx_id;


                // Ensure thresholds are in the correct order
                if (griddata.curr_thresholds[0] > griddata.curr_thresholds[1])
                    std::swap(griddata.curr_thresholds[0], griddata.curr_thresholds[1]);

                if (griddata.curr_thresholds[1] > griddata.curr_thresholds[2])
                    std::swap(griddata.curr_thresholds[1], griddata.curr_thresholds[2]);

                if (griddata.curr_thresholds[0] > griddata.curr_thresholds[1])
                    std::swap(griddata.curr_thresholds[0], griddata.curr_thresholds[1]);


                /* Add sliders for minand max values */
                ImGui::Text("Signal Min and Max");
                if (ImGui::SliderFloat("Min##slider", &griddata.curr_pxl_range[0], -300.0f, 50.0f, "%.2f dB") ||
                    ImGui::InputFloat("Min##input", &griddata.curr_pxl_range[0], -300.0f, 50.0f, "%.2f"))
                    sync.render_tx_id = render_tx_id;

                if (ImGui::SliderFloat("Max##slider", &griddata.curr_pxl_range[1], -299.0f, 50.0f, "%.2f dB") ||
                    ImGui::InputFloat("Max##input", &griddata.curr_pxl_range[1], -299.0f, 50.0f, "%.2f"))
                    sync.render_tx_id = render_tx_id;

                // Ensure minval is always less than maxval
                if (griddata.curr_pxl_range[0] >= griddata.curr_pxl_range[1])
                {
                    griddata.curr_pxl_range[0] = griddata.curr_pxl_range[1] - 1.0f;
                }
            }

            if (ImGui::CollapsingHeader("Transmitter Control", ImGuiTreeNodeFlags_DefaultOpen))
            {
                for (auto i = 0; i < curr.size(); ++i)
                {
                    if (ImGui::CollapsingHeader(tx_header[i].c_str(), ImGuiTreeNodeFlags_DefaultOpen))
                    {
                        /* position mechanism for each transmitter */
                        ImGui::Text("Placement");
                        auto& position = griddata.txdata[i].location;

                        if (ImGui::SliderInt(tx_x_slider[i].c_str(), &position.x, 0, grid_cols - 1, "%d meters") || // unsigned long -> int?
                            ImGui::InputInt(tx_x_inp[i].c_str(), &position.x, 0, grid_cols - 1))
                        {
                            curr[i].location.x = position.x - grid_tx_offsets[i].x;
                            griddata.txdata[i].setPosition();
                            sync.emplace_state(curr[i]);
                        }

                        if (ImGui::SliderInt(tx_y_slider[i].c_str(), &position.y, 0, grid_rows - 1, "%d meters") ||
                            ImGui::InputInt(tx_y_inp[i].c_str(), &position.y, 0, grid_rows - 1))
                        {
                            curr[i].location.y = position.y - grid_tx_offsets[i].y;
                            griddata.txdata[i].setPosition();
                            sync.emplace_state(curr[i]);
                        }

                        /* antenna-power mechanism for each transmitter */
                        ImGui::Text("Antenna Power");

                        if (ImGui::SliderFloat(tx_power_slider[i].c_str(), &power_dBm[i], -30.0f, +30.0f, "%.2f dBm"))
                        {
                            curr[i].settings.power = cached::dBm2watt(power_dBm[i]);
                            sync.emplace_state(curr[i]);
                        }
                        ImGui::SameLine();
                        if (ImGui::InputFloat(tx_power_inp[i].c_str(), &power_dBm[i], -30.0f, +30.0f, "%.2f"))
                        {
                            curr[i].settings.power = cached::dBm2watt(power_dBm[i]);
                            sync.emplace_state(curr[i]);
                        }

                        /* antenna-power mechanism for each transmitter */
                        ImGui::Text("Antenna Direction");

                        if (ImGui::SliderAngle (tx_dir_slider[i].c_str(), &curr[i].settings.theta_c, 0.0f, 359.9f, "%.2f deg") ||
                            ImGui::InputFloat(tx_dir_inp[i].c_str(), &curr[i].settings.theta_c, 0.0f, 359.9f, "%.2f"))
                        {
                            /* draw the placement direction indicators for each tower */
                            griddata.tx_lines(M_PIl / 2 - curr[i].settings.theta_c, griddata.txdata[i]);
                            sync.emplace_state(curr[i]);
                        }

                        /* antenna-scan angle for each transmitter */
                        ImGui::Text("Scan Angle");

                        if (ImGui::SliderAngle(tx_scan_slider[i].c_str(), &curr[i].settings.alpha, -90.0, +90.0, "%.2f deg") ||
                            ImGui::InputFloat(tx_scan_inp[i].c_str(), &curr[i].settings.alpha, -90.0, +90.0, "%.2f"))
                        {
                            sync.emplace_state(curr[i]);
                        }
                    }
                }
            }

            if (sync.got_updates(render_tx_id))
            {
                consig.notify_one(); // either have [render_tx_id] set or [compute_tx_id] set, not both
            }

            ImGui::End();
#endif

            window.clear();

            /* draw the grid */
            window.draw(griddata.grid);

            /* draw the receivers */
            for (auto& circle : griddata.rx_cicles)
            {
                window.draw(circle);
            }

            /* draw the transmitters */
            for (auto& data : griddata.txdata)
            {
                window.draw(data.transmitter);
                for (auto& element : data.indicators)
                {
                    window.draw(element);
                }
            }

            /* create a legend */
            griddata.draw_legend(window);

#ifdef CONTROLS
            ImGui::SFML::Render(window);  // Render ImGui over SFML content
#endif
            window.display();
        }

#ifdef CONTROLS
        ImGui::SFML::Shutdown();
#endif

        is_rendering = false; // Set rendering to false
        consig.notify_all(); // Notify all threads waiting on the condition variable

        heat_checker.join();

        return 0;
    }

    /*    G U I    */
    void capture_plot(
        Logger& logger,
        const string& filename,
        const state_v& tx_states,
        const placement_v& rx_locations,
        double_v& raw_cow_data,
        const unsigned& grid_rows,
        const unsigned& grid_cols,
        const double& min_color_span,
        const double& max_color_span)
    {
        float pixel_height = 1;
        float pixel_width = 1;

        sf::RenderTexture renderTexture;

        if (renderTexture.create(grid_cols * pixel_width, grid_rows * pixel_height))
        {
            renderTexture.clear();

            HeatGrid griddata(grid_cols, grid_rows, min_color_span, max_color_span, renderTexture.getSize(), rx_locations, tx_states);

            /* draw the grid */
            griddata.update_heat(raw_cow_data);
            renderTexture.draw(griddata.grid);

            /* draw the receivers */
            for (auto& circle : griddata.rx_cicles)
            {
                renderTexture.draw(circle);
            }

            /* draw the transmitters */
            for (auto& data : griddata.txdata)
            {
                renderTexture.draw(data.transmitter);
                for (auto& element : data.indicators)
                {
                    renderTexture.draw(element);
                }
            }

            renderTexture.display();

            // Save the render texture to a file
            sf::Texture texture = renderTexture.getTexture();
            sf::Image image = texture.copyToImage();
            if (!image.saveToFile(filename))
            {
                spdlog::error("Failed to save image!");
            }

            return;
        }

        spdlog::error("Failed to create render texture!");
    }
}
