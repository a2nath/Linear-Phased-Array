#pragma once
#include <iostream>
#include <csignal>
#include "common.h"
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

    struct HeatGrid
    {
        //const std::vector<vector<T>>& raw_values;
        //sf::VertexArray vertices;
        size_t data_rows, data_cols;
        unsigned pixel_height, pixel_width;

        double min_pxl, max_pxl;
        float thresholds[3];

        sf::Vector2u bounds_lower, bounds_upper, window_size;
        sf::VertexArray vertices;

        sf::Color color(const double& raw)
        {
            double color_value = (max_pxl == min_pxl) ? 0.0 : ((raw - min_pxl) / (max_pxl - min_pxl));

            // Interpolated Color : Transition from red -> yellow -> green -> blue
            if (color_value <= thresholds[0]) {
                return sf::Color(0, static_cast<sf::Uint8>(color_value * 4 * 255), 255); // Blue to Cyan
            }
            else if (color_value <= thresholds[1]) {
                return sf::Color(0, 255, static_cast<sf::Uint8>((1 - (color_value - thresholds[0]) * 4) * 255)); // Cyan to Green
            }
            else if (color_value <= thresholds[2]) {
                return sf::Color(static_cast<sf::Uint8>((color_value - thresholds[1]) * 4 * 255), 255, 0); // Green to Yellow
            }
            else {
                return sf::Color(255, static_cast<sf::Uint8>((1 - (color_value - thresholds[2]) * 4) * 255), 0); // Yellow to Red
            }
        }

        void update_panning(sf::Vector2i& moved_offset)
        {
            bounds_lower = { (unsigned)max((float)moved_offset.x, (float)0.0), (unsigned)max((float)moved_offset.y, (float)0.0) };
            bounds_upper = { (unsigned)min((float)window_size.x, (float)window_size.x + moved_offset.x), (unsigned)min((float)window_size.y, (float)window_size.y + moved_offset.y) };
            moved_offset = { 0, 0 };
        }

        void update_heat(const double_v& tx_raw_sigdata)
        {
            size_t index = 0;
            for (size_t row_idx = data_rows - 1; row_idx >= 0; --row_idx)
            {
                auto row_offset = row_idx * data_cols;
                for (size_t cols_idx = 0; cols_idx < data_cols; ++cols_idx)
                {
                    size_t v_index = (row_offset + cols_idx) << 2;

                    sf::Color v_color = color(tx_raw_sigdata[index++]);
                    vertices[v_index + 0].color = v_color;
                    vertices[v_index + 1].color = v_color;
                    vertices[v_index + 2].color = v_color;
                    vertices[v_index + 3].color = v_color;
                }

                if (row_idx == 0)
                    break;
            }
        }

        void init(const sf::Vector2u& ibounds_lower, const sf::Vector2u& ibounds_upper)
        {
            bounds_lower = ibounds_lower;
            bounds_upper = ibounds_upper;

            for (size_t row_idx = data_rows - 1; row_idx >= 0; --row_idx)
            {
                auto row_offset = row_idx * data_cols;
                for (size_t cols_idx = 0; cols_idx < data_cols; ++cols_idx)
                {
                    size_t v_index = (row_offset + cols_idx) << 2;

                    vertices[v_index + 0].position = sf::Vector2f(cols_idx * pixel_width, row_idx * pixel_height);
                    vertices[v_index + 1].position = sf::Vector2f((cols_idx + 1) * pixel_width, row_idx * pixel_height);
                    vertices[v_index + 2].position = sf::Vector2f((cols_idx + 1) * pixel_width, (row_idx + 1) * pixel_height);
                    vertices[v_index + 3].position = sf::Vector2f(cols_idx * pixel_width, (row_idx + 1) * pixel_height);
                }

                if (row_idx == 0)
                    break;
            }
        }

        HeatGrid(const size_t& irows, const size_t& icols, const double& imin, const double& imax, const sf::Vector2u& iwindow_size)
            :
            //raw_values(values),
            data_rows(irows),
            data_cols(icols),
            pixel_height(iwindow_size.y / irows),
            pixel_width(iwindow_size.x / icols),
            min_pxl(imin),
            max_pxl(imax),
            bounds_lower({0, 0}),
            bounds_upper(iwindow_size),
            window_size(iwindow_size),
            vertices(sf::Quads, irows * icols * 4)
            //vertices(sf::Quads, values.size() * 4)
        {
            thresholds[0] = 0.40; // Cyan to Green
            thresholds[1] = 0.55; // Green to Yellow
            thresholds[2] = 0.71; // Yellow to Red

            init(bounds_lower, bounds_upper);
// last good
            //thresholds[0] = 0.40; // Cyan to Green
            //thresholds[1] = 0.55; // Green to Yellow
            //thresholds[2] = 0.71; // Yellow to Red
            //thresholds[0] = 0.50; // Cyan to Green
            //thresholds[1] = 0.60; // Green to Yellow
            //thresholds[2] = 0.70; // Yellow to Red
            //thresholds[0] = 0.40; // Cyan to Green
            //thresholds[1] = 0.50; // Green to Yellow
            //thresholds[2] = 0.70; // Yellow to Red
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

    /* look for a non-inf number */
    inline void validate_ite(const double_v& raw_values, double_v::iterator& iminmax)
    {
        while (std::isinf(*iminmax))
        {
            if (iminmax != raw_values.begin())
            {
                iminmax -= 1;
            }
            else if (iminmax != raw_values.end() - 1)
            {
                iminmax += 1;
            }
        }
    }

   // inline void change_zoom(zoom()

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
        const std::vector<Placements>& rx_locations,
        const std::vector<Placements>& tx_locations,
        std::vector<double_v>& raw_values,
        const double_v& ant_direction,
        const double_v& scan_angle,
        const size_t& rows,
        const size_t& cols,
        const double& min_ptx,
        const double& max_ptx)
    {
        std::signal(SIGINT, sig_handler);


        sf::RenderWindow window(sf::VideoMode(rows, cols), "SFML Grid Plot");
        window.setVerticalSyncEnabled(true);

        auto window_size = window.getSize();


        sf::Vector2f curr_position;
        sf::Vector2i panning_view;
        sf::Vector2i moved_offset;

        sf::Clock delta_clock;
        sf::View view = window.getDefaultView();

        bool state_changed = true;
        bool panning = false;
        float zoomLevel = 1.0f;
        float mouse_delta_thresh = 0.01f;
        float zoom_change_factor = 1.1f;
        long pan_adj_factor = 10;
        int render_cow_id = 0;
        unsigned tx_count = tx_locations.size();

        /* heat data contains vertices too */
        HeatGrid griddata(rows, cols, min_ptx, max_ptx, window_size);

        /* init the heatmap to display heat from TX id */
        griddata.update_heat(raw_values[render_cow_id]);


#ifdef CONTROLS
        ImGui::SFML::Init(window);
#endif


        // Main loop
        while (window.isOpen())
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
                case sf::Event::MouseWheelScrolled:
                {
                    if (event.mouseWheelScroll.delta - mouse_delta_thresh > 0)
                    {
                        zoom_in(window, view, zoomLevel, zoom_change_factor);
                    }
                    else if (event.mouseWheelScroll.delta - mouse_delta_thresh < 0)
                    {
                        zoom_out(window, view, zoomLevel, zoom_change_factor);
                    }
                    break;
                }
                case sf::Event::MouseButtonPressed:
                {
                    if (event.mouseButton.button == sf::Mouse::Right)
                    {
                        panning_view = sf::Mouse::getPosition(window);
                        panning = true;
                    }
                    break;
                }
                case sf::Event::MouseButtonReleased:
                {
                    if (event.mouseButton.button == sf::Mouse::Right)
                    {
                        griddata.update_panning(moved_offset);

                        //auto newpos = sf::Mouse::getPosition(window);
                        //moved_offset = newpos - panning_view;
                        //
                        //if (moved_offset.x - mouse_delta_thresh != 0 || moved_offset.y != 0)
                        //{
                        //    pan_window(window, view, curr_position, moved_offset);
                        //    griddata.update_panning(moved_offset);
                        //}
                        panning = false;
                    }
                    break;
                }
                case sf::Event::MouseMoved:
                {
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
                    case sf::Keyboard::Scan::R:
                    {
                        zoomLevel = 1.0f;  // Reset zoom level
                        view = window.getDefaultView();
                        window.setView(view);
                        break;
                    }
                    case sf::Keyboard::Scan::Tab:
                    {
                        render_cow_id = (render_cow_id + 1) % tx_count;
                        griddata.update_heat(raw_values[render_cow_id]);
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

            ImGui::End();
#endif

            // Render everything

            //if (state_changed)
            //{
                window.clear();

                // Draw the grid
                //size_t index = 0;
                //for (size_t row_idx = bounds_upper.x - 1; row_idx >= bounds_lower.x; --row_idx)
                //{
                //    for (size_t cols_idx = bounds_lower.y; cols_idx < bounds_upper.y; ++cols_idx)
                //    {
                //        /* same as before with RectandleShape */
                //        sf::Color color = griddata.color(cow_data[index]);
                //        int pixel_height = window.getSize().y / rows;
                //        int pixel_width = window.getSize().x / cols;
                //
                //        /* needed if you need to fallback to this method again */
                //        //sf::RectangleShape cell(sf::Vector2f(pixel_height, pixel_width));
                //        //cell.setPosition(cols_idx * pixel_height, row_idx * pixel_width);
                //        //cell.setFillColor(color);
                //
                //        //window.draw(cell);
                //    }
                //
                //    if (row_idx == 0)
                //        break;
                //}

                window.draw(griddata.vertices);

                const auto& bounds_lower = griddata.bounds_lower;
                const auto& bounds_upper = griddata.bounds_upper;


                auto window_height = window.getSize().y;
                for (auto& loc : rx_locations)
                {
                    if (bounds_lower.x <= loc.x && loc.x <= bounds_upper.x && bounds_lower.y <= loc.y && loc.y <= bounds_upper.y)
                    {
                        sf::CircleShape sta(10.0);
                        sta.setFillColor(sf::Color(90, 90, 90));
                        sta.setOutlineColor(sf::Color::Black);
                        sta.setOutlineThickness(2.0f);
                        sta.setPosition(loc.x, window_height - loc.y);
                        window.draw(sta);
                    }
                }

                for (int i =0 ; i < tx_locations.size(); ++i)
                {
                    auto& loc = tx_locations[i];
                    //if (bounds_lower.x <= loc.x && loc.x <= bounds_upper.x && bounds_lower.y <= loc.y && loc.y <= bounds_upper.y)
                    //{
                    //    sf::RectangleShape sta(sf::Vector2f(10.0, 10.0));
                    //    sta.setFillColor(sf::Color(90, 90, 90));
                    //    sta.setOutlineColor(sf::Color::Black);
                    //    sta.setOutlineThickness(1.0f);
                    //    sta.setPosition(loc.x, window_height - loc.y);
                    //    window.draw(sta);
                    //}
                    sf::Vector2f size(10.0f, 10.0f);

                    if (bounds_lower.x <= loc.x && loc.x <= bounds_upper.x && bounds_lower.y <= loc.y && loc.y <= bounds_upper.y)
                    {
                        sf::Vector2f position(loc.x - size.x / 2, window_height - loc.y - size.y / 2);
                        auto dir_radians = M_PIl / 2 - ant_direction[i];
                        auto& angle_radians = scan_angle[i];

                        // Rectangle
                        sf::RectangleShape sta(size);
                        sta.setFillColor(sf::Color(90, 90, 90));
                        sta.setOutlineColor(sf::Color::Black);
                        sta.setOutlineThickness(1.0f);
                        sta.setPosition(position);
                        window.draw(sta);

                        // First line: perpendicular to the direction the tower is facing
                        sf::VertexArray line11(sf::Lines, 2), line12(sf::Lines, 2);


                        // Calculate the endpoint of the first line based on direction
                        float length = 30.0f / 2;  // Halfway across the rectangle

                        // Line 1's position: calculated from the antenna direction
                        float line1StartX = position.x + size.x / 2; // Middle of the rectangle;
                        float line1StartY = position.y + size.y / 2;

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

                        window.draw(line11);
                        window.draw(line12);

                        // Second line: represents the signal beam, based on scan_angle
                        sf::VertexArray line2(sf::Lines, 2);

                        // Calculate the beam angle relative to the first line (add scan_angle to direction)
                        float beam_angle_radians = ant_direction[i] + angle_radians;
                        //cout << beam_angle_radians * (180 / M_PIl) << endl;
                        float beam_length = 25.0f;  // Length of the signal beam line


                          // Line 1's position: calculated from the antenna direction
                        float line2StartX = position.x + size.x / 2; // Middle of the rectangle
                        float line2StartY = position.y + size.y / 2;

                        // Calculate the endpoint of the second line based on beam angle
                        float line2EndX = line2StartX + beam_length * cos(beam_angle_radians);
                        float line2EndY = line2StartY - beam_length * sin(beam_angle_radians);

                        line2[0].position = sf::Vector2f(line2StartX, line2StartY);  // Starts from the end of the first line
                        line2[1].position = sf::Vector2f(line2EndX, line2EndY);

                        line2[0].color = sf::Color::Blue;
                        line2[1].color = sf::Color::Blue;

                        window.draw(line2);
                    }
                }



                //state_changed = false;
            //} // end of state_changed

#ifdef CONTROLS
            ImGui::SFML::Render(window);  // Render ImGui over SFML content
#endif
            window.display();
        }

#ifdef CONTROLS
        ImGui::SFML::Shutdown();
#endif
        return 0;
    }

    /*    G U I    */
    void plot(
        Logger& logger,
        const string& filename,

        double_v& raw_values,
        const size_t& rows,
        const size_t& cols,
        const double& min_ptx,
        const double& max_ptx)
    {
        float pixel_height = 1;
        float pixel_width = 1;

        sf::RenderTexture renderTexture;

        if (renderTexture.create(cols * pixel_width, rows * pixel_height))
        {
            renderTexture.clear();

            HeatGrid griddata(rows, cols, min_ptx, max_ptx, renderTexture.getSize());

            // Draw the grid
            size_t index = 0;
            for (size_t row_idx = rows - 1; row_idx >= 0; --row_idx)
            {
                for (size_t cols_idx = 0; cols_idx < cols; ++cols_idx)
                {
                    sf::Color color = griddata.color(raw_values[index++]);

                    sf::RectangleShape cell(sf::Vector2f(pixel_height, pixel_width));
                    cell.setPosition(cols_idx * pixel_height, row_idx * pixel_width);
                    cell.setFillColor(color);

                    renderTexture.draw(cell);  // Draw the cell onto the renderTextur
                }

                if (row_idx == 0)
                    break;
            }

            renderTexture.display();

            // Save the render texture to a file
            sf::Texture texture = renderTexture.getTexture();
            sf::Image image = texture.copyToImage();
            if (!image.saveToFile(filename))
            {
                printerr(logger, "Failed to save image!");
            }

            return;
        }

        print(logger, "Failed to create render texture!");
    }
}
