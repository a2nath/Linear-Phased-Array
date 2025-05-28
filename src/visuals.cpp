#pragma once

#include "visuals.h"
using namespace rf_math;

/*    G U I    */
int graphics::render(
    Logger& logger,
    const placement_v& init_rx_locations,
    const state_v& init_states,
    state_v& curr_states,
    const std::vector<double_v>& raw_dbg_lin_data,
    std::vector<double_v>& ready_dbg_dBm_data,
    std::vector<double_v>& ready_snr_dB_data,
    const unsigned& grid_rows,
    const unsigned& grid_cols,
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
    int debounce_txid = -1;
    float debounce_timer = 0.0f;
    float debounce_delay = 0.015f; // 15ms delay or 67 FPS
    float zoom_request = 0.0f;
    float zoomLevel = 1.0f;
    bool panning = false;
    float mouse_delta_thresh = 0.01f;
    float zoom_change_factor = 1.1f;
    long pan_adj_factor = 10;
    int render_tx_id = 0;

    State debug_state;

    size_t tx_count = init_states.size();

    /* set the states, current (working variable), previous (for undo), init (for reset) */

    auto& init = init_states;
    auto& curr = curr_states;
    curr = init;

    state_v prev = init;
    txvertex* tx_dragging = nullptr;

    Pair<float> min_and_max = compute_colorspan(ready_snr_dB_data);
    std::vector<std::vector<double>>* ptr_live_data = &ready_snr_dB_data;

    HeatGrid griddata(grid_cols, grid_rows, min_and_max.first, min_and_max.second, window_size, init_rx_locations, curr);

    std::vector<std::string> tx_header, tx_x_slider, tx_x_inp, tx_y_slider, tx_y_inp, tx_power_slider, tx_power_inp, tx_dir_slider, tx_dir_inp, tx_scan_slider, tx_scan_inp;
    std::vector<float> power_dBm(init.size()), theta_deg(init.size()), scan_deg(init.size());

    for (auto i = 0; i < init.size(); ++i)
    {
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
        tx_scan_inp.emplace_back("Angle##input" + sidx);

        power_dBm[i] = watt2dBm(init[i].settings.power);
        theta_deg[i] = rad2deg(init[i].settings.theta_c);
        scan_deg[i] = rad2deg(init[i].settings.alpha);
    }

    sf::Text fpsText;
    fpsText.setFont(griddata.font);
    fpsText.setCharacterSize(16);       // Set a small size
    fpsText.setFillColor(sf::Color::White); // Set text color
    fpsText.setPosition(10.f, 10.f);    // Position at the top-left corner

    /* init the heatmap to display heat from TX id */
    sync.event_render(render_tx_id);

#ifdef CONTROLS
    ImGui::SFML::Init(window);
#endif
    /* only update the heat when INIT or moving MOVING tx on the map */
    std::thread heat_checker([&]()
        {
            while (is_rendering)
            {
                std::unique_lock<std::mutex> lock(graphics::render_mutex);  // Lock the mutex
                graphics::consig.wait(lock, [&]()
                    {
                        return sync.render_tx_id >= 0 || !is_rendering;
                    }
                );

                if (!is_rendering)
                {
                    break;
                }

                if (sync.is_debugging && sync.debug_interrupt)
                {
                    griddata.reset_span(compute_colorspan(*ptr_live_data));
                    sync.debug_interrupt = false;
                }


                griddata.update_heat((*ptr_live_data)[sync.render_tx_id]);
                sync.render_tx_id = -1;
            }
        }
    );

    FPSBench bench;
    bench.bench_start();

    float fps = 0;

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
                sync.event_resize(new_size.x, new_size.y);
                break;
            }
            case sf::Event::MouseWheelScrolled:
            {
                zoom_request = event.mouseWheelScroll.delta - mouse_delta_thresh;
                break;
            }
            case sf::Event::MouseButtonPressed:
            {
                switch (event.mouseButton.button)
                {
                case sf::Mouse::Left:
                {
                    for (auto& tx : griddata.txdata)
                    {
                        if (tx.transmitter.getGlobalBounds().contains(event.mouseButton.x, event.mouseButton.y))
                        {   // found it

                            render_tx_id = tx.id;

                            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Scan::LControl))
                            {
                                tx_dragging = &tx;
                                auto& id = tx_dragging->id;

                                auto potential_new_loc = griddata.grid_loc_2_state_loc(event.mouseButton.x, event.mouseButton.y);
                                if (potential_new_loc != curr[id].location)
                                {
                                    prev[id].location = curr[id].location;

                                    curr[id].location = potential_new_loc;
                                    griddata.set_tx_position(id, event.mouseButton.x, event.mouseButton.y);

                                    debounce_timer = 0.0f;
                                    debounce_txid = id;
                                }
                            }
                            else
                            {
                                sync.event_render(render_tx_id);
                            }

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
                        auto& id = tx_dragging->id;

                        auto potential_new_loc = griddata.grid_loc_2_state_loc(event.mouseButton.x, event.mouseButton.y);
                        if (potential_new_loc != curr[id].location)
                        {
                            curr[id].location = potential_new_loc;
                            griddata.set_tx_position(id, event.mouseButton.x, event.mouseButton.y);

                            debounce_timer = 0.0f;
                            debounce_txid = id;
                        }
                        else
                        {
                            sync.event_render(render_tx_id);
                        }

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
                    auto& id = tx_dragging->id;

                    curr[id].location = griddata.grid_loc_2_state_loc(event.mouseMove.x, event.mouseMove.y);
                    griddata.set_tx_position(id, event.mouseMove.x, event.mouseMove.y);

                    debounce_timer = 0.0f;
                    debounce_txid = id;
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
                case sf::Keyboard::Scan::S:
                {
                    if (event.key.control && event.key.shift)
                    {
                        save(griddata, *ptr_live_data, render_tx_id, &window);
                    }
                    else if (event.key.control)
                    {
                        save(griddata, *ptr_live_data, render_tx_id);
                    }
                    break;
                }
                case sf::Keyboard::Scan::Z:
                {
                    if (event.key.control)
                    {
                        for (int i = 0; i < tx_count; ++i)
                        {
                            if (prev[i] != curr[i])
                            {
                                griddata.set_state_2_grid_loc(curr[i]);

                                sync.emplace_state(curr[i]);
                            }

                            curr[i] = prev[i];
                        }

                        sync.event_render(render_tx_id); // moved the if statement to check for mainq->empty() and made it MT safe

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
                            ptr_live_data = &ready_dbg_dBm_data;
                            griddata.debug_mode = true;

                            if (ptr_live_data->empty())
                            {
                                ptr_live_data->assign(raw_dbg_lin_data.size(), double_v(griddata.data_width * griddata.data_height));

                                debug_state = curr[render_tx_id];
                                sync.event_debug(debug_state);
                            }
                            else
                            {
                                if (debug_state != curr[render_tx_id])
                                {
                                    debug_state = curr[render_tx_id];
                                    sync.event_debug(debug_state);
                                }
                                else
                                {
                                    sync.event_render(render_tx_id);
                                }

                            }

                            sync.debug_interrupt = true;
                        }
                        else
                        {
                            ptr_live_data = &ready_snr_dB_data;
                            griddata.debug_mode = false;

                            griddata.reset_span(min_and_max);
                            sync.event_render(render_tx_id);
                        }

                        sync.set_debug(griddata.debug_mode);
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

                                griddata.set_state_2_grid_loc(curr[i]);
                                griddata.rotation_update(M_PIl / 2 - curr[i].settings.theta_c, i);

                                sync.emplace_state(curr[i]);
                            }
                        }
                    }
                    else
                    {
                        sync.event_render(render_tx_id);
                    }

                    griddata.reset();

                    window.setView(view);
                    break;
                }
                case sf::Keyboard::Scan::Tab:
                {
                    render_tx_id = (render_tx_id + 1) % tx_count;
                    sync.event_render(render_tx_id);
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
                    ptr_live_data = &ready_dbg_dBm_data;
                    griddata.debug_mode = true;

                    if (ptr_live_data->empty())
                    {
                        ptr_live_data->assign(raw_dbg_lin_data.size(), double_v(griddata.data_width * griddata.data_height));

                        debug_state = curr[render_tx_id];
                        sync.event_debug(debug_state);
                    }
                    else
                    {
                        if (debug_state != curr[render_tx_id])
                        {
                            debug_state = curr[render_tx_id];
                            sync.event_debug(debug_state);
                        }
                        else
                        {
                            sync.event_render(render_tx_id);
                        }

                    }

                    sync.debug_interrupt = true;
                }
                else
                {
                    ptr_live_data = &ready_snr_dB_data;
                    griddata.debug_mode = false;

                    griddata.reset_span(min_and_max);
                    sync.event_render(render_tx_id);
                }

                sync.set_debug(griddata.debug_mode);
            }
        }

        if (ImGui::CollapsingHeader("Heatmap Control", ImGuiTreeNodeFlags_DefaultOpen))
        {
            /* Signal Threshold Sliders */
            ImGui::Text("Signal Threshold Control");

            if (ImGui::SliderFloat("Low##slider", &griddata.curr_thresholds[0], 0.0f, 1.0f, "%.2f") ||
                ImGui::InputFloat("Low##input", &griddata.curr_thresholds[0], 0.1f, 1.0f, "%.2f"))
                sync.event_render(render_tx_id);

            if (ImGui::SliderFloat("Mid##slider", &griddata.curr_thresholds[1], 0.0f, 1.0f, "%.2f") ||
                ImGui::InputFloat("Mid##input", &griddata.curr_thresholds[1], 0.1f, 1.0f, "%.2f"))
                sync.event_render(render_tx_id);

            if (ImGui::SliderFloat("High##slider", &griddata.curr_thresholds[2], 0.0f, 1.0f, "%.2f") ||
                ImGui::InputFloat("High##input", &griddata.curr_thresholds[2], 0.1f, 1.0f, "%.2f"))
                sync.event_render(render_tx_id);


            // Ensure thresholds are in the correct order
            if (griddata.curr_thresholds[0] > griddata.curr_thresholds[1])
                std::swap(griddata.curr_thresholds[0], griddata.curr_thresholds[1]);

            if (griddata.curr_thresholds[1] > griddata.curr_thresholds[2])
                std::swap(griddata.curr_thresholds[1], griddata.curr_thresholds[2]);

            if (griddata.curr_thresholds[0] > griddata.curr_thresholds[1])
                std::swap(griddata.curr_thresholds[0], griddata.curr_thresholds[1]);


            /* Add sliders for minand max values */
            ImGui::Text("Signal Min and Max");
            if (ImGui::SliderFloat("Min##slider", &griddata.curr_pxl_range[0], min_and_max.first, min_and_max.second, "%.2f dB") ||
                ImGui::InputFloat("Min##input", &griddata.curr_pxl_range[0], min_and_max.first, min_and_max.second, "%.2f"))
                sync.event_render(render_tx_id);

            if (ImGui::SliderFloat("Max##slider", &griddata.curr_pxl_range[1], min_and_max.first, min_and_max.second, "%.2f dB") ||
                ImGui::InputFloat("Max##input", &griddata.curr_pxl_range[1], min_and_max.first, min_and_max.second, "%.2f"))
                sync.event_render(render_tx_id);

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

                    if (ImGui::SliderInt(tx_x_slider[i].c_str(), &position.x, 0, grid_cols - 1, "%d meters")) // unsigned long -> int?
                    {
                        curr[i].location.x = griddata.grid_2_state_x(i);
                        griddata.update_tx_vertex(i);

                        debounce_timer = 0.0f;
                        debounce_txid = i;
                    }

                    ImGui::SameLine();
                    if (ImGui::InputInt(tx_x_inp[i].c_str(), &position.x, 0, grid_cols - 1, ImGuiInputTextFlags_EnterReturnsTrue))
                    {
                        curr[i].location.x = griddata.grid_2_state_x(i);
                        griddata.update_tx_vertex(i);

                        sync.emplace_state(curr[i]);
                    }

                    if (ImGui::SliderInt(tx_y_slider[i].c_str(), &position.y, 0, grid_rows - 1, "%d meters"))
                    {

                        curr[i].location.y = griddata.grid_2_state_y(i);
                        griddata.update_tx_vertex(i);

                        debounce_timer = 0.0f;
                        debounce_txid = i;
                    }

                    ImGui::SameLine();
                    if (ImGui::InputInt(tx_y_inp[i].c_str(), &position.y, 0, grid_rows - 1, ImGuiInputTextFlags_EnterReturnsTrue))
                    {
                        curr[i].location.y = griddata.grid_2_state_y(i);
                        griddata.update_tx_vertex(i);

                        sync.emplace_state(curr[i]);
                    }

                    /* antenna-power mechanism for each transmitter */
                    ImGui::Text("Antenna Power");

                    if (ImGui::SliderFloat(tx_power_slider[i].c_str(), &power_dBm[i], -30.0f, +30.0f, "%.2f dBm"))
                    {
                        curr[i].settings.power = dBm2watt(power_dBm[i]);

                        debounce_timer = 0.0f;
                        debounce_txid = i;
                    }

                    ImGui::SameLine();
                    if (ImGui::InputFloat(tx_power_inp[i].c_str(), &power_dBm[i], -30.0f, +30.0f, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue))
                    {
                        curr[i].settings.power = dBm2watt(power_dBm[i]);

                        sync.emplace_state(curr[i]);
                    }

                    ImGui::Text("Antenna Direction");

                    if (ImGui::SliderFloat(tx_dir_slider[i].c_str(), &theta_deg[i], 0.0f, 359.9f, "%.2f deg"))
                    {
                        curr[i].settings.theta_c = deg2rad(theta_deg[i]);
                        griddata.rotation_update(M_PIl / 2 - curr[i].settings.theta_c, i);

                        debounce_timer = 0.0f;
                        debounce_txid = i;

                    }

                    ImGui::SameLine();
                    if (ImGui::InputFloat(tx_dir_inp[i].c_str(), &theta_deg[i], 0.0f, 359.9f, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue))
                    {
                        curr[i].settings.theta_c = deg2rad(theta_deg[i]);
                        griddata.rotation_update(M_PIl / 2 - curr[i].settings.theta_c, i);

                        sync.emplace_state(curr[i]);
                    }


                    /* antenna-scan angle for each transmitter */
                    ImGui::Text("Scan Angle");

                    if (ImGui::SliderFloat(tx_scan_slider[i].c_str(), &scan_deg[i], -90.0, +90.0, "%.2f deg"))
                    {
                        curr[i].settings.alpha = deg2rad(scan_deg[i]);
                        griddata.scan_angle_update(i);

                        debounce_timer = 0.0f;
                        debounce_txid = i;
                    }

                    ImGui::SameLine();
                    if (ImGui::InputFloat(tx_scan_inp[i].c_str(), &scan_deg[i], -90.0, +90.0, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue))
                    {
                        curr[i].settings.alpha = deg2rad(scan_deg[i]);
                        griddata.scan_angle_update(i);
                        sync.emplace_state(curr[i]);
                    }


                } // end TX header (if)
            } // end TX header section (for loop)
        }

        debounce_timer += ImGui::GetIO().DeltaTime;
        if (debounce_txid != -1 && debounce_timer >= debounce_delay)
        {
            sync.emplace_state(curr[debounce_txid]);

            debounce_timer = 0.0f; // Reset timer after update
            debounce_txid = -1;
        }

        if (zoom_request != 0 && !ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow)) // menu scrolls, grid zooms.
        {
            if (zoom_request > 0)
                zoom_in(window, view, zoomLevel, zoom_change_factor);
            else
                zoom_out(window, view, zoomLevel, zoom_change_factor);
        }

        zoom_request = 0; // reset and forget if tried to "scroll" inside the grid

        if (sync.got_updates())
        {
            consig.notify_one(); // either have [render_tx_id] set or [is_computing] set, not both
        }

        ImGui::End();
#endif
        window.clear();
        griddata.draw(window);

        /* show fps on the top left of the screen */
        bench.mark();
        fps = 1.0f / bench.get();
        fpsText.setString("FPS: " + str(static_cast <int>(fps)));

        window.draw(fpsText);

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