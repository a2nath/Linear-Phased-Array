cmake_minimum_required(VERSION 3.18)
project(PhaseArray VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(PHASEARRAY_ENABLE_ASAN "Enable AddressSanitizer in Debug builds" ON)

include(FetchContent)
find_package(Threads REQUIRED)

# The application uses the SFML 2 API. Keep the binding and ImGui versions
# compatible with it, and build only the SFML modules used by PhaseArray.
set(SFML_BUILD_AUDIO OFF CACHE BOOL "Build SFML audio" FORCE)
set(SFML_BUILD_NETWORK OFF CACHE BOOL "Build SFML network" FORCE)
set(SFML_BUILD_EXAMPLES OFF CACHE BOOL "Build SFML examples" FORCE)
set(SFML_BUILD_DOC OFF CACHE BOOL "Build SFML documentation" FORCE)
FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 2.6.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(SFML)

FetchContent_Declare(ImGui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.89.9
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ImGui)
set(IMGUI_DIR "${imgui_SOURCE_DIR}" CACHE PATH "ImGui source directory" FORCE)
set(IMGUI_SFML_FIND_SFML OFF CACHE BOOL "Use the fetched SFML targets" FORCE)
set(IMGUI_SFML_IMGUI_DEMO OFF CACHE BOOL "Build the ImGui demo" FORCE)
FetchContent_Declare(ImGui-SFML
    GIT_REPOSITORY https://github.com/SFML/imgui-sfml.git
    GIT_TAG v2.6.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ImGui-SFML)

FetchContent_Declare(argparse
    GIT_REPOSITORY https://github.com/morrisfranken/argparse.git
    GIT_TAG b1479bf9f2f44010ad79efff00bb9bf8ec56dbab
)
FetchContent_MakeAvailable(argparse)

set(SPDLOG_BUILD_EXAMPLE OFF CACHE BOOL "Build spdlog examples" FORCE)
set(SPDLOG_BUILD_TESTS OFF CACHE BOOL "Build spdlog tests" FORCE)
FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.14.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(spdlog)

set(JSONCPP_WITH_TESTS OFF CACHE BOOL "Build JsonCpp tests" FORCE)
set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF CACHE BOOL "Run JsonCpp tests" FORCE)
set(JSONCPP_WITH_PKGCONFIG_SUPPORT OFF CACHE BOOL "Generate JsonCpp pkg-config files" FORCE)
set(JSONCPP_WITH_CMAKE_PACKAGE OFF CACHE BOOL "Generate JsonCpp CMake package files" FORCE)
FetchContent_Declare(jsoncpp
    GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
    GIT_TAG 1.9.5
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(jsoncpp)

# RapidJSON is header-only; skip its legacy CMake project and package registry.
FetchContent_Declare(rapidjson
    GIT_REPOSITORY https://github.com/Tencent/rapidjson.git
    GIT_TAG v1.1.0
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR include
)
FetchContent_MakeAvailable(rapidjson)

add_executable(PhaseArray
    src/main.cpp
    src/antenna.cpp
    src/visuals.cpp
)
set_target_properties(PhaseArray PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
target_include_directories(PhaseArray PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/inc"
    "${CMAKE_CURRENT_SOURCE_DIR}/deps"
    "${CMAKE_CURRENT_SOURCE_DIR}/deps/include"
    "${argparse_SOURCE_DIR}/include/argparse"
    "${rapidjson_SOURCE_DIR}/include"
)
target_link_libraries(PhaseArray PRIVATE
    sfml-graphics
    sfml-window
    sfml-system
    ImGui-SFML::ImGui-SFML
    argparse
    spdlog::spdlog
    jsoncpp_static
    Threads::Threads
)
target_compile_definitions(PhaseArray PRIVATE "$<$<CONFIG:Debug>:_DEBUG>")

if(MSVC)
    target_compile_options(PhaseArray PRIVATE /W4)
else()
    target_compile_options(PhaseArray PRIVATE -Wall -Wextra -Wpedantic
        "$<$<CONFIG:Debug>:-O0>"
        "$<$<CONFIG:Debug>:-g>"
    )
    if(PHASEARRAY_ENABLE_ASAN)
        target_compile_options(PhaseArray PRIVATE
            "$<$<CONFIG:Debug>:-fsanitize=address>"
            "$<$<CONFIG:Debug>:-fno-omit-frame-pointer>"
        )
        target_link_options(PhaseArray PRIVATE "$<$<CONFIG:Debug>:-fsanitize=address>")
    endif()
endif()
