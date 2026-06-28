cmake_minimum_required(VERSION 3.14)
project(PhaseArray VERSION 1.0 LANGUAGES CXX)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set optimization flags for different build types
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0 -g)
	message(STATUS "Debug Build")
else() #(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Release build flags
    add_compile_options(-O3)
	message(STATUS "Release Build")
endif()

# FetchContent for SFML, ImGui, and argparse
include(FetchContent)

# SFML
FetchContent_Declare(
    SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 3.0.0
)

# set(SFML_DIR ${sfml_SOURCE_DIR}/build)
# file(MAKE_DIRECTORY ${SFML_DIR}) # if it doesnt't exist

FetchContent_MakeAvailable(SFML)

# ImGui and ImGui-SFML integration
FetchContent_Declare(
    ImGui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.91.7
)
FetchContent_MakeAvailable(ImGui)

FetchContent_Declare(
    ImGui-SFML
    GIT_REPOSITORY https://github.com/eliasdaler/imgui-sfml.git
    GIT_TAG v3.0
)
FetchContent_MakeAvailable(ImGui-SFML)

# argparse
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/morrisfranken/argparse.git
)
FetchContent_MakeAvailable(argparse)

# Add source files
set(SOURCES
    src/main.cpp
    src/antenna.cpp
    src/visuals.cpp
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME}
    PRIVATE
    sfml-system sfml-window sfml-graphics
    ImGui-SFML::ImGui-SFML
    argparse
)

# Include directories for header files
target_include_directories(${PROJECT_NAME}
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/source/inc
)

# OpenGL support for future use
# find_package(OpenGL)
# target_link_libraries(${PROJECT_NAME} PRIVATE OpenGL::GL)

# Optionally, set compiler warnings
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic)
endif()

execute_process(COMMAND date OUTPUT_VARIABLE time)

message(STATUS "-----------------------------------------------")
message(STATUS "\n${CMAKE_BUILD_TYPE} Build Done ${time}")
