cmake_minimum_required(VERSION 3.14)
project(Wireless_Systems_Design VERSION 1.0)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set optimization flags for different build types
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Debug build flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")  # No optimization, include debug information
	message(STATUS "Debug Build")
else() #(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Release build flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")  # Highest level of optimization
	message(STATUS "Release Build")
endif()

# FetchContent is available from CMake 3.11 onwards
include(FetchContent)

# Fetch and build jsoncpp if not found
FetchContent_Declare(
    jsoncpp
    GIT_REPOSITORY https://github.com/open-source-parsers/jsoncpp.git
    GIT_TAG master
)

FetchContent_GetProperties(jsoncpp)
if(NOT jsoncpp_POPULATED)
    FetchContent_Populate(jsoncpp)
    add_subdirectory(${jsoncpp_SOURCE_DIR} ${jsoncpp_BINARY_DIR})
endif()

# Fetch and make available argparse
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/morrisfranken/argparse.git
)

FetchContent_MakeAvailable(argparse)

# Add your source files here
set(SOURCES
    source/src/project.cpp
    # Add other source files here if any
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME}
    PRIVATE
    jsoncpp_lib
    argparse
)

# Include directories for header files
target_include_directories(${PROJECT_NAME}
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/source/inc
)

# Optionally, set compiler flags
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic)
endif()

execute_process(COMMAND date OUTPUT_VARIABLE time)

message(STATUS "-----------------------------------------------")
message(STATUS "\n${CMAKE_BUILD_TYPE} Build Done ${time}")
