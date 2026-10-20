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
option(PHASEARRAY_BUILD_CUDA "Build the CudaSolution GPU executable" OFF)

include(FetchContent)
find_package(Threads REQUIRED)

if(PHASEARRAY_BUILD_CUDA)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(phasearray_default_cuda_architectures ON)
    endif()
    include(CheckLanguage)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR "PHASEARRAY_BUILD_CUDA=ON requires a CUDA 12.8/12.9 compiler. Install CUDA 12.9 and set CMAKE_CUDA_COMPILER to nvcc if needed.")
    endif()
    enable_language(CUDA)
    find_package(CUDAToolkit 12.8 REQUIRED)
    if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL 13)
        message(FATAL_ERROR "CudaSolution requires CUDA 12.8 or 12.9 to retain GTX 1080 support; install CUDA 12.9.")
    endif()
    if(phasearray_default_cuda_architectures)
        set(CMAKE_CUDA_ARCHITECTURES "61;75;86;89;120" CACHE STRING "CudaSolution GPU architectures" FORCE)
    endif()
endif()

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
add_library(phasearray_dependencies INTERFACE)
target_include_directories(phasearray_dependencies INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/inc"
    "${CMAKE_CURRENT_SOURCE_DIR}/deps"
    "${CMAKE_CURRENT_SOURCE_DIR}/deps/include"
    "${argparse_SOURCE_DIR}/include/argparse"
    "${rapidjson_SOURCE_DIR}/include"
)
target_link_libraries(phasearray_dependencies INTERFACE
    sfml-graphics
    sfml-window
    sfml-system
    ImGui-SFML::ImGui-SFML
    argparse
    spdlog::spdlog
    jsoncpp_static
    Threads::Threads
)
target_link_libraries(PhaseArray PRIVATE phasearray_dependencies)
target_compile_definitions(PhaseArray PRIVATE "$<$<CONFIG:Debug>:_DEBUG>")

if(PHASEARRAY_BUILD_CUDA)
    add_executable(CudaSolution src/main.cu src/compute.cu src/visuals.cpp)
    set_target_properties(CudaSolution PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_EXTENSIONS OFF
        CUDA_RUNTIME_LIBRARY Static
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    target_link_libraries(CudaSolution PRIVATE phasearray_dependencies CUDA::cudart_static)
    target_compile_definitions(CudaSolution PRIVATE "$<$<CONFIG:Debug>:_DEBUG>")
    target_compile_options(CudaSolution PRIVATE
        "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Debug>>:-G>"
    )
    if(MSVC)
        target_compile_options(CudaSolution PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:/W4>")
    else()
        target_compile_options(CudaSolution PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>")
    endif()
endif()

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
