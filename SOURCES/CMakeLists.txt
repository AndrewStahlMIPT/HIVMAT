cmake_minimum_required(VERSION 3.24)
project(HIVMAT)

include(CMakePrintHelpers)
set(CMAKE_CXX_STANDARD 23)


set(CMAKE_CXX_FLAGS_RELEASE -Wall -Wextra -Wpedantic -O3)
set(CMAKE_CXX_FLAGS_DEBUG -g)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMY_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -DMY_RELEASE")

file(GLOB_RECURSE source *.hpp *.h *.cpp)
add_library(HIVMAT INTERFACE
        integrator.hpp)
target_include_directories(HIVMAT INTERFACE "${PROJECT_SOURCE_DIR}")