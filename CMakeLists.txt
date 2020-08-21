cmake_minimum_required(VERSION 3.15)
project(cpp)

set(CMAKE_CXX_STANDARD 14)

add_executable(cpp
        staircase.cpp rings.h rings.cpp combinatorics.h)