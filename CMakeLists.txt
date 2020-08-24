cmake_minimum_required(VERSION 3.15)
project(extCatalan)

set(CMAKE_CXX_STANDARD 14)

add_executable(catalan
        staircase.cpp rings.h rings.cpp combinatorics.h)