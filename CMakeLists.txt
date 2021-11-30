cmake_minimum_required(VERSION 3.20)
project(test)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(test
        main.cpp
        matrix.cpp
        matrix.h)
