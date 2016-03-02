cmake_minimum_required(VERSION 3.3)
project(CS267_HW2_geohash)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    common.cpp
    common.h
    GeoRegion.cpp
    GeoRegion.h
    Region.cpp
    Region.h
    serial.cpp)

add_executable(CS267_HW2_geohash ${SOURCE_FILES})