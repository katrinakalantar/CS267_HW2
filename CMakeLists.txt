cmake_minimum_required(VERSION 3.3)
project(CS267_HW2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    autograder.cpp
    common.cpp
    common.h
    Geohash.cpp
    GeoRegion.cpp
    GeoRegion.h
    Makefile
    mpi.cpp
    openmp.cpp
    pthreads.cpp
    quadtree.cpp
    quadtree.h
    Region.cpp
    Region.h
    serial.cpp)

add_executable(CS267_HW2 ${SOURCE_FILES})