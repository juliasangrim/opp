cmake_minimum_required(VERSION 3.19)
project(lab5)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER "/usr/bin/mpicxx")

add_executable(lab5 main.cpp)
target_link_libraries(lab5 pthread)