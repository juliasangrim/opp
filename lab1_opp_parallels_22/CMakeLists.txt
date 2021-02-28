cmake_minimum_required(VERSION 3.17)
project(lab1_opp_parallels_22 C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
add_executable(lab1_opp_parallels_22 main.c)
target_link_libraries(lab1_opp_parallels_22 m)