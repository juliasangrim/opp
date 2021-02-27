cmake_minimum_required(VERSION 3.17)
project(lab1_opp_parallels_1 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
add_executable(lab1_opp_parallels_1 main.c)
target_link_libraries(lab1_opp_parallels_1 m)