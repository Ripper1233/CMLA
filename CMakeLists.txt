cmake_minimum_required(VERSION 3.16)
project(cmla)

set(CMAKE_CXX_STANDARD 17)

set(SOURCES main.cpp)

set(LIBRARIES cmla_1_1.cpp cmla_1_2.cpp cmla_1_3.cpp cmla_2.cpp cmla_3.cpp cmla_4.cpp cmla_5.cpp cmla_6.cpp cmla_7_1.cpp cmla_7_2.cpp)

add_library(algorithms ${LIBRARIES})

add_executable(cmla ${SOURCES})

target_link_libraries(cmla algorithms)
