cmake_minimum_required(VERSION 3.5)
project(dparp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -std=c++0x ")

include_directories(/opt/gurobi1101/linux64/include /opt/boost_1_84_0 /opt/lemon-1.3.1)

link_directories(/opt/gurobi1101/linux64/lib /opt/boost_1_84_0 /opt/lemon-1.3.1/lemon)

set(
    sources
    main.cpp
    headers/Graph.h
    src/Graph.cpp
    headers/Include.h
    src/Arc.cpp
    headers/Arc.h
    headers/Model.h
    src/Model.cpp
    headers/WarmStart.h
    src/WarmStart.cpp
    headers/Scenario.h
    src/Scenario.cpp
    headers/StochasticModel.h
    src/StochasticModel.cpp
)

add_executable(dparp ${sources})

target_link_libraries(dparp gurobi_c++ gurobi110 pthread m emon)
