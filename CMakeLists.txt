cmake_minimum_required(VERSION 3.5)
project(dparp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(/opt/gurobi912/linux64/include /opt/boost_1_76_0 /opt/lemon-1.3.1/build) # /opt/lemon-1.3.1)

link_directories(/opt/gurobi912/linux64/lib /opt/boost_1_76_0 /opt/lemon-1.3.1/build/lemon)

set(SOURCE_FILES main.cpp headers/Graph.h src/Graph.cpp headers/Include.h src/Arc.cpp headers/Arc.h headers/Model.h src/Model.cpp)

add_executable(dparp ${SOURCE_FILES})

target_link_libraries(dparp gurobi_c++ gurobi91 pthread m emon)
