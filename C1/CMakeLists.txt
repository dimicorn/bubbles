cmake_minimum_required(VERSION 3.20)
project(C1)

set(CMAKE_CXX_STANDARD 14)

add_executable(C1 main.cpp)

include_directories(odeint-v2/include)