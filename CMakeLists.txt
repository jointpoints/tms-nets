cmake_minimum_required(VERSION 3.8)
project(tms-nets VERSION 1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES *.cpp *.hpp)		

include_directories(include)

add_subdirectory(tests)
