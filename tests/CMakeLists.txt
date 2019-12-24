cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)

get_filename_component(CURRDIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" CURRDIR ${CURRDIR})

project(${PROJECT_NAME}_${CURRDIR})

file(GLOB SOURCES *.cpp *.hpp)		

add_executable(${PROJECT_NAME} "${SOURCES}")
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${CURRDIR}")
target_include_directories("${PROJECT_NAME}" PRIVATE ../irrpoly/include)