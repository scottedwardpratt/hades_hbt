cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})
set(GSL_INCLUDE_DIR $ENV{GSL_INCLUDE_DIR})
project(hades_hbt_project CXX C)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
include_directories(SYSTEM include ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}")

file(GLOB hades_hbt_src ./src/*.cc)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(hades_hbt STATIC ${hades_hbt_src})
target_link_libraries(hades_hbt ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES})



