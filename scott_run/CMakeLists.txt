cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})
project(hades_hbt_project CXX)

add_custom_target(extern_msu_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_msu_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)
add_custom_target(extern_hades_hbt COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/hades_hbt/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wmissing-declarations -std=c++11")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/coral/software/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/hades_hbt/software/include)

link_directories(${GITHOME_MSU}/coral/software/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/hades_hbt/software/lib)

find_library(${GITHOME_MSU}/coral/software/lib msu_coral)
find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
find_library(${GITHOME_MSU}/hades_hbt/software/lib hades_hbt)

set(SOURCE_FILES ./hades_hbt_fromOSCAR.cc)
add_executable(hades_hbt_fromOSCAR ${SOURCE_FILES})
SET_SOURCE_FILES_PROPERTIES(
	hades_hbt_fromOSCAR.cc PROPERTIES OBJECT_DEPENDS
    ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a
 )
SET_SOURCE_FILES_PROPERTIES(
 	hades_hbt_fromOSCAR.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a
)
SET_SOURCE_FILES_PROPERTIES(
 	hades_hbt_fromOSCAR.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/hades_hbt/software/lib/libhades_hbt.a
)
add_dependencies(hades_hbt_fromOSCAR extern_hades_hbt extern_msu_coral extern_msu_commonutils)
target_link_libraries(hades_hbt_fromOSCAR hades_hbt msu_commonutils msu_coral ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a ${GITHOME_MSU}/hades_hbt/software/lib/libhades_hbt.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})


set(SOURCE_FILES ./hades_hbt_fromGAUSS.cc)
add_executable(hades_hbt_fromGAUSS ${SOURCE_FILES})
SET_SOURCE_FILES_PROPERTIES(
	hades_hbt_fromGAUSS.cc PROPERTIES OBJECT_DEPENDS
    ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a
 )
SET_SOURCE_FILES_PROPERTIES(
 	hades_hbt_fromGAUSS.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a
)
SET_SOURCE_FILES_PROPERTIES(
 	hades_hbt_fromGAUSS.cc PROPERTIES OBJECT_DEPENDS ${GITHOME_MSU}/hades_hbt/software/lib/libhades_hbt.a
)
add_dependencies(hades_hbt_fromGAUSS extern_hades_hbt extern_msu_coral extern_msu_commonutils)
target_link_libraries(hades_hbt_fromGAUSS hades_hbt msu_coral msu_commonutils ${GITHOME_MSU}/commonutils/software/lib/libmsu_commonutils.a ${GITHOME_MSU}/coral/software/lib/libmsu_coral.a ${GITHOME_MSU}/hades_hbt/software/lib/libhades_hbt.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})