cmake_minimum_required(VERSION 3.7 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/common.cmake)

project("SimpleLBM")

# Define the project's flags:
#set(${PROJECT_NAME}_C_FLAGS "-O3" CACHE STRING "${PROJECT_NAME} C FLAGS")

# Enable the project's flags:
#set(CMAKE_C_FLAGS "-O3" CACHE STRING "cflags")
#set(CMAKE_C_FLAGS "-O3")
message(STATUS "CMAKE_C_FLAGS=${CMAKE_C_FLAGS} is loaded...")
include_directories(${PAPI_ROOT}/include)
#find_library(PAPILIB papi ${PAPI_ROOT}/lib)
find_library(PAPILIB libpapi.a ${PAPI_ROOT}/lib)

# find openmp
find_package(OpenMP)
if (NOT OPENMP_FOUND)
    #message("openmp found")
    message(FATAL_ERROR "openmp not found")
endif(NOT OPENMP_FOUND)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

# Enter the subdirectories:
add_subdirectory(origin)
add_subdirectory(origin_papi)
add_subdirectory(origin_openmp)
add_subdirectory(origin_openmp_papi)

add_subdirectory(combine)
add_subdirectory(combine_papi)
add_subdirectory(combine_openmp)
add_subdirectory(combine_openmp_papi)

add_subdirectory(tight)
add_subdirectory(tight_papi)
add_subdirectory(tight_openmp)
add_subdirectory(tight_openmp_papi)

add_subdirectory(tight_block)
add_subdirectory(tight_block_papi)
add_subdirectory(tight_block_openmp)
add_subdirectory(tight_block_openmp_papi)

add_subdirectory(panel)
add_subdirectory(panel_papi)
add_subdirectory(panel_openmp)
add_subdirectory(panel_openmp_papi)
