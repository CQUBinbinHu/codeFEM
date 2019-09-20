cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_COMPILER "/usr/bin/g++")
set (CMAKE_C_COMPILER "/usr/bin/gcc")
#project info
project(Hello)
#set version
set (OPTIMIZATION_VERSION_MAJOR 0)
set (OPTIMIZATION_VERSION_MINOR 0)

#set (CMAKE_BUILD_TYPE "Release")
set (CMAKE_BUILD_TYPE "Debug")
# add sub libs
add_subdirectory(./source/NumericalMesh)
add_subdirectory(./source/NumericalMethods)
#add_subdirectory(./source/OpenMesh)
# Target
aux_source_directory(./source/Test/Optimization/ TEST_PROJ)
add_executable(optimization.o ${TEST_PROJ})
# link
target_link_libraries(optimization.o 
    NumericalMesh
    NumericalMethods)
# set output directories
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/lib)