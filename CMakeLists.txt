cmake_minimum_required(VERSION 3.10)
#project info
project(Hello)
#set version
set (FEM_VERSION_MAJOR 0)
set (FEM_VERSION_MINOR 0)
#set (CMAKE_BUILD_TYPE "Release")
set (CMAKE_BUILD_TYPE "Debug")
# add source files
aux_source_directory(./source/ TEST_PROJ)
# Target
add_executable(hello.o ${TEST_PROJ})
# set output directories
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/lib)