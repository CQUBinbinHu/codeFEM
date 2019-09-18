project(Hello)
add_executable(hello.o helloworld.cpp)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/lib)