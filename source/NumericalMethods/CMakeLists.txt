set (CMAKE_CXX_COMPILER "/usr/bin/g++")
set (CMAKE_C_COMPILER "/usr/bin/gcc")
aux_source_directory(./core NUMERICALMETHODS_DIR)
add_library (NumericalMethods SHARED ${NUMERICALMETHODS_DIR})