set (CMAKE_CXX_COMPILER "/usr/bin/g++")
set (CMAKE_C_COMPILER "/usr/bin/gcc")
aux_source_directory(./core NUMERICALMESH_DIR)
add_library (NumericalMesh SHARED ${NUMERICALMESH_DIR})