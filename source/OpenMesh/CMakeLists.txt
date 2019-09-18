set (CMAKE_CXX_COMPILER "/usr/bin/g++")
set (CMAKE_C_COMPILER "/usr/bin/gcc")
aux_source_directory(./core OPENMESH_DIR)
add_library (OpenMesh SHARED ${OPENMESH_DIR})