# This is a CMake example for Python
cmake_minimum_required(VERSION 2.8.8)
project(clipper-wrapper C CXX)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

find_package(CCP4 REQUIRED mmdb2 clipper-core clipper-mmdb clipper-ccp4
                           clipper-cif clipper-minimol clipper-contrib ccp4c)

include_directories(${CCP4_INCLUDE_DIRS}
                    "${CCP4_INCLUDE_DIRS}/clipper"
                    "${CCP4_INCLUDE_DIRS}/clipper/core"
                    "${CMAKE_SOURCE_DIR}")

SET_SOURCE_FILES_PROPERTIES(clipper.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(clipper python clipper.i )

SWIG_LINK_LIBRARIES(clipper ${PYTHON_LIBRARIES} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CCP4_LIBRARY}
    ${CLIPPER-CONTRIB_LIBRARY} ${CLIPPER-CORE_LIBRARIES} ${CLIPPER-MINIMOL_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY} m)

