# This is a CMake example for Python
cmake_minimum_required(VERSION 2.8.8)
project(clipper-wrapper C CXX)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "\${ORIGIN}/..")

set(CMAKE_CXX_FLAGS "-std=c++11")
 
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

if(APPLE)
find_program(PYTHON NAMES python3.5 python3.5m python.exe
                    HINTS /home/tic20/apps/chimerax/bin
                    PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
else()
find_program(PYTHON NAMES python3.5 python3.5m python.exe
                    HINTS /home/tic20/apps/chimerax/bin
                    PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
endif()
message(STATUS "Python interpreter: ${PYTHON}")

#if(APPLE OR WIN32)
  # avoid FindPythonLibs as it can find libpython from different interpreter
  execute_process(COMMAND ${PYTHON} -c "import sys; sys.stdout.write(sys.prefix)"
                  OUTPUT_VARIABLE _PY_PREFIX
                  RESULT_VARIABLE failure)
  if (failure)
    message(FATAL_ERROR "Failed to run ${PYTHON}: ${failure}")
  endif()
  message(STATUS "Python sys.prefix: ${_PY_PREFIX}")
  find_library(PYTHON_LIBRARY NAMES python3.5m python35
               HINTS /home/tic20/apps/chimerax/lib ${_PY_PREFIX}/lib ${_PY_PREFIX}/lib64 ${_PY_PREFIX}/libs
               NO_DEFAULT_PATH)
  if (NOT PYTHON_LIBRARY)
    message(FATAL_ERROR "Python library not found.")
  endif()
  message(STATUS "Python library: ${PYTHON_LIBRARY}")
  find_path(PYTHON_INCLUDE_DIR NAMES Python.h
            HINTS ${_PY_PREFIX}/include
                  ${_PY_PREFIX}/include/python2.7
                  ${_PY_PREFIX}/include/python3.5m
            NO_DEFAULT_PATH)
  message(STATUS "Python.h is in: ${PYTHON_INCLUDE_DIR}")
  if (NOT PYTHON_INCLUDE_DIR)
    message(FATAL_ERROR "Python header not found.")
  endif()
  set(PYTHON_LIBRARIES ${PYTHON_LIBRARY})
#endif()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET(CCP4_DIR /home/tic20/ccp4/devtools/install)

find_package(CCP4 REQUIRED mmdb2 clipper-core clipper-mmdb clipper-ccp4
                           clipper-cif clipper-minimol clipper-contrib ccp4c)

find_package(NumPy REQUIRED)

include_directories(${CCP4_INCLUDE_DIRS}
                    "${CCP4_INCLUDE_DIRS}/clipper"
                    "${CCP4_INCLUDE_DIRS}/clipper/core"
                    "${CMAKE_SOURCE_DIR}"
                    ${PYTHON_NUMPY_INCLUDE_PATH}
                    "${PYTHON_INCLUDE_DIR}")

SET_SOURCE_FILES_PROPERTIES(clipper.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(clipper python clipper.i )

IF(WIN32)
    SET(SWIG_CLIPPER_MODULE_SUFFIX "pyd")
ELSE()
    SET(SWIG_CLIPPER_MODULE_SUFFIX "so")
ENDIF()

string(REGEX REPLACE "clipper-cif" "clipper-phs" CLIPPER-PHS_LIBRARY ${CLIPPER-CIF_LIBRARY})

SWIG_LINK_LIBRARIES(clipper ${PYTHON_LIBRARIES} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CCP4_LIBRARY}
    ${CLIPPER-CONTRIB_LIBRARY} ${CLIPPER-PHS_LIBRARY} ${CLIPPER-CORE_LIBRARIES} ${CLIPPER-MINIMOL_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY} m)

set ( PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/lib/clipper_python_core)
 
message(STATUS "Python module path: ${PYTHON_MODULE_PATH}")

INSTALL ( FILES ${CMAKE_CURRENT_BINARY_DIR}/_clipper.${SWIG_CLIPPER_MODULE_SUFFIX} DESTINATION ${PYTHON_MODULE_PATH} )
INSTALL ( FILES ${CMAKE_CURRENT_BINARY_DIR}/clipper.py DESTINATION ${PYTHON_MODULE_PATH} )
INSTALL ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/clipper_tools DESTINATION ${PYTHON_MODULE_PATH} )
