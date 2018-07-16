# - Find CCP4 libraries
# Find one or more of CCP4 libraries: ccp4c, ccp4f, mmdb, ccif,
# clipper-core, clipper-ccp4, clipper-contrib, clipper-minimol, clipper-mmdb,
# or cctbx.
#
# Once done this will define
#  CCP4_INCLUDE_DIRS - all include directories
#  <name>_LIBRARY - library, name is one of CCP4C, CCP4F, MMDB, MMDB2, CCIF,
#                            SSM, CCP4SRS,
#                            CLIPPER-CORE, CLIPPER-CCP4, CLIPPER-CONTRIB,
#                            CLIPPER-MINIMOL, CLIPPER-MMDB, CLIPPER-CIF,
#                            RFFTW2, FFTW2,
#                            HKLFILE, CCTBX
#  CCIF_LIBRARIES - CCIF_LIBRARY with (if needed) regex library
#  CCP4SRS_LIBRARIES - CCP4SRS_LIBRARY with (if needed) mmdb and zlib libraries
#  CLIPPER-CORE_LIBRARIES - CLIPPER-CORE_LIBRARY with (if needed)
#                           fftw2 and thread library
#  CCP4_LIBRARIES - all requested libraries with necessary dependencies
#  <name>_FOUND - TRUE if library and header files are found.
#
#  CCP4_PY2_DIR - when building the CCP4 suite, CCP4_PY2_DIR is set to
#                 preferred installation directory for Python2 modules.
#                 (defaults to ${LIB_INSTALL_DIR}/py2)
#
# You need to name libraries that you will use as components:
# FIND_PACKAGE(CCP4 COMPONENTS mmdb ccp4c)
# or
# FIND_PACKAGE(CCP4 REQUIRED mmdb ccp4c)
#
# Checking for ccp4srs in some cases (static linking) sets also ZLIB_LIBRARY.
#
# Checking for clipper-core sets also variables for FFTW2 library used
# by clipper: FFTW2_LIBRARY, RFFTW2_LIBRARY, FFTW2_INCLUDE_DIRS,
# and if fftw2 is prefixed with 's' FFTW2_DEFINITIONS=-DFFTW2_PREFIX_S.
#
# Sample usage:
#   FIND_PACKAGE(CCP4 REQUIRED mmdb ccp4)
#   INCLUDE_DIRECTORIES(${CCP4_INCLUDE_DIRS})
#   TARGET_LINK_LIBRARIES(<YourTarget> ${CCP4_LIBRARIES})

if(NOT CCP4_FIND_COMPONENTS)
    message(FATAL_ERROR "Please specify CCP4 libraries/components.")
endif()

include(FindPackageHandleStandardArgs)
include(CheckFunctionExists)
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)

set(_ccp4 $ENV{CCP4})
if(_ccp4)
    set(_clib "${_ccp4}/lib")
    set(_cincl "${_ccp4}/include")
endif()

foreach(_component ${CCP4_FIND_COMPONENTS})
    string(TOUPPER ${_component} _upper)
    string(TOLOWER ${_component} _lower)
    set(_lib_var "${_upper}_LIBRARY")
    find_library(${_lib_var} NAMES ${_lower}
                 HINTS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
                 PATHS ${_clib})
    #message("DEBUG: ${_lib_var} (${_lower}) = ${${_upper}_LIBRARY}")

    if (${_upper} STREQUAL "MMDB")
        set(_header "mmdb/mmdb_defs.h")
    elseif (${_upper} STREQUAL "MMDB2")
        set(_header "mmdb2/mmdb_defs.h")
    elseif (${_upper} STREQUAL "CCP4C")
        set(_header "ccp4/ccp4_general.h")
    elseif (${_upper} STREQUAL "CCP4F")
        # fortran library, no headers
        set(_header)
    elseif (${_upper} STREQUAL "CCIF")
        set(_header "ccif/ccif_defines.h")
    elseif (${_upper} STREQUAL "SSM")
        set(_header "ssm/ssm_vxedge.h")
    elseif (${_upper} STREQUAL "CCP4SRS")
        set(_header "ccp4srs/ccp4srs_manager.h")
    elseif (${_upper} STREQUAL "CLIPPER-CORE")
        find_package(Threads)
        set(_header "clipper/clipper.h")
    elseif (${_upper} STREQUAL "CLIPPER-CCP4")
        set(_header "clipper/clipper-ccp4.h")
    elseif (${_upper} STREQUAL "CLIPPER-CONTRIB")
        set(_header "clipper/clipper-contrib.h")
    elseif (${_upper} STREQUAL "CLIPPER-MINIMOL")
        set(_header "clipper/clipper-minimol.h")
    elseif (${_upper} STREQUAL "CLIPPER-MMDB")
        set(_header "clipper/clipper-mmdb.h")
    elseif (${_upper} STREQUAL "CLIPPER-CIF")
        set(_header "clipper/clipper-cif.h")
    elseif (${_upper} STREQUAL "HKLFILE")
        set(_header "hklfile/columngroups.hh")
    elseif (${_upper} STREQUAL "CCTBX")
        set(_header "cctbx/crystal/symmetry.h")
    else()
        message(FATAL_ERROR "'${_component}' is not a valid CCP4 library.")
    endif()

    if (_header)
        set(_incl_var "${_upper}_INCLUDE_DIR")
        find_path(${_incl_var} ${_header}
                  HINTS ${CMAKE_INSTALL_PREFIX}/include
                  HINTS ${CMAKE_INSTALL_PREFIX}/include/cctbx
                  PATHS ${_cincl})
        message(STATUS "Looking for ${_header} - ${${_incl_var}}")
        find_package_handle_standard_args(CCP4 DEFAULT_MSG
                                          ${_lib_var} ${_incl_var})
        set(CCP4_INCLUDE_DIRS ${CCP4_INCLUDE_DIRS} ${${_incl_var}})
    else()
        find_package_handle_standard_args(CCP4 DEFAULT_MSG
                                          ${_lib_var})
    endif()

    # check if libccif needs -lregex
    if (${_upper} STREQUAL "CCIF")
        set(_CCIF_SRC "int main() {zzs_undump(0,0,0);}")
        foreach (_REGEX_LIB "" regex)
            message (STATUS "Trying _REGEX_LIB=${_REGEX_LIB}")
            set(CMAKE_REQUIRED_LIBRARIES ${CCIF_LIBRARY} ${_REGEX_LIB} m)
            set(_VAR _LINKING_WITH_CCIF_${_REGEX_LIB})
            check_c_source_compiles("${_CCIF_SRC}" ${_VAR})
            if (${_VAR})
                set(CCIF_LIBRARIES ${CCIF_LIBRARY} ${_REGEX_LIB} m)
                break()
            endif()
        endforeach()
        if (NOT ${_VAR})
            message(FATAL_ERROR "Linking with libccif failed.")
        endif()
    endif()

    # check if libccp4srs needs mmdb and zlib
    if (${_upper} STREQUAL "CCP4SRS")
        set(_SAVE ${CMAKE_REQUIRED_INCLUDES})
        set(CMAKE_REQUIRED_INCLUDES "${_SAVE};${CCP4SRS_INCLUDE_DIR}")
        set(_SRS_SRC "#include <ccp4srs/memio_.h>\n"
                      "int main() { ccp4srs::MemIO m; m.read(\"f\");}")
        set(CMAKE_REQUIRED_LIBRARIES ${CCP4SRS_LIBRARY})
        set(_VAR _LINKING_WITH_CCP4SRS)
        check_cxx_source_compiles("${_SRS_SRC}" ${_VAR})
        if (${_VAR})
            set(CCP4SRS_LIBRARIES ${CCP4SRS_LIBRARY})
        else()
            find_library(ZLIB_LIBRARY NAMES z zlib zlib1 zdll
                 HINTS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
                 PATHS ${_clib} c:/MinGW/lib)
            if (NOT ZLIB_LIBRARY)
             message(FATAL_ERROR "Linking with ccp4srs fails. zlib not found.")
            endif()
            if (NOT MMDB2_LIBRARY)
                find_library(MMDB2_LIBRARY NAMES mmdb2
                         HINTS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
                         PATHS ${_clib})
            endif()
            set(_ADD_LIBS ${MMDB2_LIBRARY} ${ZLIB_LIBRARY})
            set(CMAKE_REQUIRED_LIBRARIES ${CCP4SRS_LIBRARY} ${_ADD_LIBS})
            set(_VAR _LINKING_WITH_CCP4SRS_MMDB_AND_ZLIB)
            check_cxx_source_compiles("${_SRS_SRC}" ${_VAR})
            if (NOT ${_VAR})
                message(FATAL_ERROR "Linking with ccp4srs failed.")
            endif()
            set(CCP4SRS_LIBRARIES ${CCP4SRS_LIBRARY} ${_ADD_LIBS})
        endif()
        SET(CMAKE_REQUIRED_INCLUDES ${_SAVE})
    endif()

    # check if clipper-core needs rfftw fftw
    if (${_upper} STREQUAL "CLIPPER-CORE")
        set(_fftw_names sfftw fftw)
        if (${${_lib_var}} MATCHES ".a$")
            # fftw should also be static
            set(_fftw_names libsfftw.a libfftw.a ${_fftw_names})
        endif()
        # first search in the install path, then in system and $CCP4
        find_library(FFTW2_LIBRARY NAMES ${_fftw_names}
                     HINTS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
                     NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
        find_library(FFTW2_LIBRARY NAMES ${_fftw_names}
                     PATHS ${_clib})
        message(STATUS "FFTW2 libraries - ${FFTW2_LIBRARY}")
        if (${FFTW2_LIBRARY} MATCHES "sfftw")
            set (_fftw_prefix s)
            set(FFTW2_DEFINITIONS "-DFFTW2_PREFIX_S")
        endif()
        set(_rfftw_names ${_fftw_prefix}rfftw)
        if (${FFTW2_LIBRARY} MATCHES ".a$")
            set(_rfftw_names lib${_fftw_prefix}rfftw.a ${_rfftw_names})
        else()
        endif()
        find_library(RFFTW2_LIBRARY NAMES ${_rfftw_names}
                     HINTS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib
                     PATHS ${_clib})
        message(STATUS "                - ${RFFTW2_LIBRARY}")
        find_path(FFTW2_INCLUDE_DIRS ${_fftw_prefix}rfftw.h
                  HINTS ${CMAKE_INSTALL_PREFIX}/include
                  PATHS ${_cincl})
        message(STATUS "FFTW2 header directory - ${FFTW2_INCLUDE_DIRS}")

        set(_SAVE ${CMAKE_REQUIRED_INCLUDES})
        set(CMAKE_REQUIRED_INCLUDES "${_SAVE};${CLIPPER-CORE_INCLUDE_DIR}")
        set(_CLIP_SRC "#include <clipper/clipper.h>\n"
                      "int main() {clipper::FFTmap().fft_h_to_x();}")
            set(CMAKE_REQUIRED_LIBRARIES ${CLIPPER-CORE_LIBRARY})
            set(_VAR _LINKING_WITH_CLIPPER_CORE)
            check_cxx_source_compiles("${_CLIP_SRC}" ${_VAR})
        if (NOT ${_VAR})
            set(_ADD_LIBS ${CMAKE_THREAD_LIBS_INIT})
            set(CMAKE_REQUIRED_LIBRARIES ${CLIPPER-CORE_LIBRARY} ${_ADD_LIBS})
            set(_VAR _LINKING_WITH_CLIPPER_CORE_AND_THREADS)
            check_cxx_source_compiles("${_CLIP_SRC}" ${_VAR})
        endif()
        if ((NOT ${_VAR}) AND FFTW2_LIBRARY)
            set(_ADD_LIBS ${RFFTW2_LIBRARY} ${FFTW2_LIBRARY}
                          ${CMAKE_THREAD_LIBS_INIT})
            set(CMAKE_REQUIRED_LIBRARIES ${CLIPPER-CORE_LIBRARY} ${_ADD_LIBS})
            set(_VAR _LINKING_WITH_CLIPPER_CORE_AND_FFTW)
            check_cxx_source_compiles("${_CLIP_SRC}" ${_VAR})
        endif()
        if (NOT ${_VAR})
            message(FATAL_ERROR "Linking with clipper-core failed.")
        endif()
        set(CLIPPER-CORE_LIBRARIES ${CLIPPER-CORE_LIBRARY} ${_ADD_LIBS})
        SET(CMAKE_REQUIRED_INCLUDES ${_SAVE})
    endif()
endforeach()


set(CCP4_LIBRARIES ${CCTBX_LIBRARY}
                   ${HKLFILE_LIBRARY}
                   ${CLIPPER-CCP4_LIBRARY}
                   ${CLIPPER-CIF_LIBRARY}
                   ${CLIPPER-CONTRIB_LIBRARY}
                   ${CLIPPER-MINIMOL_LIBRARY}
                   ${CLIPPER-MMDB_LIBRARY}
                   ${CLIPPER-CORE_LIBRARIES}
                   ${CCP4SRS_LIBRARIES}
                   ${SSM_LIBRARY}
                   ${CCIF_LIBRARIES}
                   ${CCP4F_LIBRARY}
                   ${CCP4C_LIBRARY}
                   ${MMDB2_LIBRARY}
                   ${MMDB_LIBRARY})

if(CCP4_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES CCP4_INCLUDE_DIRS)
endif()
message(STATUS "CCP4 include directory: ${CCP4_INCLUDE_DIRS}")

if(CCP4_LIBRARIES)
    # reverse twice, because REMOVE_DUPLICATES leaves the first
    list(REVERSE CCP4_LIBRARIES)
    list(REMOVE_DUPLICATES CCP4_LIBRARIES)
    list(REVERSE CCP4_LIBRARIES)
endif()

if(DEFINED LIB_INSTALL_DIR)
    set(CCP4_PY2_DIR ${LIB_INSTALL_DIR}/py2)
elseif(DEFINED CMAKE_INSTALL_PREFIX)
    set(CCP4_PY2_DIR ${CMAKE_INSTALL_PREFIX}/lib/py2)
endif()

mark_as_advanced(CCP4C_INCLUDE_DIR MMDB_INCLUDE_DIR MMDB2_INCLUDE_DIR
                 CCIF_INCLUDE_DIR SSM_INCLUDE_DIR
                 CCP4SRS_INCLUDE_DIR
                 CLIPPER-CORE_INCLUDE_DIR CLIPPER-CCP4_INCLUDE_DIR
                 CLIPPER-CIF_INCLUDE_DIR
                 CLIPPER-CONTRIB_INCLUDE_DIR CLIPPER-MINIMOL_INCLUDE_DIR
                 CLIPPER-MMDB_INCLUDE_DIR FFTW2_INCLUDE_DIRS
                 CCTBX_INCLUDE_DIR HKLFILE_INCLUDE_DIR CCP4_PY2_DIR)

