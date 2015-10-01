
# Requires at least CMake 2.8 to configure the package.
cmake_minimum_required (VERSION 2.8) 

include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/pkg-config.cmake)

set(PROJECT_NAME CubicHermiteSpline)
set(PROJECT_VERSION 1.3)
set(PROJECT_DESCRIPTION "A library to compute cubic Hermite spline")
set(PROJECT_URL "http://github.com/Galdeano/CubicHermiteSpline")

set(HEADERS
  ${PROJECT_SOURCE_DIR}/include/CHSpline/CHSpline.h
)

option(VALGRIND_CHECK "VALGRIND_CHECK" OFF)
if (VALGRIND_CHECK)
    FIND_PROGRAM( VALGRIND_PATH valgrind )
    if (NOT VALGRIND_PATH)
        MESSAGE(FATAL_ERROR "valgrind not found! Aborting...")
    endif () # NOT VALGRIND_PATH
    set(VALGRIND_OPTIONS "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
    MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_MEMCHECK
    CMAKE_C_FLAGS_MEMCHECK
    CMAKE_EXE_LINKER_FLAGS_MEMCHECK
    CMAKE_SHARED_LINKER_FLAGS_MEMCHECK )
endif (VALGRIND_CHECK)

option(ERROR_CHECK "Activate all errors" ON)
if (ERROR_CHECK)
    if (CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
    endif (CMAKE_COMPILER_IS_GNUCXX)
endif (ERROR_CHECK)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(CMAKE_C_FLAGS="${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
set(CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")

INCLUDE(CTest)

setup_project()

add_subdirectory (src) 
add_subdirectory (test) 

if (VALGRIND_CHECK)
    ADD_CUSTOM_TARGET(valgrind_check
    
    COMMAND ${VALGRIND_PATH} --trace-children=yes --quiet --tool=memcheck
                             --leak-check=yes --show-reachable=yes
                             --num-callers=50 --verbose --demangle=yes
                             CHSplineDemo
    COMMENT "Valgrind run complete for CHSplineDemo"
    )
endif (VALGRIND_CHECK)

header_install(${HEADERS})
setup_project_finalize()
