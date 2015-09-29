
# Requires at least CMake 2.8 to configure the package.
cmake_minimum_required (VERSION 2.8) 

include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/pkg-config.cmake)


set(PROJECT_NAME CubicHermiteSpline)
set(PROJECT_DESCRIPTION "A library to compute cubic Hermite spline")
set(PROJECT_URL "http://github.com/Galdeano/CubicHermiteSpline")

set(HEADERS
  ${PROJECT_SOURCE_DIR}/include/CHSpline/CHSpline.h
)

IF(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(CMAKE_C_FLAGS="${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
set(CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")

INCLUDE(CTest)

setup_project()

add_subdirectory (src) 
add_subdirectory (test) 

header_install(${HEADERS})
setup_project_finalize()
