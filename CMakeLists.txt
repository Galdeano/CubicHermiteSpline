
# Requires at least CMake 2.8 to configure the package.
cmake_minimum_required (VERSION 2.8) 
project (CHSpline) 

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
SET(CMAKE_C_FLAGS="${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
SET(CMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")

SET(PROJECT_DESCRIPTION "A library to compute cubic Hermite spline")
SET(PROJECT_URL "http://github.com/Galdeano/CubicHermiteSpline")

INCLUDE(CTest)

add_subdirectory (src) 
add_subdirectory (test) 
