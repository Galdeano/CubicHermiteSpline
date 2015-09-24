
# Requires at least CMake 2.8 to configure the package.
cmake_minimum_required (VERSION 2.8) 
project (CHSpline) 

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

SET(PROJECT_DESCRIPTION "A library to compute cubic Hermite spline")
SET(PROJECT_URL "http://github.com/Galdeano/CubicHermiteSpline")

INCLUDE(CTest)

add_subdirectory (src) 
add_subdirectory (test) 
