cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(bicycles CXX C)

enable_testing()

#search for libcloudphxx
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

# search for libmpdata
find_package(libmpdata++)
if (NOT libmpdataxx_FOUND)
  message(FATAL_ERROR "libmpdata++ not found")
endif()

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_subdirectory(drawbicyc)
add_subdirectory(src) 
add_subdirectory(tests)
