cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(bicycles CXX C)

enable_testing()

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
  add_custom_target(git_revision.h ALL
   git log -1 "--format=format:#define UWLCM_GIT_REVISION \"%H\"%n" HEAD > git_revision.h
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} VERBATIM
  )
endif()

# handle the mpdata options
if(NOT MPDATA_OPTS)
  set(MPDATA_OPTS "IGA;FCT") # iga, fct is the default
endif()
message(STATUS "UWLCM will use the following libmpdata++ options: ${MPDATA_OPTS}")

add_subdirectory(src) 
add_subdirectory(tests)
