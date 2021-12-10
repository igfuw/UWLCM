cmake_minimum_required(VERSION 3.0)
project(UWLCM_tests)
cmake_policy(SET CMP0028 NEW)


enable_testing()

find_package(HDF5 COMPONENTS CXX HL REQUIRED)

add_subdirectory(unit)

#################################################
# find UWLCM_plotters needed by the moist_thermal test
find_package(UWLCM_plotters)
if (NOT UWLCM_plotters_FOUND)
  message(WARNING "UWLCM_plotters not found, the moist_thermal test will not be built")
else()
  add_subdirectory(moist_thermal)
endif()
#################################################
#add_subdirectory(DYCOMS2D)
