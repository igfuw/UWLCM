cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0028 NEW)

find_package(HDF5 COMPONENTS CXX HL REQUIRED)

add_subdirectory(unit)

#################################################
# find drawbicyc needed by the moist_thermal test
find_package(drawbicyc)
if (NOT drawbicyc_FOUND)
  message(WARNING "drawbicyc not found, the moist_thermal test will not be built")
else()
  add_subdirectory(moist_thermal)
endif()
#################################################
#add_subdirectory(DYCOMS2D)
