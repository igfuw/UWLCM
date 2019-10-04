cmake_minimum_required(VERSION 3.0)
cmake_policy(SET CMP0028 NEW)

find_package(HDF5 COMPONENTS CXX HL REQUIRED)
add_subdirectory(unit)
add_subdirectory(moist_thermal)
#add_subdirectory(DYCOMS2D)
