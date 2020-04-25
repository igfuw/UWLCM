
if(NOT MPDATA_OPTS STREQUAL "ABS")
  message(WARNING "The moist_thermal test should be built with the ABS libmpdata++ option, test is most likely to fail") 
endif()

add_executable(moist_thermal moist_thermal.cpp)
target_compile_features(moist_thermal PRIVATE cxx_std_11)

if(USE_MPI)
  add_test(NAME moist_thermal COMMAND moist_thermal "${libmpdataxx_MPIRUN} -np 2 ${CMAKE_BINARY_DIR}") 
else()
  add_test(moist_thermal moist_thermal ${CMAKE_BINARY_DIR})
endif()

add_executable(moist_thermal_diff moist_thermal_diff.cpp)
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
if(TARGET Boost::system AND TARGET Boost::program_options AND TARGET Boost::filesystem)
  target_link_libraries(moist_thermal_diff PRIVATE Boost::program_options Boost::filesystem Boost::system)
else()
  # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
  target_link_libraries(moist_thermal_diff PRIVATE ${Boost_LIBRARIES})
endif()
target_link_libraries(moist_thermal_diff PRIVATE ${HDF5_LIBRARIES})
target_include_directories(moist_thermal_diff PUBLIC ${HDF5_INCLUDE_DIRS})
target_compile_features(moist_thermal_diff PRIVATE cxx_std_11)

add_test(moist_thermal_diff moist_thermal_diff)
