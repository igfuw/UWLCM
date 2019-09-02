if(NOT MPDATA_OPTS STREQUAL "ABS")
  message(WARNING "The moist_thermal test should be built with the ABS libmpdata++ option, test is most likely to fail") 
endif()

find_package(Boost COMPONENTS program_options filesystem system REQUIRED)

#regular moist thermal test
add_executable(moist_thermal moist_thermal.cpp)
target_compile_features(moist_thermal PRIVATE cxx_std_11)
target_link_libraries(moist_thermal PRIVATE ${HDF5_LIBRARIES})
target_include_directories(moist_thermal PUBLIC ${HDF5_INCLUDE_DIRS})
if(TARGET Boost::system AND TARGET Boost::program_options AND TARGET Boost::filesystem)
  target_link_libraries(moist_thermal PRIVATE Boost::program_options Boost::filesystem Boost::system)
else()
  # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
  target_link_libraries(moist_thermal PRIVATE ${Boost_LIBRARIES})
endif()

if(USE_MPI)
  add_test(NAME moist_thermal COMMAND moist_thermal "${libmpdataxx_MPIRUN} -np 2 ${CMAKE_BINARY_DIR}") 
else()
  add_test(moist_thermal moist_thermal ${CMAKE_BINARY_DIR})
endif()

# simplified moist thermal test for MPI tests on travis
#add_executable(moist_thermal_travis EXCLUDE_FROM_ALL moist_thermal_travis.cpp)
#target_compile_features(moist_thermal_travis PRIVATE cxx_std_11)
#target_link_libraries(moist_thermal_travis PRIVATE ${HDF5_LIBRARIES})
#target_include_directories(moist_thermal_travis PUBLIC ${HDF5_INCLUDE_DIRS})
#if(TARGET Boost::system AND TARGET Boost::program_options AND TARGET Boost::filesystem)
#  target_link_libraries(moist_thermal_travis PRIVATE Boost::program_options Boost::filesystem Boost::system)
#else()
#  # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
#  target_link_libraries(moist_thermal_travis PRIVATE ${Boost_LIBRARIES})
#endif()
#
#if(USE_MPI)
#  add_test(NAME moist_thermal_travis COMMAND moist_thermal_travis "${libmpdataxx_MPIRUN} -np 2 ${CMAKE_BINARY_DIR}") 
#else()
#  add_test(moist_thermal_travis moist_thermal_travis ${CMAKE_BINARY_DIR})
#endif()
