if(NOT MPDATA_OPTS STREQUAL "ABS")
  message(WARNING "The moist_thermal test should be built with the ABS libmpdata++ option, test is most likely to fail") 
endif()

add_executable(moist_thermal moist_thermal.cpp)
target_compile_features(moist_thermal PRIVATE cxx_std_11)
add_test(moist_thermal moist_thermal ${CMAKE_BINARY_DIR})

target_link_libraries(moist_thermal PRIVATE ${HDF5_LIBRARIES})
target_include_directories(moist_thermal PUBLIC ${HDF5_INCLUDE_DIRS})

find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
target_link_libraries(moist_thermal PRIVATE Boost::program_options Boost::filesystem Boost::system)

