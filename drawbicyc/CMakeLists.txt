cmake_minimum_required(VERSION 2.8.6) # not sure, simply copied from UWLCM
find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)

#search for libcloudphxx
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

foreach(exe drawbicyc plot_series_compare average)
  add_executable(${exe} ${exe}.cpp)
  target_link_libraries(${exe} ${HDF5_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${exe} ${Boost_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries(${exe} ${libmpdataxx_LIBRARIES})
  target_link_libraries(${exe} PRIVATE clphxx::cloudphxx_lgrngn)
  target_include_directories(${exe} PUBLIC ${libmpdataxx_INCLUDE_DIRS})
endforeach()

