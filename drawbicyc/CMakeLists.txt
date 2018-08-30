find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)

foreach(exe drawbicyc plot_series_compare average)
  add_executable(${exe} ${exe}.cpp)
  target_link_libraries(${exe} ${HDF5_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${exe} ${Boost_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries(${exe} ${libmpdataxx_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${libcloudphxx_INCLUDE_DIRS})
  target_include_directories(${exe} PUBLIC ${libmpdataxx_INCLUDE_DIRS})
endforeach()

