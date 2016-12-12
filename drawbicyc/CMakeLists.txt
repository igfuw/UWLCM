find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(drawbicyc drawbicyc.cpp)
target_link_libraries(drawbicyc ${libmpdataxx_LIBRARIES})
target_include_directories(drawbicyc PUBLIC ${libcloudphxx_INCLUDE_DIRS})
target_include_directories(drawbicyc PUBLIC ${libmpdataxx_INCLUDE_DIRS})

add_executable(plot_series_compare plot_series_compare.cpp)
target_link_libraries(plot_series_compare ${libmpdataxx_LIBRARIES})
target_include_directories(plot_series_compare PUBLIC ${libcloudphxx_INCLUDE_DIRS})
target_include_directories(plot_series_compare PUBLIC ${libmpdataxx_INCLUDE_DIRS})

