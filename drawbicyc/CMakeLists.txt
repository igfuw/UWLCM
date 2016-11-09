find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(drawbicyc drawbicyc.cpp)
target_link_libraries(drawbicyc ${HDF5_LIBRARIES})
target_include_directories(drawbicyc PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(drawbicyc ${Boost_LIBRARIES})
target_include_directories(drawbicyc PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(drawbicyc PUBLIC ${libcloudphxx_INCLUDE_DIRS})

