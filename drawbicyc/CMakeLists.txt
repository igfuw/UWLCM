find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_subdirectory(src)

