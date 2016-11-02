cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(bicycles CXX C)

enable_testing()

#search for libcloudphxx
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

add_subdirectory(src) 
add_subdirectory(drawbicyc)
add_subdirectory(tests)
