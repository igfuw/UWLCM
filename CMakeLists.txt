cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(bicycles CXX C)

enable_testing()

add_subdirectory(src) 
add_subdirectory(tests)
