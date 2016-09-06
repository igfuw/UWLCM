add_executable(unit_2d unit_2d.cpp)
add_executable(unit_3d unit_3d.cpp)
add_test(unit_2d unit_2d ${CMAKE_BINARY_DIR})
add_test(unit_3d unit_3d ${CMAKE_BINARY_DIR})
