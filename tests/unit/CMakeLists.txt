add_executable(api_test api_test.cpp)
add_test(api_test_iles api_test ${CMAKE_BINARY_DIR})
add_test(api_test_sgs  api_test ${CMAKE_BINARY_DIR} " --sgs=1 ")
