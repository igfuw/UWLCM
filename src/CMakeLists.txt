# search for libmpdata
find_package(libmpdata++)
if (NOT libmpdataxx_FOUND)
  message(FATAL_ERROR "libmpdata++ not found")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  # make the Release flags the default
  set(CMAKE_CXX_FLAGS ${libmpdataxx_CXX_FLAGS_RELEASE})
  # ignore CMake defaults
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "${libmpdataxx_CXX_FLAGS_DEBUG} -O1")
endif()

# std::future
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <future>\n int main() { std::future<void> f; }\n" STD_FUTURE_WORKS)
if(${STD_FUTURE_WORKS})
  add_definitions(-DSTD_FUTURE_WORKS)
endif()

add_executable(bicycles bicycles.cpp)

target_link_libraries(bicycles ${libmpdataxx_LIBRARIES})
target_include_directories(bicycles PUBLIC ${libmpdataxx_INCLUDE_DIRS})

target_link_libraries(bicycles ${libcloudphxx_LIBRARIES})
target_include_directories(bicycles PUBLIC ${libcloudphxx_INCLUDE_DIRS})

find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)
target_link_libraries(bicycles ${Boost_LIBRARIES})

install(TARGETS bicycles DESTINATION bin)
