cmake_policy(SET CMP0057 NEW) # for the in_list operator


# Select flags.
#-g0 negates -g; -UBZ_DEBUG negates blitz debug; NOTE: don't use -Ofast as it makes std::isfinite(sum(x)) return true even if some elements of x are nan/inf
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${libmpdataxx_CXX_FLAGS_DEBUG} -O3 -g0 -UBZ_DEBUG -march=native") 

SET(CMAKE_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE}")
SET(CMAKE_CXX_FLAGS_DEBUG  "${libmpdataxx_CXX_FLAGS_DEBUG}")

add_executable(bicycles opts/opts_common.cpp detail/get_uwlcm_git_revision.cpp bicycles.cpp)

# std::future test
SET(CMAKE_REQUIRED_FLAGS "-std=c++11") # c++11 needed by the std_future test
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <future>\n int main() { std::future<void> f; }\n" STD_FUTURE_WORKS)
if(${STD_FUTURE_WORKS})
  target_compile_definitions(bicycles PRIVATE STD_FUTURE_WORKS)
endif()

# set flags for libmpdata++ options
set(MPDATA_OPTS_PREFIX "MPDATA_OPTS_")
foreach(opts ${MPDATA_OPTS})
  target_compile_definitions(bicycles PRIVATE ${MPDATA_OPTS_PREFIX}${opts})
endforeach()

# handle the disable compilation options
if(UWLCM_DISABLE)
  set(UWLCM_DISABLE_PREFIX "UWLCM_DISABLE_")
  foreach(disable ${UWLCM_DISABLE})
    target_compile_definitions(bicycles PRIVATE ${UWLCM_DISABLE_PREFIX}${disable})
  endforeach()

  # add optional sources
  if(NOT "2D_BLK_1M" IN_LIST UWLCM_DISABLE)
    target_sources(bicycles PRIVATE run_hlpr_2d_blk_1m.cpp)
  endif()
  if(NOT "3D_BLK_1M" IN_LIST UWLCM_DISABLE)
    target_sources(bicycles PRIVATE run_hlpr_3d_blk_1m.cpp)
  endif()
  if(NOT "2D_LGRNGN" IN_LIST UWLCM_DISABLE)
    target_sources(bicycles PRIVATE run_hlpr_2d_lgrngn.cpp)
  endif()
  if(NOT "3D_LGRNGN" IN_LIST UWLCM_DISABLE)
    target_sources(bicycles PRIVATE run_hlpr_3d_lgrngn.cpp)
  endif()
else()
  # if no options were disabled, add all sources
  target_sources(bicycles PRIVATE run_hlpr_2d_blk_1m.cpp run_hlpr_3d_blk_1m.cpp run_hlpr_2d_lgrngn.cpp run_hlpr_3d_lgrngn.cpp)
endif()

# ensure that bicycles is not built before current git_revision file is created
add_dependencies(bicycles git_revision.h) 

target_link_libraries(bicycles PRIVATE ${libmpdataxx_LIBRARIES})
target_include_directories(bicycles PRIVATE ${libmpdataxx_INCLUDE_DIRS})

#search for libcloudphxx
message(STATUS "Looking for libcloudph++")
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

target_link_libraries(bicycles
  PRIVATE
    clphxx::cloudphxx_lgrngn
)

# enabling c++14, but not gnu++14
set_target_properties(bicycles PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(bicycles PRIVATE cxx_std_14)

# search for Boost
find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)
if(TARGET Boost::thread AND TARGET Boost::iostreams AND TARGET Boost::system AND TARGET Boost::timer AND TARGET Boost::program_options AND TARGET Boost::filesystem)
  message("linkiing boost targets")
  target_link_libraries(bicycles PRIVATE Boost::thread Boost::iostreams Boost::system Boost::timer Boost::program_options Boost::filesystem) 
else()
  message("linkiing boost without targets")
  # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
  target_link_libraries(bicycles PRIVATE ${Boost_LIBRARIES})
endif()


install(TARGETS bicycles DESTINATION bin)
