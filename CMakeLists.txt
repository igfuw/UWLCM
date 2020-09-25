cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
cmake_policy(SET CMP0057 NEW) # for the in_list operator

project(uwlcm CXX C)

#enable_testing()

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
  add_custom_target(git_revision.h ALL
   git log -1 "--format=format:#define UWLCM_GIT_REVISION \"%H\"%n" HEAD > include/UWLCM/git_revision.h
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} VERBATIM
  )
endif()

# handle the mpdata options
if(NOT MPDATA_OPTS)
  set(MPDATA_OPTS "IGA;FCT") # iga, fct is the default
endif()
message(STATUS "UWLCM will use the following libmpdata++ options: ${MPDATA_OPTS}")

# search for libmpdata
# done here to set USE_MPI if it is a mpi build
message(STATUS "Looking for libmpdata++")
find_package(libmpdata++)
if (NOT libmpdataxx_FOUND)
  message(FATAL_ERROR "libmpdata++ not found")
endif()

# Select flags.
#-g0 negates -g; -UBZ_DEBUG negates blitz debug; NOTE: don't use -Ofast as it makes std::isfinite(sum(x)) return true even if some elements of x are nan/inf
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${libmpdataxx_CXX_FLAGS_DEBUG} -O3 -g0 -UBZ_DEBUG -march=native") 

SET(CMAKE_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE}")
SET(CMAKE_CXX_FLAGS_DEBUG  "${libmpdataxx_CXX_FLAGS_DEBUG}")

add_executable(uwlcm src/opts/opts_common.cpp src/detail/get_uwlcm_git_revision.cpp src/uwlcm.cpp)
add_library(uwlcm_includes INTERFACE)

# pointing to UWLCM headers
target_include_directories(uwlcm
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)
target_include_directories(uwlcm_includes
  INTERFACE
    $<INSTALL_INTERFACE:include>
)

# std::future test
SET(CMAKE_REQUIRED_FLAGS "-std=c++11") # c++11 needed by the std_future test
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <future>\n int main() { std::future<void> f; }\n" STD_FUTURE_WORKS)
if(${STD_FUTURE_WORKS})
  target_compile_definitions(uwlcm PRIVATE STD_FUTURE_WORKS)
endif()

# set flags for libmpdata++ options
set(MPDATA_OPTS_PREFIX "MPDATA_OPTS_")
foreach(opts ${MPDATA_OPTS})
  target_compile_definitions(uwlcm PRIVATE ${MPDATA_OPTS_PREFIX}${opts})
endforeach()

# handle the disable compilation options
if(UWLCM_DISABLE)
  set(UWLCM_DISABLE_PREFIX "UWLCM_DISABLE_")
  foreach(disable ${UWLCM_DISABLE})
    target_compile_definitions(uwlcm PRIVATE ${UWLCM_DISABLE_PREFIX}${disable})
  endforeach()

  # add optional sources
  if(NOT "2D_BLK_1M" IN_LIST UWLCM_DISABLE)
    target_sources(uwlcm PRIVATE src/run_hlpr_2d_blk_1m.cpp)
  endif()
  if(NOT "3D_BLK_1M" IN_LIST UWLCM_DISABLE)
    target_sources(uwlcm PRIVATE src/run_hlpr_3d_blk_1m.cpp)
  endif()
  if(NOT "2D_BLK_2M" IN_LIST UWLCM_DISABLE)
    target_sources(uwlcm PRIVATE src/run_hlpr_2d_blk_2m.cpp)
  endif()
  if(NOT "3D_BLK_2M" IN_LIST UWLCM_DISABLE)
    target_sources(uwlcm PRIVATE src/run_hlpr_3d_blk_2m.cpp)
  endif()
  if(NOT "2D_LGRNGN" IN_LIST UWLCM_DISABLE)
    target_sources(uwlcm PRIVATE src/run_hlpr_2d_lgrngn.cpp)
  endif()
  if(NOT "3D_LGRNGN" IN_LIST UWLCM_DISABLE)
    target_sources(uwlcm PRIVATE src/run_hlpr_3d_lgrngn.cpp)
  endif()
else()
  # if no options were disabled, add all sources
  target_sources(uwlcm PRIVATE src/run_hlpr_2d_blk_1m.cpp src/run_hlpr_3d_blk_1m.cpp src/run_hlpr_2d_blk_2m.cpp src/run_hlpr_3d_blk_2m.cpp src/run_hlpr_2d_lgrngn.cpp src/run_hlpr_3d_lgrngn.cpp)
endif()

# ensure that uwlcm is not built before current git_revision file is created
add_dependencies(uwlcm git_revision.h)

target_link_libraries(uwlcm PRIVATE ${libmpdataxx_LIBRARIES})
target_include_directories(uwlcm PRIVATE ${libmpdataxx_INCLUDE_DIRS})

#search for libcloudphxx
message(STATUS "Looking for libcloudph++")
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

target_link_libraries(uwlcm
  PRIVATE
    clphxx::cloudphxx_lgrngn
)

# enabling c++14, but not gnu++14
set_target_properties(uwlcm PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(uwlcm PRIVATE cxx_std_14)

# search for Boost
find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)
if(TARGET Boost::thread AND TARGET Boost::iostreams AND TARGET Boost::system AND TARGET Boost::timer AND TARGET Boost::program_options AND TARGET Boost::filesystem)
  message("linkiing boost targets")
  target_link_libraries(uwlcm PRIVATE Boost::thread Boost::iostreams Boost::system Boost::timer Boost::program_options Boost::filesystem)
else()
  message("linkiing boost without targets")
  # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
  target_link_libraries(uwlcm PRIVATE ${Boost_LIBRARIES})
endif()



#add_subdirectory(tests)
add_subdirectory(include) 

############################################################################################
# installation

install(TARGETS uwlcm uwlcm_includes EXPORT uwlcm-targets 
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include)

install(EXPORT uwlcm-targets NAMESPACE uwlcm:: DESTINATION share/UWLCM)

# generate and install a config file
include(CMakePackageConfigHelpers)

configure_package_config_file(
  uwlcm-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/uwlcm-config.cmake"
  INSTALL_DESTINATION
    share/UWLCM
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/uwlcm-config.cmake"
  DESTINATION
    share/UWLCM
)
############################################################################################

