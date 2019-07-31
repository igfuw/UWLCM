cmake_minimum_required(VERSION 2.8.6) # not sure, simply copied from UWLCM
find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)

############################################################################################
#libcloudphxx
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

############################################################################################
# Blitz++
find_package(PkgConfig)
pkg_check_modules(BLITZ QUIET blitz>=0.10)
if (NOT BLITZ_FOUND)
  message(FATAL_ERROR "Blitz++ library not found or Blitz++ version requirement not met (>=0.10)

* To insall Blitz++, please try:
*   Debian/Ubuntu: sudo apt-get install libblitz0-dev
*   Fedora: sudo yum install blitz-devel
*   Homebrew: brew install blitz
  ")
endif()
############################################################################################


foreach(exe drawbicyc plot_series_compare average)
  add_executable(${exe} ${exe}.cpp)
  target_link_libraries(${exe} PRIVATE ${HDF5_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${exe} PRIVATE ${Boost_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries(${exe} PRIVATE clphxx::cloudphxx_lgrngn)
  target_link_libraries(${exe} PRIVATE ${libmpdataxx_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${libmpdataxx_INCLUDE_DIRS})
  target_link_libraries(${exe} PRIVATE ${BLITZ_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${BLITZ_INCLUDE_DIRS})
endforeach()

