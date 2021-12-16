############################################################################################
# MPI - detecting if the C++ compiler is actually an MPI wrapper

set(msg "Detecting if the compiler is an MPI wrapper...")
message(STATUS "${msg}")
execute_process(COMMAND ${CMAKE_CXX_COMPILER} "-show" RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_QUIET)
if (status EQUAL 0 AND output MATCHES "mpi")
  set(USE_MPI TRUE)
  set(MPIRUN ${CMAKE_CXX_COMPILER})
  string(REPLACE "mpic++" "mpirun" MPIRUN ${MPIRUN})
  string(REPLACE "mpicxx" "mpirun" MPIRUN ${MPIRUN})
  string(REPLACE "mpiXX"  "mpirun" MPIRUN ${MPIRUN})
else()
  set(USE_MPI FALSE)
endif()
message(STATUS "${msg} - ${USE_MPI}")
unset(msg)
unset(status)
unset(output)

############################################################################################

add_executable(moist_thermal moist_thermal.cpp)
target_compile_features(moist_thermal PRIVATE cxx_std_11)
target_link_libraries(moist_thermal PRIVATE UWLCM_plotters)

if(USE_MPI)
  add_test(NAME moist_thermal COMMAND moist_thermal "${MPIRUN} -np 2 ${CMAKE_BINARY_DIR}") 
else()
  add_test(moist_thermal moist_thermal ${CMAKE_BINARY_DIR})
endif()

find_package(Boost COMPONENTS program_options filesystem system REQUIRED)

foreach(exe moist_thermal_diff moist_thermal_create_refdata)
  add_executable(${exe} ${exe}.cpp)
  if(TARGET Boost::system AND TARGET Boost::program_options AND TARGET Boost::filesystem)
    target_link_libraries(${exe} PRIVATE Boost::program_options Boost::filesystem Boost::system)
  else()
    # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
    target_link_libraries(${exe} PRIVATE ${Boost_LIBRARIES})
  endif()
  target_link_libraries(${exe} PRIVATE ${HDF5_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_compile_features(${exe} PRIVATE cxx_std_11)
  
  target_link_libraries(${exe} PRIVATE UWLCM_plotters)
endforeach()

add_test(moist_thermal_diff moist_thermal_diff)


