add_executable(drawbicyc drawbicyc.cpp)
target_link_libraries(drawbicyc ${HDF5_LIBRARIES})
target_include_directories(drawbicyc PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(drawbicyc ${Boost_LIBRARIES})
target_include_directories(drawbicyc PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(drawbicyc PUBLIC ${libcloudphxx_INCLUDE_DIRS})

