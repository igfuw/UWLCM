# Select flags.
#-g0 negates -g; -UBZ_DEBUG negates blitz debug; NOTE: don't use -Ofast as it makes std::isfinite(sum(x)) return true even if some elements of x are nan/inf
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${libmpdataxx_CXX_FLAGS_DEBUG} -O3 -g0 -UBZ_DEBUG -march=native") 

SET(CMAKE_CXX_FLAGS_RELEASE "${libmpdataxx_CXX_FLAGS_RELEASE}")
SET(CMAKE_CXX_FLAGS_DEBUG  "${libmpdataxx_CXX_FLAGS_DEBUG}")

add_executable(bicycles bicycles.cpp)

# std::future
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
endif()


# ensure that bicycles is not built before current git_revision file is created
add_dependencies(bicycles git_revision.h) 

target_link_libraries(bicycles PRIVATE ${libmpdataxx_LIBRARIES})
target_include_directories(bicycles PRIVATE ${libmpdataxx_INCLUDE_DIRS})

target_link_libraries(bicycles
  PRIVATE
    clphxx::cloudphxx_lgrngn
)

# enabling c++14, but not gnu++14
set_target_properties(bicycles PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(bicycles PRIVATE cxx_std_14)

find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)
target_link_libraries(bicycles PRIVATE ${Boost_LIBRARIES})

install(TARGETS bicycles DESTINATION bin)
