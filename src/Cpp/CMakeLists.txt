cmake_minimum_required(VERSION 3.14)
project(ProteoformNetworks)

set(CMAKE_CXX_STANDARD 17)
set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

set(SOURCE_FILES overlap_analysis.cpp overlap_analysis.hpp
        create_modules.hpp create_modules.cpp
        Interactome.hpp Interactome.cpp
        Module.cpp Module.hpp
        networks_lib/types.hpp networks_lib/bimap_str_int.hpp networks_lib/bimap_str_int.cpp networks_lib/types.cpp)
add_executable(ProteoformNetworks_run ${SOURCE_FILES})

include_directories(networks_lib)
add_subdirectory(networks_lib)
target_link_libraries(ProteoformNetworks_run networks_lib)

include_directories(tools)
add_subdirectory(tools)
target_link_libraries(ProteoformNetworks_run tools)

include_directories(base)
add_subdirectory(base)
target_link_libraries(ProteoformNetworks_run base)

add_subdirectory(Google_tests)