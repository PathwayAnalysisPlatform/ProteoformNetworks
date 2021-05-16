cmake_minimum_required(VERSION 3.14)
project(ProteoformNetworks)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp)
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