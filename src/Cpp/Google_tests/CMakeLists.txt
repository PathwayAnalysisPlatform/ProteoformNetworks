project(Google_tests)
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# adding the Google_Tests_run target
add_executable(Google_Tests_run conversions_tests.cpp
        bimap_tests.cpp
        phegeni_tests.cpp
        scores_tests.cpp
        networks_tests.cpp
        )

# linking Google_Tests_run with DateConverter_lib which will be tested
target_link_libraries(Google_Tests_run networks_lib)

target_link_libraries(Google_Tests_run gtest gtest_main)