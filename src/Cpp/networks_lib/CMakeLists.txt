project(networks_lib)

set(HEADER_FILES
        reactome.hpp
        bimap_str_int.hpp
        conversions.hpp
        scores.hpp
        phegeni.hpp
        types.hpp
        )

set(SOURCE_FILES
        reactome.cpp
        bimap_str_int.cpp
        conversions.cpp
        scores.cpp
        phegeni.cpp
        networks.hpp networks.cpp maps.hpp maps.cpp)

add_library(networks_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})