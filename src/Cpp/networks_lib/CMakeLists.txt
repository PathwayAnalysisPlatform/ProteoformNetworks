project(networks_lib)

set(HEADER_FILES
        bimap_str_int.hpp
        scores.hpp
        types.hpp
        maps.hpp
        )

set(SOURCE_FILES
        bimap_str_int.cpp
        scores.cpp
        types.cpp)

add_library(networks_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})