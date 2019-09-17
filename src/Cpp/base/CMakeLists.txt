project(base)

set(HEADER_FILES
        algorithm.h
        array.h
        atomic.h
        bits.h
        bitset.h
        coding.h
        compressed.h
        functional.h
        integer.h
        intrinsic.h
        iterator.h
        memory.h
        numeric.h
        process.h
        simd.h
        struct.h
        vector.h
        )

set(SOURCE_FILES
        example_bitset.cpp
        example_bitset.hpp)

add_library(base STATIC ${SOURCE_FILES} ${HEADER_FILES})