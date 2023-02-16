cmake_minimum_required(VERSION 3.0)
project(rng)

set(RNG_SOURCE_FILES include/rng/rng.h src/rng.c)
source_group("rng" FILES ${RNG_SOURCE_FILES})

add_library(${PROJECT_NAME} ${RNG_SOURCE_FILES})
target_compile_features(${PROJECT_NAME} PRIVATE c_std_90)
target_include_directories(${PROJECT_NAME} PUBLIC include)
target_include_directories(${PROJECT_NAME} PUBLIC src)
