cmake_minimum_required(VERSION 3.0)
project(rng)

set(RNG_SOURCE_FILES include/rng/rng.h src/rng.c)
source_group("rng" FILES ${RNG_SOURCE_FILES})

add_library(${PROJECT_NAME} ${RNG_SOURCE_FILES})
set_target_properties(${PROJECT_NAME} C_STANDARD 90)
