# rng
A C89 header-only pseudorandom number library

## Features
* LCG for integers and floats
* Random Gaussians
* 2D Simplex-like noise

## Building
There are multiple ways to use rng in your project:

* You can just add rng.c to your build system
* You can use provided CMake files
* You can use rng in header-only fashion. In exactly one source file, define RNG_IMPLEMENTATION before including rng.h. Do not build rng.c at all in this case.

