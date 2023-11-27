#ifndef LIBRNG_INCLUDE_RNG_H_
#define LIBRNG_INCLUDE_RNG_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Provide a seed for random number generation */
void rng_srand(int64_t seed);

/* Returns a pseudorandom integer between 0 and RNG_RAND_MAX (32767) */
int rng_int(void);

/* Returns a pseudorandom float evenly distributed between 0 and 1 */
float rng_float(void);

/* Returns a normally distributed random number with a given arithmetic mean and standard deviation */
float rng_normal(float mean, float std);

/* Provide a seed for simplex-like noise generation */
void rng_simplex_srand(int64_t seed);

/* Samples 2-dimensional simplex-like noise at a given location */
float rng_simplex2d(float x, float y);

#ifdef LIBRNG_IMPLEMENTATION
#undef LIBRNG_IMPLEMENTATION
/* Prevent tools like dependency checkers that don't evaluate
   macros from detecting a cyclic dependency. */
#define LIBRNG_SOURCE "rng.c"
#include LIBRNG_SOURCE
#endif /* LIBRNG_IMPLEMENTATION */

#ifdef __cplusplus
}
#endif

#endif /* LIBRNG_INCLUDE_RNG_H_ */
