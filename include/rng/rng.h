#ifndef LIBRNG_INCLUDE_RNG_H_
#define LIBRNG_INCLUDE_RNG_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Provide a seed for random number generation */
void rng_srand(int seed);

/* Returns a pseudorandom integer between 0 and RNG_RAND_MAX (32767) */
int rng_int();

/* Returns a normally distributed random number with a given arithmetic mean and standard deviation */
float rng_normal(float mean, float std);

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
