#ifndef LIBRNG_INCLUDE_RNG_H_
#define LIBRNG_INCLUDE_RNG_H_

#ifdef __cplusplus
extern "C" {
#endif

void rng_srand(int seed);
int rng_int();

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
