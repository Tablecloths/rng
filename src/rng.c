#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "rng/rng.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RNG_RAND_MAX 0x7fff

static int seed;

void rng_srand(int s) {
  seed = s;
}

/* LCG based on the C standard sample implementation */
int rng_int(void) {
  seed = seed * 1103515245 + 12345;
  return (unsigned int)(seed / 65536) % 0x7fff;
}

float rng_float(void) {
  return (float)rng_int()/(float)(RNG_RAND_MAX);
}

/* Box-Muller */
float rng_normal(float mean, float std) {
  static float cached = 0.0f;
  float x, y, r, res;
  
  if (cached == 0.0f) {
    do {
      x = 2.0f * rng_int() / RNG_RAND_MAX - 1;
      y = 2.0f * rng_int() / RNG_RAND_MAX - 1;
      r = x * x + y * y;
    } while (r == 0.0f || r > 1.0f);
    
    float d = (float)sqrt(-2.0f * (float)log(r) / r);
    
    float n1 = x * d;
    float n2 = y * d;
    
    res = n1 * std + mean;
    cached = n2;
  }
  else {
    res = cached * std + mean;
    cached = 0.0f;
  }
  
  return res;
}

#ifdef __cplusplus
}
#endif
