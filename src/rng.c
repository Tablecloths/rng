#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "rng/rng.h"

#ifdef __cplusplus
extern "C" {
#endif

static int seed;

void rng_srand(int s) {
  seed = s;
}

int rng_int() {
	return seed = (seed * 1103515245 + 12345);
}

#ifdef __cplusplus
}
#endif
