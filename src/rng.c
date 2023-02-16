#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "rng/rng.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RNG_RAND_MAX 0x7fff

static long lcg_seed;

void rng_srand(long seed) {
  lcg_seed = seed;
}

/* Random numbers */

/* LCG based on the C standard sample implementation */
int rng_int(void) {
  lcg_seed = lcg_seed * 1103515245 + 12345;
  return (unsigned int)(lcg_seed / 65536) % 0x7fff;
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

/* Simplex Noise */

static long simplex_seed;

void rng_simplex_srand(long seed) {
  simplex_seed = seed;
}

#define PRIME_X 0x5205402B9270C86FL
#define PRIME_Y 0x598CD327003817B5L
#define PRIME_Z 0x5BCC226E9FA0BACBL
#define PRIME_W 0x56CC5227E58F554BL
#define HASH_MULTIPLIER 0x53A3F72DEEC546F5L
#define SEED_FLIP_3D -0x52D547B2E96ED629L
#define SEED_OFFSET_4D 0xE83DC3E0DA7164DL

#define ROOT2OVER2 0.7071067811865476f
#define SKEW_2D 0.366025403784439f
#define UNSKEW_2D -0.21132486540518713f

#define ROOT3OVER3 0.577350269189626f
#define FALLBACK_ROTATE_3D (2.0f / 3.0f)
#define ROTATE_3D_ORTHOGONALIZER UNSKEW_2D

#define SKEW_4D -0.138196601125011f
#define UNSKEW_4D 0.309016994374947f
#define LATTICE_STEP_4D 0.2f

#define N_GRADS_2D_EXPONENT 7
#define N_GRADS_3D_EXPONENT 8
#define N_GRADS_4D_EXPONENT 9
#define N_GRADS_2D (1 << N_GRADS_2D_EXPONENT)
#define N_GRADS_3D (1 << N_GRADS_3D_EXPONENT)
#define N_GRADS_4D (1 << N_GRADS_4D_EXPONENT)

#define NORMALIZER_2D 0.01001634121365712f
#define NORMALIZER_3D 0.07969837668935331f
#define NORMALIZER_4D 0.0220065933241897f

#define RSQUARED_2D 0.5f
#define RSQUARED_3D 0.6f
#define RSQUARED_4D 0.6f

static float GRADIENTS_2D[N_GRADS_2D * 2] = {
   0.38268343236509f / NORMALIZER_2D,   0.923879532511287f / NORMALIZER_2D,
   0.923879532511287f / NORMALIZER_2D,  0.38268343236509f / NORMALIZER_2D,
   0.923879532511287f / NORMALIZER_2D, -0.38268343236509f / NORMALIZER_2D,
   0.38268343236509f / NORMALIZER_2D,  -0.923879532511287f / NORMALIZER_2D,
  -0.38268343236509f / NORMALIZER_2D,  -0.923879532511287f / NORMALIZER_2D,
  -0.923879532511287f / NORMALIZER_2D, -0.38268343236509f / NORMALIZER_2D,
  -0.923879532511287f / NORMALIZER_2D,  0.38268343236509f / NORMALIZER_2D,
  -0.38268343236509f / NORMALIZER_2D,   0.923879532511287f / NORMALIZER_2D,
   0.130526192220052f / NORMALIZER_2D,  0.99144486137381f / NORMALIZER_2D,
   0.608761429008721f / NORMALIZER_2D,  0.793353340291235f / NORMALIZER_2D,
   0.793353340291235f / NORMALIZER_2D,  0.608761429008721f / NORMALIZER_2D,
   0.99144486137381f / NORMALIZER_2D,   0.130526192220051f / NORMALIZER_2D,
   0.99144486137381f / NORMALIZER_2D,  -0.130526192220051f / NORMALIZER_2D,
   0.793353340291235f / NORMALIZER_2D, -0.60876142900872f / NORMALIZER_2D,
   0.608761429008721f / NORMALIZER_2D, -0.793353340291235f / NORMALIZER_2D,
   0.130526192220052f / NORMALIZER_2D, -0.99144486137381f / NORMALIZER_2D,
  -0.130526192220052f / NORMALIZER_2D, -0.99144486137381f / NORMALIZER_2D,
  -0.608761429008721f / NORMALIZER_2D, -0.793353340291235f / NORMALIZER_2D,
  -0.793353340291235f / NORMALIZER_2D, -0.608761429008721f / NORMALIZER_2D,
  -0.99144486137381f / NORMALIZER_2D,  -0.130526192220052f / NORMALIZER_2D,
  -0.99144486137381f / NORMALIZER_2D,   0.130526192220051f / NORMALIZER_2D,
  -0.793353340291235f / NORMALIZER_2D,  0.608761429008721f / NORMALIZER_2D,
  -0.608761429008721f / NORMALIZER_2D,  0.793353340291235f / NORMALIZER_2D,
  -0.130526192220052f / NORMALIZER_2D,  0.99144486137381f / NORMALIZER_2D,
};

static int fast_floor(float x) {
  int xi = (int)x;
  return x < xi ? xi - 1 : xi;
}

static float grad(long xsvp, long ysvp, float dx, float dy) {
  long hash = simplex_seed ^ xsvp ^ ysvp;
  hash *= HASH_MULTIPLIER;
  hash ^= hash >> (64 - N_GRADS_2D_EXPONENT + 1);
  int gi = (int)hash & ((N_GRADS_2D - 1) << 1);
  return GRADIENTS_2D[gi | 0] * dx + GRADIENTS_2D[gi | 1] * dy;
}

 /* 2D Simplex noise base */
static float simplex2d_unskewed_base(float xs, float ys) {

  // Get base points and offsets
  int xsb = fast_floor(xs);
  int ysb = fast_floor(ys);
  float xi = (float)(xs - xsb);
  float yi = (float)(ys - ysb);
  
  // Prime pre-multiplication for hash
  long xsbp = xsb * PRIME_X;
  long ysbp = ysb * PRIME_Y;
  
  // Unskew
  float t = (xi + yi) * (float)UNSKEW_2D;
  float dx0 = xi + t, dy0 = yi + t;
  
  // First vertex
  float value = 0;
  float a0 = RSQUARED_2D - dx0 * dx0 - dy0 * dy0;
  if (a0 > 0) {
    value = (a0 * a0) * (a0 * a0) * grad(xsbp, ysbp, dx0, dy0);
  }
  
  // Second vertex
  float a1 = (float)(2 * (1 + 2 * UNSKEW_2D) * (1 / UNSKEW_2D + 2)) * t + ((float)(-2 * (1 + 2 * UNSKEW_2D) * (1 + 2 * UNSKEW_2D)) + a0);
  if (a1 > 0) {
    float dx1 = dx0 - (float)(1 + 2 * UNSKEW_2D);
    float dy1 = dy0 - (float)(1 + 2 * UNSKEW_2D);
    value += (a1 * a1) * (a1 * a1) * grad(xsbp + PRIME_X, ysbp + PRIME_Y, dx1, dy1);
  }
  
  // Third vertex
  if (dy0 > dx0) {
    float dx2 = dx0 - (float)UNSKEW_2D;
    float dy2 = dy0 - (float)(UNSKEW_2D + 1);
    float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
    if (a2 > 0) {
      value += (a2 * a2) * (a2 * a2) * grad(xsbp, ysbp + PRIME_Y, dx2, dy2);
    }
  }
  else {
    float dx2 = dx0 - (float)(UNSKEW_2D + 1);
    float dy2 = dy0 - (float)UNSKEW_2D;
    float a2 = RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
    if (a2 > 0) {
      value += (a2 * a2) * (a2 * a2) * grad(xsbp + PRIME_X, ysbp, dx2, dy2);
    }
  }
  
  return value;
}

float rng_simplex2d(float x, float y) {
  /* Get points for A2* lattice */
  float s = SKEW_2D * (x + y);
  float xs = x + s;
  float ys = y + s;
  
  return simplex2d_unskewed_base(xs, ys);
}

#ifdef __cplusplus
}
#endif
