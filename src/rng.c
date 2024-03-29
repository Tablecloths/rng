#include <math.h>

#include "rng/rng.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RNG_RAND_MAX 0x7fff

static int64_t seed = 0;
static int64_t lcg_seed = 0;

void rng_srand(int64_t s) {
  seed = lcg_seed = s;
}

/* LCG based on the C standard sample implementation */
int rng_int(void) {
  lcg_seed = lcg_seed * 1103515245 + 12345;
  return (unsigned int)(lcg_seed / 65536) % RNG_RAND_MAX;
}

float rng_float(void) {
  return (float)rng_int() / (float)(RNG_RAND_MAX);
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

static const float F2 = 0.36602540378f; /* 0.5*(sqrt(3.0)-1.0) */
static const float G2 = 0.2113248654f;  /* (3.0-sqrt(3.0))/6.0 */
static const float F3 = 0.33333333333f; /* 1.0/3.0 */
static const float G3 = 0.16666666666f; /* 1.0/6.0; */
static const float F4 = 0.30901699437f; /* (sqrt(5.0)-1.0)/4.0; */
static const float G4 = 0.13819660112f; /* (5.0-sqrt(5.0))/20.0; */
static const short PERM[] = { 151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
  190,6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,88,237,149,56,87,174,20,125,136,171,168,68,175,74,
  165,71,134,139,48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,102,143,54,65,25,63,161,
  1,216,80,73,209,76,132,187,208,89,18,169,200,196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,52,217,226,250,124,123,
  5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,223,183,170,213,119,248,152,2,44,154,163,70,221,
  153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113,224,232,178,185,112,104,218,246,97,228,251,34,242,193,238,210,
  144,12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181,199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,
  138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,
  225,140,36,103,30,69,142,8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,88,
  237,149,56,87,174,20,125,136,171,168,68,175,74,165,71,134,139,48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,220,105,
  92,41,55,46,245,40,244,102,143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200,196,135,130,116,188,159,86,164,
  100,109,198,173,186,3,64,52,217,226,250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
  223,183,170,213,119,248,152,2,44,154,163,70,221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113,224,232,178,185,
  112,104,218,246,97,228,251,34,242,193,238,210,144,12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181,199,106,
  157,184,84,204,176,115,121,50,45,127,4,150,254,138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180, };
static const short PERM_MOD_12[] = { 7,4,5,7,6,3,11,1,9,11,0,5,2,5,7,9,8,0,7,6,9,10,8,3,1,0,9,10,11,10,6,4,7,0,6,3,0,2,5,2,10,
  0,3,11,9,11,11,8,9,9,9,4,9,5,8,3,6,8,5,4,3,0,8,7,2,9,11,2,7,0,3,10,5,2,2,3,11,3,1,2,0,7,1,2,4,9,8,5,7,10,5,4,4,6,11,6,5,1,3,
  5,1,0,8,1,5,4,0,7,4,5,6,1,8,4,3,10,8,8,3,2,8,4,1,6,5,6,3,4,4,1,10,10,4,3,5,10,2,3,10,6,3,10,1,8,3,2,11,11,11,4,10,5,2,9,4,6,
  7,3,2,9,11,8,8,2,8,10,7,10,5,9,5,11,11,7,4,9,9,10,3,1,7,2,0,2,7,5,8,4,10,5,4,8,2,6,1,0,11,10,2,1,10,6,0,0,11,11,6,1,9,3,1,7,
  9,2,11,11,1,0,10,7,1,7,10,1,4,0,0,8,7,1,2,9,7,4,6,2,6,8,1,9,6,6,7,5,0,0,3,9,8,3,6,6,11,1,0,0,7,4,5,7,6,3,11,1,9,11,0,5,2,5,
  7,9,8,0,7,6,9,10,8,3,1,0,9,10,11,10,6,4,7,0,6,3,0,2,5,2,10,0,3,11,9,11,11,8,9,9,9,4,9,5,8,3,6,8,5,4,3,0,8,7,2,9,11,2,7,0,3,
  10,5,2,2,3,11,3,1,2,0,7,1,2,4,9,8,5,7,10,5,4,4,6,11,6,5,1,3,5,1,0,8,1,5,4,0,7,4,5,6,1,8,4,3,10,8,8,3,2,8,4,1,6,5,6,3,4,4,1,
  10,10,4,3,5,10,2,3,10,6,3,10,1,8,3,2,11,11,11,4,10,5,2,9,4,6,7,3,2,9,11,8,8,2,8,10,7,10,5,9,5,11,11,7,4,9,9,10,3,1,7,2,0,2,
  7,5,8,4,10,5,4,8,2,6,1,0,11,10,2,1,10,6,0,0,11,11,6,1,9,3,1,7,9,2,11,11,1,0,10,7,1,7,10,1,4,0,0,8,7,1,2,9,7,4,6,2,6,8,1,9,6,
  6,7,5,0,0,3,9,8,3,6,6,11,1,0,0, };
static const float GRAD3[] = { 1.0f,1.0f,0.0f,-1.0f,1.0f,0.0f,1.0f,-1.0f,0.0f,-1.0f,-1.0f,0.0f,1.0f,0.0f,1.0f,-1.0f,0.0f,1.0f,
   1.0f,0.0f,-1.0f,-1.0f,0.0f,-1.0f,0.0f,1.0f,1.0f,0.0f,-1.0f,1.0f,0.0f,1.0f,-1.0f,0.0f,-1.0f,-1.0f, };

static inline int fast_floor(float x) {
  int xi = (int)x;
  return xi - (x < xi);
}

static inline float dot2d(float gx, float gy, float x, float y) {
  return gx * x + gy * y;
}

float rng_simplex2d(float xin, float yin) {
  float n0, n1, n2;
  float s = (xin + yin) * F2;
  int i = fast_floor(xin + s);
  int j = fast_floor(yin + s);
  float t = (i + j) * G2;
  float x0 = xin - (i - t);
  float y0 = yin - (j - t);
  int i1, j1; 
  if(x0 > y0) {
    i1 = 1; j1 = 0;
  }
  else {
    i1 = 0; j1 = 1;
  }
  float x1 = x0 - i1 + G2;
  float y1 = y0 - j1 + G2;
  float x2 = x0 - 1.0f + 2.0f * G2;
  float y2 = y0 - 1.0f + 2.0f * G2;
  int ii = (i + seed) & 255;
  int jj = (j + seed) & 255;
  int gi0 = PERM_MOD_12[ii + PERM[jj]];
  int gi1 = PERM_MOD_12[ii + i1 + PERM[jj + j1]];
  int gi2 = PERM_MOD_12[ii + 1 + PERM[jj + 1]];
  float t0 = 0.5f - x0 * x0 - y0 * y0;
  if (t0 < 0) {
    n0 = 0.0;
  }
  else {
    t0 *= t0;
    n0 = t0 * t0 * dot2d(GRAD3[0 + gi0 * 3], GRAD3[1 + gi0 * 3], x0, y0);
  }
  float t1 = 0.5f - x1 * x1 - y1 * y1;
  if (t1 < 0) {
    n1 = 0.0;
  }
  else {
    t1 *= t1;
    n1 = t1 * t1 * dot2d(GRAD3[0 + gi1 * 3], GRAD3[1 + gi1 * 3], x1, y1);
  }
  float t2 = 0.5f - x2 * x2 - y2 * y2;
  if (t2 < 0) {
    n2 = 0.0f;
  }
  else {
    t2 *= t2;
    n2 = t2 * t2 * dot2d(GRAD3[0 + gi2 * 3], GRAD3[1 + gi2 * 3], x2, y2);
  }
  return 70.0f * (n0 + n1 + n2);
}

#ifdef __cplusplus
}
#endif
