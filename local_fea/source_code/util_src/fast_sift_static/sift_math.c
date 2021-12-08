/*
 * sift_math.c
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#include "sift_math.h"

int g_expn_tab_init = 0;  // fast_expn table initialize flag
float g_expn_tab[EXPN_SZ + 1]
    __attribute__((aligned(32)));  // fast_expn table, 0~25.0
// static int g_expn_tab_init = 0; // fast_expn table initialize flag
// static float g_expn_tab[EXPN_SZ+1] __attribute__((aligned (32))); //
// fast_expn table, 0~25.0

inline float *fast_expn_init() {
  int k;
  if (g_expn_tab_init) return g_expn_tab;

  g_expn_tab_init = 1;
  for (k = 0; k < EXPN_SZ + 1; ++k)                          // 0~256
    g_expn_tab[k] = expf(-(float)k * (EXPN_MAX / EXPN_SZ));  // 0~25.0
  return g_expn_tab;
}

inline float fast_expn_f(float x) {
  float a, b, r;
  int i;

  x *= EXPN_SZ / EXPN_MAX;

  // floor
  int xi = (int)x;
  if (x >= 0 || (float)xi == x)
    i = xi;
  else
    i = xi - 1;

  // linear interpolation
  r = x - i;
  a = g_expn_tab[i];
  b = g_expn_tab[i + 1];
  return a + r * (b - a);
}

// http://www.lomont.org/Math/Papers/2003/InvSqrt.pdf
inline float sift_fast_resqrt_f(const float x) {
  /* 32-bit version */
  union {
    float x;
    int i;
  } u;

  float xhalf = (float)0.5 * x;

  /* convert floating point value in RAW integer */
  u.x = x;

  /* gives initial guess y0 */
  u.i = 0x5f3759df - (u.i >> 1);
  /*u.i = 0xdf59375f - (u.i>>1);*/

  /* two Newton steps */
  u.x = u.x * ((float)1.5 - xhalf * u.x * u.x);
  u.x = u.x * ((float)1.5 - xhalf * u.x * u.x);
  return u.x;
}

// SIFT_INLINE float sift_fast_sqrt_f(float x)
inline float sift_fast_sqrt_f(const float x) {
  return (x < 1e-8) ? 0 : x * sift_fast_resqrt_f(x);
}

// f(r)=c_0 + c_1 r + c_2 r^2 + c_3 r^3
// To fit the polynomial we impose the constraints
// f(+1) &=& c_0 + c_1 + c_2 + c_3  = atan(0)       = 0,
// f(-1) &=& c_0 - c_1 + c_2 - c_3  = atan(\infty)  = \pi/2,
// f(0)  &=& c_0                    = atan(1)       = \pi/4.
// c_0=\pi/4,
// c_1=-0.9675,
// c_2=0,
// c_3=0.1821,
inline float sift_fast_atan2_f(float y, float x) {
  float angle, r;
  float const c3 = 0.1821F;
  float const c1 = 0.9675F;
  float abs_y = fabsf(y) + EPSILON_F;

  if (x >= 0) {
    r = (x - abs_y) / (x + abs_y);
    angle = (float)(PI / 4);
  } else {
    r = (x + abs_y) / (abs_y - x);
    angle = (float)(3 * PI / 4);  // +PI/2
  }
  angle += (c3 * r * r - c1) * r;
  return (y < 0) ? -angle : angle;
}

inline float sift_mod_2pi_f(float x) {
  while (x > (float)(2 * PI)) x -= (float)(2 * PI);
  while (x < 0.0F) x += (float)(2 * PI);
  return x;
}

inline long int sift_floor_f(const float x) {
  long int xi = (long int)x;
  if (x >= 0 || (float)xi == x)
    return xi;
  else
    return xi - 1;
}
