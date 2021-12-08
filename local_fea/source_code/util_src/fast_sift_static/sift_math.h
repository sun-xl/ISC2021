/*
 * sift_math.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_MATH_H_
#define FAST_SIFT_SIFT_MATH_H_

#include <math.h>

// fast_expn
#define EXPN_SZ 256    // fast_expn table size
#define EXPN_MAX 25.0  // fast_expn table max
#define EXPN_SZ_EXPN_MAX (EXPN_SZ / EXPN_MAX)
inline float *fast_expn_init();
inline float fast_expn_f(float x);

inline float sift_fast_sqrt_f(const float x);
inline float sift_fast_atan2_f(float y, float x);

inline float sift_mod_2pi_f(float x);

inline long int sift_floor_f(const float x);

// EPSILON_F is equal to 2^{-23}
#define EPSILON_F 1.19209290E-07F

// Compute the minimum between two values
#define SIFT_MIN(x, y) (((x) < (y)) ? (x) : (y))
// Compute the maximum between two values
#define SIFT_MAX(x, y) (((x) > (y)) ? (x) : (y))

// Signed left shift operation
#define SIFT_SHIFT_LEFT(x, n) (((n) >= 0) ? ((x) << (n)) : ((x) >> -(n)))

#define PI 3.141592653589793f
#define PI2 (PI * 2)

#endif /* FAST_SIFT_SIFT_MATH_H_ */
