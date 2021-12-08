/*
 * sift_simd_typedef.h
 *
 *  Created on: Nov 5, 2014
 *      Author: peter
 */

#ifndef SIFT_SIMD_TYPEDEF_H_
#define SIFT_SIMD_TYPEDEF_H_

//#include "xmmintrin.h" // sse
//#include "emmintrin.h" // sse2
//#include "pmmintrin.h" // sse3
//#include "tmmintrin.h" // ssse3
//#include "smmintrin.h" // sse4_1
#include "nmmintrin.h"  // sse4_2
//#include "immintrin.h" // avx, avx2

#include <stdint.h>

typedef struct {
  union {
    float f32;
    int i32;
  };
} UF32I32;

typedef struct {
  union {
    __m128 v;
    float f32[4];
  };
} UV128F32;

typedef struct {
  union {
    __m128i v;
    int i32[4];
  };
} UV128I32;

#endif /* SIFT_SIMD_TYPEDEF_H_ */
