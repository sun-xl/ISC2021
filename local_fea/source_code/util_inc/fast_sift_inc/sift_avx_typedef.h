/*
 * sift_avx_typedef.h
 *
 *  Created on: Nov 13, 2014
 *      Author: peter
 */

#ifndef SIFT_AVX_TYPEDEF_H_
#define SIFT_AVX_TYPEDEF_H_

#include "immintrin.h"  // avx, avx2

#include <stdint.h>

typedef struct {
  union {
    __m256 v;
    float f32[8];
  };
} UV256F32;

typedef struct {
  union {
    __m256i v;
    int i32[8];
  };
} UV256I32;

typedef struct {
  union {
    __m256 v;
    float f32[8];
    unsigned int i32[8];
    unsigned long long int i64[4];
  };
} UV256F32I32I64;

#endif /* SIFT_AVX_TYPEDEF_H_ */
