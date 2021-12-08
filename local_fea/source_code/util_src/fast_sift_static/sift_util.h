/*
 * sift_util.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_UTIL_H_
#define FAST_SIFT_SIFT_UTIL_H_

#include "sift_typedef.h"

// aligned malloc & free
SIFT_EXPORT inline void *sift_aligned_malloc(const unsigned int size,
                                             const unsigned int align);
inline void sift_aligned_free(void **ptr);
inline void *sift_aligned_realloc(void **oldptr, const unsigned int newsize,
                                  const unsigned int align);

inline unsigned int sift_padding(const unsigned int size,
                                 const unsigned int align);

inline int cpu_simd_instruction();

#endif /* FAST_SIFT_SIFT_UTIL_H_ */
