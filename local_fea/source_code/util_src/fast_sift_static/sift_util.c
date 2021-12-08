/*
 * sift_util.c
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#include "sift_util.h"
#include "sift_typedef.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

static int g_optimize_type;
static int g_optimize_type_init = 0;

// aligned malloc & free
SIFT_EXPORT inline void *sift_aligned_malloc(const unsigned int size,
                                             const unsigned int align) {
  const unsigned int size_real = size + align + 3;
  unsigned char *ptr = (unsigned char *)calloc(size_real, 1);
  unsigned int align_offset4, offset;
  if (ptr == NULL) {
    fprintf(stderr, "sift_aligned_malloc: fail to allocate %d bytes\n",
            size_real);
    return NULL;
  }
  align_offset4 = ((uintptr_t)(ptr + 4)) % align;
  offset = align_offset4 == 0 ? 4 : 4 + align - align_offset4;
  ptr += offset;
  *(unsigned int *)(ptr - 4) = offset;
  return ptr;
}

inline void sift_aligned_free(void **ptr) {

  if ((*ptr) != NULL) {
    unsigned char *p = (unsigned char *)(*ptr);
    free(p - *(unsigned int *)(p - 4));
    (*ptr) = NULL;
  }
}

inline void *sift_aligned_realloc(void **oldptr, const unsigned int newsize,
                                  const unsigned int align) {
  unsigned char *oldptr_real = (unsigned char *)(*oldptr);
  oldptr_real -= *(unsigned int *)(oldptr_real - 4);
  (*oldptr) = NULL;

  const unsigned int size_real = newsize + align + 3;
  unsigned char *ptr = (unsigned char *)realloc(oldptr_real, size_real);
  unsigned int align_offset4, offset;
  if (ptr == NULL) {
    fprintf(stderr, "sift_aligned_realloc: fail to allocate %d bytes\n",
            size_real);
    return NULL;
  }
  align_offset4 = ((uintptr_t)(ptr + 4)) % align;
  offset = align_offset4 == 0 ? 4 : 4 + align - align_offset4;
  ptr += offset;
  *(unsigned int *)(ptr - 4) = offset;
  return ptr;
}

// padding
inline unsigned int sift_padding(const unsigned int size,
                                 const unsigned int align) {
  const unsigned int tailing = size % align;
  const unsigned int padding = tailing == 0 ? 0 : align - tailing;
  return (size + padding);
}

// optimization type detection
inline int cpu_simd_instruction() {
  int reg_ebx, reg_ecx, reg_edx;

  if (g_optimize_type_init) return g_optimize_type;

  g_optimize_type_init = 1;
  g_optimize_type = NO_OPTIMIZE;
  __asm__ __volatile__(
      "movl $1, %%eax\n\t"
      "cpuid\n\t"
      "movl %%ebx, %0\n\t"
      "movl %%ecx, %1\n\t"
      "movl %%edx, %2\n\t"
      : "=r"(reg_ebx), "=r"(reg_ecx), "=r"(reg_edx)
      :
      : "%eax", "%ebx", "%ecx", "%edx");
  if (reg_edx & 0x2000000) g_optimize_type |= OPTIMIZE_SSE;
  if (reg_edx & 0x4000000) g_optimize_type |= OPTIMIZE_SSE2;
  if (reg_ecx & 0x01) g_optimize_type |= OPTIMIZE_SSE3;
  if (reg_ecx & 0x200) g_optimize_type |= OPTIMIZE_SSSE3;
  if (reg_ecx & 0x80000) g_optimize_type |= OPTIMIZE_SSE4_1;
  if (reg_ecx & 0x100000) {
    g_optimize_type |= OPTIMIZE_SSE4_2;

    if (reg_ecx & 0x18000000) g_optimize_type |= OPTIMIZE_AVX;
  }

  __asm__ __volatile__(
      "movl $7, %%eax\n\t"
      "movl $0, %%ecx\n\t"
      "cpuid\n\t"
      "movl %%ebx, %0\n\t"
      "movl %%ecx, %1\n\t"
      "movl %%edx, %2\n\t"
      : "=r"(reg_ebx), "=r"(reg_ecx), "=r"(reg_edx)
      :
      : "%eax", "%ebx", "%ecx", "%edx");
  if (reg_ebx & 0x20) g_optimize_type |= OPTIMIZE_AVX2;

  fprintf(stderr, "\nSIMD Instruction: ");
  if (g_optimize_type) {
    if (g_optimize_type & OPTIMIZE_SSE) fprintf(stderr, "SSE ");
    if (g_optimize_type & OPTIMIZE_SSE2) fprintf(stderr, "SSE2 ");
    if (g_optimize_type & OPTIMIZE_SSE3) fprintf(stderr, "SSE3 ");
    if (g_optimize_type & OPTIMIZE_SSSE3) fprintf(stderr, "SSSE3 ");
    if (g_optimize_type & OPTIMIZE_SSE4_1) fprintf(stderr, "SSE4_1 ");
    if (g_optimize_type & OPTIMIZE_SSE4_2) fprintf(stderr, "SSE4_2 ");
    if (g_optimize_type & OPTIMIZE_AVX) fprintf(stderr, "AVX ");
    if (g_optimize_type & OPTIMIZE_AVX2) fprintf(stderr, "AVX2 ");
  } else
    fprintf(stderr, "NONE ");
  fprintf(stderr, "\n");

  return g_optimize_type;
}
