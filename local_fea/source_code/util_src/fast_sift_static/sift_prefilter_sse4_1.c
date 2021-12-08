/*
 * sift_prefilter.c
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#include "sift_prefilter.h"
#include "sift_util.h"
#include "sift_math.h"
#include "sift_sse4_2_typedef.h"
#include <string.h>

SIFT_EXPORT void sharp_filter_sse4_1(const unsigned char *const input,
                                     unsigned char *const output, const int w,
                                     const int h) {
  int i, j;
  const unsigned short kernel[13 * 16] __attribute__((aligned(32))) = {
      13, 18, 24, 30, 35,  39,  40,  39,  35,  30, 24, 18, 13, 0, 0, 0,
      18, 26, 34, 43, 50,  55,  57,  55,  50,  43, 34, 26, 18, 0, 0, 0,
      24, 34, 46, 57, 66,  73,  75,  73,  66,  57, 46, 34, 24, 0, 0, 0,
      30, 43, 57, 70, 82,  91,  93,  91,  82,  70, 57, 43, 30, 0, 0, 0,
      35, 50, 66, 82, 96,  106, 109, 106, 96,  82, 66, 50, 35, 0, 0, 0,
      39, 55, 73, 91, 106, 116, 120, 116, 106, 91, 73, 55, 39, 0, 0, 0,
      40, 57, 75, 93, 109, 120, 124, 120, 109, 93, 75, 57, 40, 0, 0, 0,
      39, 55, 73, 91, 106, 116, 120, 116, 106, 91, 73, 55, 39, 0, 0, 0,
      35, 50, 66, 82, 96,  106, 109, 106, 96,  82, 66, 50, 35, 0, 0, 0,
      30, 43, 57, 70, 82,  91,  93,  91,  82,  70, 57, 43, 30, 0, 0, 0,
      24, 34, 46, 57, 66,  73,  75,  73,  66,  57, 46, 34, 24, 0, 0, 0,
      18, 26, 34, 43, 50,  55,  57,  55,  50,  43, 34, 26, 18, 0, 0, 0,
      13, 18, 24, 30, 35,  39,  40,  39,  35,  30, 24, 18, 13, 0, 0, 0, };

  unsigned char *tempOutput =
      (unsigned char *)sift_aligned_malloc(w * h * sizeof(unsigned char), 32);
  memcpy(tempOutput, input, w * h * sizeof(unsigned char));

  for (j = 6; j < h - 6; j++) {
    for (i = 6; i < w - 9; i++) {
      int r, new_pixel;
      __m128i sum32;
      sum32 = _mm_setzero_si128();
      for (r = 0; r < 13; r++) {
        const unsigned char *const pixel_ptr =
            input + (j + r - 6) * w + (i - 6);
        __m128i pixel8, pixel16_l, pixel16_h;
        __m128i coeff16_l, coeff16_h;

        if (((uintptr_t)(pixel_ptr)) & 0xF)
          pixel8 = _mm_loadu_si128((__m128i const *)pixel_ptr);
        else
          pixel8 = _mm_load_si128((__m128i const *)pixel_ptr);

        pixel16_l = _mm_cvtepu8_epi16(pixel8);
        pixel16_h = _mm_cvtepu8_epi16(_mm_srli_si128(pixel8, 8));  // >>8bytes

        coeff16_l = _mm_load_si128((__m128i const *)(kernel + r * 16));
        coeff16_h = _mm_load_si128((__m128i const *)(kernel + r * 16 + 8));

        sum32 = _mm_add_epi32(sum32, _mm_madd_epi16(pixel16_l, coeff16_l));
        sum32 = _mm_add_epi32(sum32, _mm_madd_epi16(pixel16_h, coeff16_h));
      }

      sum32 = _mm_add_epi32(sum32, _mm_srli_si128(sum32, 8));  // >>8bytes
      sum32 = _mm_add_epi32(sum32, _mm_srli_si128(sum32, 4));  // >>4bytes
      new_pixel = _mm_cvtsi128_si32(sum32);
      new_pixel = new_pixel / 10000;  // round, shift, not very high performance
      // new_pixel = new_pixel >> shift_bits;
      new_pixel = (((int)input[j * w + i]) << 1) - new_pixel;
      new_pixel = new_pixel > 255 ? 255 : new_pixel;
      new_pixel = new_pixel < 0 ? 0 : new_pixel;
      tempOutput[j * w + i] = (unsigned char)new_pixel;
    }

    for (; i < w - 6; i++) {
      int m, n, new_pixel;
      int sum = 0;
      for (m = 0; m < 13; m++)
        for (n = 0; n < 13; n++)
          sum += ((int)kernel[m * 16 + n]) * input[(j + m - 6) * w + i + n - 6];
      // sum = sum >> shift_bits;
      new_pixel = sum / 10000;
      new_pixel = (((int)input[j * w + i]) << 1) - new_pixel;
      new_pixel = new_pixel > 255 ? 255 : new_pixel;
      new_pixel = new_pixel < 0 ? 0 : new_pixel;
      tempOutput[j * w + i] = (unsigned char)new_pixel;
    }
  }

  memcpy(output, tempOutput, w * h * sizeof(unsigned char));
  sift_aligned_free((void **)(&tempOutput));
}

SIFT_EXPORT void add_weight_sse4_1(unsigned char *in1, const int width,
                                   const int stride1, const int height,
                                   const int factor1, unsigned char *in2out,
                                   const int stride2, const int factor2) {
  int x, y;
  unsigned char tmpval1;
  unsigned char tmpval2;
  unsigned char *ptrY1, *ptrY2;

  if (factor1 == 1 && factor2 == 1) {
    for (y = 0; y < height; y++) {
      ptrY1 = in1 + stride1 * y;
      ptrY2 = in2out + stride2 * y;

      __m128i pixel1_u8, pixel2_u8;
      for (x = 0; x <= width - 16; x += 16) {
        if (((uintptr_t)(ptrY1)) & 0xF)
          pixel1_u8 = _mm_loadu_si128((__m128i const *)ptrY1);
        else
          pixel1_u8 = _mm_load_si128((__m128i const *)ptrY1);

        if (((uintptr_t)(ptrY2)) & 0xF)
          pixel2_u8 = _mm_loadu_si128((__m128i const *)ptrY2);
        else
          pixel2_u8 = _mm_load_si128((__m128i const *)ptrY2);

        // u8 -> u16 zero extend
        __m128i pixel1_u16_l, pixel1_u16_h;
        __m128i pixel2_u16_l, pixel2_u16_h;
        pixel1_u16_l = _mm_cvtepu8_epi16(pixel1_u8);
        pixel1_u16_h = _mm_cvtepu8_epi16(_mm_srli_si128(pixel1_u8, 8));

        pixel2_u16_l = _mm_cvtepu8_epi16(pixel2_u8);
        pixel2_u16_h = _mm_cvtepu8_epi16(_mm_srli_si128(pixel2_u8, 8));

        // filter
        pixel2_u16_l = _mm_add_epi16(_mm_srli_epi16(pixel1_u16_l, 1),
                                     _mm_srli_epi16(pixel2_u16_l, 1));
        pixel2_u16_h = _mm_add_epi16(_mm_srli_epi16(pixel1_u16_h, 1),
                                     _mm_srli_epi16(pixel2_u16_h, 1));

        // xx-07-xx-06-xx-05-xx-04-xx-03-xx-02-xx-01-xx-00
        // xx-15-xx-14-xx-13-xx-12-xx-11-xx-10-xx-09-xx-08

        // shift
        // xx-07-xx-06-xx-05-xx-04-xx-03-xx-02-xx-01-xx-00
        // 15-xx-14-xx-13-xx-12-xx-11-xx-10-xx-09-xx-08-xx

        // or
        // 15-07-14-06-13-05-12-04-11-03-10-02-09-01-08-00
        pixel2_u16_l =
            _mm_or_si128(pixel2_u16_l, _mm_slli_epi16(pixel2_u16_h, 8));

        // shuffle
        // 15-07-14-06-13-05-12-04-11-03-10-02-09-01-08-00
        // 15-14-13-12-11-10-09-08-07-06-05-04-03-02-01-00
        // 15-13-11-09-07-05-03-01-14-12-10-08-06-04-02-00 (pattern)

        pixel2_u16_l = _mm_shuffle_epi8(
            pixel2_u16_l,
            _mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0));

        // store
        if (!(((uintptr_t)(ptrY2)) & 0xF))
          _mm_store_si128((__m128i *)ptrY2, pixel2_u16_l);
        else
          _mm_storeu_si128((__m128i *)ptrY2, pixel2_u16_l);

        ptrY1 += 16;
        ptrY2 += 16;
      }

      for (; x < width; x++) {
        tmpval1 = (*ptrY1) >> 1;
        tmpval2 = (*ptrY2) >> 1;

        *ptrY2 = tmpval1 + tmpval2;

        ptrY1++;
        ptrY2++;
      }
    }  // for(y=0; y< height; y++)
  }    // if(factor1==1 && factor2==1)
}
