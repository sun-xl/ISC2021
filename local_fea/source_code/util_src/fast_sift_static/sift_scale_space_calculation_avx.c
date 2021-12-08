/*
 * sift_scale_space_calculation_avx.c
 *
 *  Created on: Nov 10, 2014
 *      Author: peter
 */

#include "sift_scale_space_calculation.h"
#include "sift_math.h"
#include "sift_avx_typedef.h"
#include "sift_util.h"
#include <string.h>

// to do
static void copy_and_downsample_avx(float *dst, float const *src,
                                    const int width, const int height, int d) {
  int x, y;
  d = 1 << d;  // d = 2^d
  for (y = 0; y < height; y += d) {
    float const *srcrowp = src + y * width;
    for (x = 0; x < width - (d - 1); x += d) {
      *dst++ = *srcrowp;
      srcrowp += d;
    }
  }
}

static void copy_and_upsample_rows_avx(float *dst, float const *src,
                                       const int width, const int height) {
  int x, y;
  float a, b;

  UV256F32 a256, interp256;
  __m256 b256;
  const __m256 half_f32 = _mm256_set1_ps(0.5f);

  for (y = 0; y < height; y++) {
    for (x = 0; x <= width - 9; x += 8)  // avx 256
    {
      if (((uintptr_t)(src)) & 0x1F)
        a256.v = _mm256_loadu_ps(src);
      else
        a256.v = _mm256_load_ps(src);

      src++;
      if (((uintptr_t)(src)) & 0x1F)
        b256 = _mm256_loadu_ps(src);
      else
        b256 = _mm256_load_ps(src);
      src += 7;

      interp256.v = _mm256_mul_ps(
          half_f32, _mm256_add_ps(a256.v, b256));  // linear interpolation

#define copy_and_upsample_rows_avx_store(i) \
  *dst = a256.f32[i];                       \
  dst += height;                            \
  *dst = interp256.f32[i];                  \
  dst += height;
      copy_and_upsample_rows_avx_store(0)
      copy_and_upsample_rows_avx_store(1)
      copy_and_upsample_rows_avx_store(2)
      copy_and_upsample_rows_avx_store(3)
      copy_and_upsample_rows_avx_store(4)
      copy_and_upsample_rows_avx_store(5)
      copy_and_upsample_rows_avx_store(6)
      copy_and_upsample_rows_avx_store(7)
#undef copy_and_upsample_rows_avx_store
    }

    // tailing
    b = a = *src++;  // in case x=width-1, so also set b = a
    for (; x < width - 1; x++) {
      b = *src++;
      *dst = a;
      dst += height;
      *dst = 0.5 * (a + b);
      dst += height;
      a = b;
    }
    // last pixel
    *dst = b;
    dst += height;
    *dst = b;
    dst += height;
    dst += 1 - width * 2 * height;  // next column
  }
}

static void sift_covolution_column_f_avx(
    float *dst, const int dst_stride, float const *src, const int src_width,
    const int src_height, const int src_stride, float const *filt,
    const int filt_begin, const int filt_end, int step, unsigned int flags) {
  int x, y;
  const int dst_height = (src_height - 1) / step + 1;
  const int flag_transp = flags & FILTER_TRANSPOSE_OUT;
  const int flag_zeropad = (flags & FILTER_PAD_MASK) == FILTER_PAD_BY_ZERO;

  // let filt point to the last coefficient
  filt +=
      filt_end - filt_begin;  // so filt_begin will always start from index 0!!!
  for (x = 0; x < src_width;)  // along column
  {
    // calculate dst[x,y] = sum_p (src[x,p]*filt[y-p])
    // where supp(filt) = [filt_begin, filt_end] = [fb, fe]

    // chunk a: y-fe <= p < 0
    //			completes MAX(fe-y, 0) samples
    // chunk b: MAX(y-fe, 0) <= p < MIN(y-fb, height-1)
    //			completes fe - MAX(fb, height-y)+1 samples
    // chunk c: completes all samples
    // ______
    // |_fe_|
    // |    |
    // |    |
    // |____|
    // |_fb_|

    if (x + 8 < src_width)  // avx 256
    {
      for (y = 0; y < src_height; y += step) {
        int stop = filt_end - y;
        float const *srci = src + x - stop * src_stride;
        float const *filti = filt;

        UV256F32 acc256;
        //				__m256 v256, c256; // value, coefficient
        __m256 v256 = _mm256_setzero_ps(), c256;  // value, coefficient

        acc256.v = _mm256_setzero_ps();  // accumulation, set to 0s
        // chunk a
        if (stop > 0) {
          if (flag_zeropad)
            v256 = _mm256_setzero_ps();  // zero padding, set to 0s
          else {
            float const *srcx = src + x;  // continuity padding
            if (((uintptr_t)(srcx)) & 0x1F)
              v256 = _mm256_loadu_ps(srcx);
            else
              v256 = _mm256_load_ps(srcx);
          }
          while (filti > filt - stop) {
            c256 = _mm256_broadcast_ss(
                filti--);  // load 1 coefficient, copying to others elements
            acc256.v = _mm256_add_ps(acc256.v, _mm256_mul_ps(v256, c256));
            srci += src_stride;
          }
        }

        // chunk b
        stop = filt_end - SIFT_MAX(filt_begin, y - src_height + 1) + 1;
        while (filti > filt - stop) {
          // value
          if (((uintptr_t)(srci)) & 0x1F)
            v256 = _mm256_loadu_ps(srci);
          else
            v256 = _mm256_load_ps(srci);
          c256 = _mm256_broadcast_ss(
              filti--);  // load 1 coefficient, copying to others elements
          acc256.v = _mm256_add_ps(acc256.v, _mm256_mul_ps(v256, c256));
          srci += src_stride;
        }

        // chunk c
        if (flag_zeropad)              // else v = last pixel in current column
          v256 = _mm256_setzero_ps();  // zero padding, set to 0s
        stop = filt_end - filt_begin + 1;
        while (filti > filt - stop) {
          c256 = _mm256_broadcast_ss(
              filti--);  // load 1 coefficient, copying to others elements
          acc256.v = _mm256_add_ps(acc256.v, _mm256_mul_ps(v256, c256));
        }

        if (flag_transp)  // transpose
        {
#define sift_covolution_column_f_avx_transpose_store(i) \
  *dst = acc256.f32[i];                                 \
  dst += dst_stride;
          sift_covolution_column_f_avx_transpose_store(0)
          sift_covolution_column_f_avx_transpose_store(1)
          sift_covolution_column_f_avx_transpose_store(2)
          sift_covolution_column_f_avx_transpose_store(3)
          sift_covolution_column_f_avx_transpose_store(4)
          sift_covolution_column_f_avx_transpose_store(5)
          sift_covolution_column_f_avx_transpose_store(6)
          sift_covolution_column_f_avx_transpose_store(7)
#undef sift_covolution_column_f_avx_transpose_store
          dst += 1 - 8 * dst_stride;
        } else {
#define sift_covolution_column_f_avx_store(i) \
  *dst = acc256.f32[i];                       \
  dst++;
          sift_covolution_column_f_avx_store(0)
          sift_covolution_column_f_avx_store(1)
          sift_covolution_column_f_avx_store(2)
          sift_covolution_column_f_avx_store(3)
          sift_covolution_column_f_avx_store(4)
          sift_covolution_column_f_avx_store(5)
          sift_covolution_column_f_avx_store(6)
          sift_covolution_column_f_avx_store(7)
#undef sift_covolution_column_f_avx_store
          dst += dst_stride - 8;
        }
      }  // for (y=0; y<src_height; y+=step)

      // calculate dst_height pixels in y direction, goto next column
      if (flag_transp)  // transpose
        dst += 8 * dst_stride - dst_height;
      else
        dst += 8 - dst_height * dst_stride;
      x += 8;
    } else  // pure c
    {
      for (y = 0; y < src_height; y += step) {
        int stop = filt_end - y;
        float const *srci = src + x - stop * src_stride;
        float const *filti = filt;

        float acc = 0;   // accumulation
        float v = 0, c;  // value, coefficient

        // chunk a
        if (stop > 0) {
          if (flag_zeropad)
            v = 0;  // zero padding
          else
            v = *(src + x);  // continuity padding
          while (filti > filt - stop) {
            c = *filti--;  // coefficient
            acc += v * c;
            srci += src_stride;
          }
        }

        // chunk b
        stop = filt_end - SIFT_MAX(filt_begin, y - src_height + 1) + 1;
        while (filti > filt - stop) {
          v = *srci;     // value
          c = *filti--;  // coefficient
          acc += v * c;
          srci += src_stride;
        }

        // chunk c
        if (flag_zeropad)  // else v = last pixel in current column
          v = 0;
        stop = filt_end - filt_begin + 1;
        while (filti > filt - stop) {
          c = *filti--;  // coefficient
          acc += v * c;
        }

        *dst = acc;
        if (flag_transp)  // transpose
          dst += 1;
        else
          dst += dst_stride;
      }  // for (y=0; y<src_height; y+=step)

      // calculate dst_height pixels in y direction, goto next column
      if (flag_transp)  // transpose
        dst += dst_stride - dst_height;
      else
        dst += 1 - dst_height * dst_stride;
      x++;
    }  // if (x + 8 < src_width) // avx 256
  }    // for (x=0; x<src_width;)
}

// Gaussian smooth an image
static void sift_gaussian_smooth_f_avx(SIFT_FILTER *f, float *outputImage,
                                       float *tempImage,
                                       float const *inputImage, const int width,
                                       const int height, const float sigma) {
  // calculate Gaussian filter coefficients
  if (f->gaussFilterSigma != sigma) {
    int j;
    int gaussFilterSize;
    float gaussFilterAcc = 0;
    if (f->gaussFilter)  // free old coefficients space
      sift_aligned_free((void **)(&f->gaussFilter));
    // allocate new coefficients space
    f->gaussFilterSigma = sigma;
    f->gaussFilterWidth = SIFT_MAX(ceilf(4.0 * f->gaussFilterSigma), 1);
    gaussFilterSize = 2 * f->gaussFilterWidth + 1;
    f->gaussFilter =
        (float *)sift_aligned_malloc(gaussFilterSize * sizeof(float), 32);
    // calculate new coefficients 1d
    for (j = 0; j < gaussFilterSize; j++) {
      float d = ((float)(j - f->gaussFilterWidth)) / f->gaussFilterSigma;
      f->gaussFilter[j] = expf(-0.5 * (d * d));
      gaussFilterAcc += f->gaussFilter[j];
    }
    // normalize
    for (j = 0; j < gaussFilterSize; j++) f->gaussFilter[j] /= gaussFilterAcc;
  }

  if (f->gaussFilterWidth == 0)  // Peter: little possible
  {
    memcpy(outputImage, inputImage, sizeof(float) * width * height);
    return;
  }

  // convolution along column
  sift_covolution_column_f_avx(tempImage, height, inputImage, width, height,
                               width, f->gaussFilter, -f->gaussFilterWidth,
                               f->gaussFilterWidth, 1,
                               FILTER_PAD_BY_CONTINUITY | FILTER_TRANSPOSE_OUT);

  // convolution along column
  sift_covolution_column_f_avx(outputImage, width, tempImage, height, width,
                               height, f->gaussFilter, -f->gaussFilterWidth,
                               f->gaussFilterWidth, 1,
                               FILTER_PAD_BY_CONTINUITY | FILTER_TRANSPOSE_OUT);
}

int sift_process_first_octave_f_avx(SIFT_FILTER *f) {
  int o, s;
  float sa, sb, sd;
  float *octave_p, *octave_c;

  if (f->O == 0)  // need at least 1 octave
    return ERR_EOF;

  // compute the 1st level of the 1st octave
  // if the 1st octave has negative index, we up-scale the image;
  // if the 1st octave has positive index, we down-scale the image;
  // if the 1st octave has 0 index, we just copy the image
  if (f->o_cur < 0)  // up-scale
  {
    // f->o_cur == -1, double
    copy_and_upsample_rows_avx(f->temp, f->fimg, f->width,
                               f->height);  // w * h -> h * 2w
    copy_and_upsample_rows_avx(f->octave, f->temp, f->height,
                               f->width << 1);  // h * 2w -> 2w * 2h

    // f->o_cur < -1, double more
    for (o = -1; o > f->o_cur; --o) {
      const int w = f->width << -o;   // (2^-o)w
      const int h = f->height << -o;  // (2^-o)h

      copy_and_upsample_rows_avx(f->temp, f->octave, w, h);  // w * h -> h * 2w
      copy_and_upsample_rows_avx(f->octave, f->temp, h,
                                 w << 1);  // // h * 2w -> 2w * 2h
    }
  } else if (f->o_cur > 0)  // down-scale
    copy_and_downsample_avx(f->octave, f->fimg, f->width, f->height, f->o_cur);
  else
    // f->o_cur == 0, copy
    memcpy(f->octave, f->fimg, sizeof(float) * f->fimg_size);

  // we adjust the smoothing of the 1st level of the 1st octave.
  // the input image is assumed to have nominal smoothing equal to f->sigman
  sa = f->sigma0 *
       powf(f->sigmak, f->s_min);  // Peter: wrong??? assume f->o_min == 0???
  sb = f->sigman * powf(2.0, -f->o_min);
  if (sa > sb) {
    sd = sqrtf(sa * sa - sb * sb);
    sift_gaussian_smooth_f_avx(f, f->octave, f->temp, f->octave,
                               f->octave_width, f->octave_height, sd);
  }

  // continue calculating the 1st octave
  octave_p = f->octave;
  octave_c = octave_p + f->octave_size;
  for (s = f->s_min + 1; s <= f->s_max; s++) {
    sd = f->dsigma0 * powf(f->sigmak, s);
    sift_gaussian_smooth_f_avx(f, octave_c, f->temp, octave_p, f->octave_width,
                               f->octave_height, sd);
    octave_p = octave_c;
    octave_c += f->octave_size;
  }

  return ERR_OK;
}

int sift_process_next_octave_f_avx(SIFT_FILTER *f) {
  float sa, sb, sd;
  int s, s_base;
  float *octave_p, *octave_c;

  // already processed last octave
  if (f->o_cur == f->o_min + f->O - 1) return ERR_EOF;

  // retrieve base for the next octave
  s_base = SIFT_MIN(f->s_min + f->S,
                    f->s_max);  // Peter: is there something wrong???
  // s_base = MIN(f->s_min + f->S + 1, f->s_max);

  // downscale by 2, the base of the next octave
  copy_and_downsample_avx(
      f->octave,
      f->octave + f->octave_width * f->octave_height * (s_base - f->s_min),
      f->octave_width, f->octave_height, 1);

  // update f-> o_cur, update f-> octave_width, f-> octave_height due to
  // downscale
  f->o_cur += 1;
  f->xper = powf(2.0, f->o_cur);  // current octave
  f->octave_width = SIFT_SHIFT_LEFT(f->width, -f->o_cur);
  f->octave_height = SIFT_SHIFT_LEFT(f->height, -f->o_cur);
  f->octave_size = f->octave_width * f->octave_height;

  // 1st level of the current octave
  sa = f->sigma0 * powf(f->sigmak, f->s_min);
  sb = f->sigma0 * powf(f->sigmak, s_base - f->S);  // Peter: why sub f->S
  if (sa > sb) {
    sd = sqrtf(sa * sa - sb * sb);
    sift_gaussian_smooth_f_avx(f, f->octave, f->temp, f->octave,
                               f->octave_width, f->octave_height, sd);
  }

  // continue calculating the current octave
  octave_p = f->octave;
  octave_c = octave_p + f->octave_size;
  for (s = f->s_min + 1; s <= f->s_max; s++) {
    sd = f->dsigma0 * powf(f->sigmak, s);
    sift_gaussian_smooth_f_avx(f, octave_c, f->temp, octave_p, f->octave_width,
                               f->octave_height, sd);
    octave_p = octave_c;
    octave_c += f->octave_size;
  }

  return ERR_OK;
}
