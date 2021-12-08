/*
 * sift_feature_extraction.c
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#include "sift_feature_extraction.h"

#include "sift_util.h"
#include "sift_math.h"
#include "sift_scale_space_calculation.h"
#include "sift_keypoint_detection.h"
#include "sift_gradient_calculation.h"
#include "sift_orientation_assignment.h"
#include "sift_descriptor_calculation.h"

#include <stdio.h>
#include <stdlib.h>
// SIFT_FILTER g_filter; // sift filter
// static SIFT_FILTER g_filter; // sift filter

// buffer using aligned malloc & free
// SIFT_FEATURE::keypoint_ptr
// SIFT_FEATURE::descriptor_ptr
SIFT_EXPORT int new_sift_feature_f(
    SIFT_FEATURE *feature, const unsigned char *img, const int width,
    const int stride, const int height, const int octaves, const int levels,
    const int o_min, const float peak_thresh, const float edge_thresh,
    const float norm_thresh, const float magnif, const float window_size,
    const int keypoint_min, const int keypoint_max) {
  if (!feature || !img) return ERR_NULL_PTR;
  if (width <= 0 || stride <= 0 || height <= 0 ||
      levels <= 0)  // if levels==0, divide zero error
    return ERR_BAD_ARG;

  // copy parameters to feature structure
  feature->params.img = img;
  feature->params.width = width;
  feature->params.stride = stride;
  feature->params.height = height;
  feature->params.octaves = octaves;
  feature->params.levels = levels;
  feature->params.o_min = o_min;
  feature->params.peak_thresh = peak_thresh;
  feature->params.edge_thresh = edge_thresh;
  feature->params.norm_thresh = norm_thresh;
  feature->params.magnif = magnif;
  feature->params.window_size = window_size;
  feature->params.keypoint_min = keypoint_min;
  feature->params.keypoint_max = keypoint_max;

  feature->keypoint_res = feature->params.width * feature->params.height >>
                          2;  // initial capacity, (width*height)/4
  feature->keypoint_num = 0;  // initial keypoint number

  feature->keypoint_ptr = (OUTPUT_KEYPOINT_TYPE *)sift_aligned_malloc(
      KEYPOINT_SIZE * feature->keypoint_res * sizeof(OUTPUT_KEYPOINT_TYPE), 32);
  feature->descriptor_ptr = (int *)sift_aligned_malloc(
      DESCRIPTOR_SIZE * feature->keypoint_res * sizeof(int), 32);
  if (feature->keypoint_ptr == NULL || feature->descriptor_ptr == NULL) {
    printf("KEYPOINT_SIZE: %d\n", KEYPOINT_SIZE);
    printf("keypoint_res: %d\n", feature->keypoint_res);
    printf("DESCRIPTOR_SIZE: %d\n", DESCRIPTOR_SIZE);
  }
  return ERR_OK;
}

SIFT_EXPORT int new_sift_feature_container_f(
    SIFT_FEATURE *feature, int initFeatNum, const int width, const int stride,
    const int height, const int octaves, const int levels, const int o_min,
    const float peak_thresh, const float edge_thresh, const float norm_thresh,
    const float magnif, const float window_size, const int keypoint_min,
    const int keypoint_max) {
  if (!feature) return ERR_NULL_PTR;
  if (width <= 0 || stride <= 0 || height <= 0 ||
      levels <= 0)  // if levels==0, divide zero error
    return ERR_BAD_ARG;

  feature->params.img = NULL;
  feature->params.width = width;
  feature->params.stride = stride;
  feature->params.height = height;
  feature->params.octaves = octaves;
  feature->params.levels = levels;
  feature->params.o_min = o_min;
  feature->params.peak_thresh = peak_thresh;
  feature->params.edge_thresh = edge_thresh;
  feature->params.norm_thresh = norm_thresh;
  feature->params.magnif = magnif;
  feature->params.window_size = window_size;
  feature->params.keypoint_min = keypoint_min;
  feature->params.keypoint_max = keypoint_max;

  feature->keypoint_res = initFeatNum;  // initial capacity, (width*height)/4
  feature->keypoint_num = initFeatNum;  // initial keypoint number
  feature->keypoint_ptr = (OUTPUT_KEYPOINT_TYPE *)sift_aligned_malloc(
      KEYPOINT_SIZE * feature->keypoint_res * sizeof(OUTPUT_KEYPOINT_TYPE), 32);
  feature->descriptor_ptr = (int *)sift_aligned_malloc(
      DESCRIPTOR_SIZE * feature->keypoint_res * sizeof(int), 32);
  if (feature->keypoint_ptr == NULL || feature->descriptor_ptr == NULL) {
    printf("KEYPOINT_SIZE: %d\n", KEYPOINT_SIZE);
    printf("keypoint_res: %d\n", feature->keypoint_res);
    printf("DESCRIPTOR_SIZE: %d\n", DESCRIPTOR_SIZE);
  }
  return ERR_OK;
}

SIFT_EXPORT void delete_sift_feature_f(SIFT_FEATURE *feature) {
  if (feature) {
    feature->params.img = NULL;
    feature->keypoint_res = 0;
    feature->keypoint_num = 0;
    if (feature->keypoint_ptr != NULL)
      sift_aligned_free((void **)(&feature->keypoint_ptr));
    if (feature->descriptor_ptr != NULL)
      sift_aligned_free((void **)(&feature->descriptor_ptr));

    // free (feature);
    // feature = NULL;
  }
}

// buffer using aligned malloc & free
// SIFT_FILTER::fimg
// SIFT_FILTER::temp
// SIFT_FILTER::octave
// SIFT_FILTER::dog
// SIFT_FILTER::gaussFilter
// SIFT_FILTER::grad
int init_sift_filter_f(SIFT_FEATURE *feature, SIFT_FILTER *g_filter) {
  int i, j;

  if (!feature) return ERR_NULL_PTR;

  if (g_filter->fimg)  // image data
    sift_aligned_free((void **)(g_filter->fimg));

  // transpose
  g_filter->width = feature->params.height;
  g_filter->height = feature->params.width;
  g_filter->fimg_size = g_filter->width * g_filter->height;
  g_filter->fimg =
      (float *)sift_aligned_malloc(g_filter->fimg_size * sizeof(float), 32);

  for (i = 0; i < feature->params.height; i++)
    for (j = 0; j < feature->params.width; j++)
      g_filter->fimg[j * g_filter->width + i] =
          (float)feature->params.img[i * feature->params.stride + j];

  g_filter->sigman = 0.5f;  // nominal image smoothing
  g_filter->sigmak = powf(2.0, 1.0 / feature->params.levels);  // k-smoothing
  g_filter->sigma0 = 1.6 * g_filter->sigmak;  // smoothing of pyramid base
  g_filter->dsigma0 = g_filter->sigma0 *
                      sqrtf(1.0 - 1.0 / (g_filter->sigmak *
                                         g_filter->sigmak));  // delta-smoothing

  g_filter->O = feature->params.octaves;  // number of octaves
  if (g_filter->O < 0)  // in case the number of the octave is invalid
  {
    g_filter->O =
        SIFT_MAX(floorf(log2f(SIFT_MIN(g_filter->width, g_filter->height))) -
                     g_filter->o_min - 3,
                 1);
    printf("init_sift_filter_f: new octave number is %d\n", g_filter->O);
  }
  g_filter->S = feature->params.levels;          // number of levels per octave
  g_filter->o_min = feature->params.o_min;       // minimum octave index
  g_filter->s_min = -1;                          // minimum level index
  g_filter->s_max = feature->params.levels + 1;  // maximum level index
  g_filter->o_cur = g_filter->o_min;  // current octave is the 1st octave
  g_filter->xper = powf(2.0, g_filter->o_cur);  // current octave

  if (g_filter->temp)  // temporary buffer
    sift_aligned_free((void **)(g_filter->temp));
  if (g_filter->octave)  // current GSS
    sift_aligned_free((void **)(g_filter->octave));
  if (g_filter->dog)  // current DoG
    sift_aligned_free((void **)(g_filter->dog));
  // current is the 1st octave
  g_filter->octave_width = SIFT_SHIFT_LEFT(g_filter->width, -g_filter->o_cur);
  g_filter->octave_height = SIFT_SHIFT_LEFT(g_filter->height, -g_filter->o_cur);
  g_filter->octave_size = g_filter->octave_width * g_filter->octave_height;
  g_filter->temp =
      (float *)sift_aligned_malloc(g_filter->octave_size * sizeof(float), 32);
  g_filter->octave =
      (float *)sift_aligned_malloc(g_filter->octave_size * sizeof(float) *
                                       (g_filter->s_max - g_filter->s_min + 1),
                                   32);  // S+3
  g_filter->dog =
      (float *)sift_aligned_malloc(g_filter->octave_size * sizeof(float) *
                                       (g_filter->s_max - g_filter->s_min),
                                   32);  // S+2

  if (g_filter->gaussFilter)  // current Gaussian filter
    sift_aligned_free((void **)(g_filter->gaussFilter));
  g_filter->gaussFilterSigma =
      0.0f;                        // current Gaussian filter standard deviation
  g_filter->gaussFilterWidth = 0;  // current Gaussian filter width (radius)

  if (g_filter->keys)  // detected keypoints
  {
    free(g_filter->keys);
    g_filter->keys = NULL;
  }
  // allocate keypoints buffer
  g_filter->keys_res = 1000;
  g_filter->keys =
      (SIFT_KEYPOINT *)calloc(g_filter->keys_res * sizeof(SIFT_KEYPOINT), 1);
  g_filter->nkeys = 0;

  // default
  g_filter->peak_thresh = 0.0f;   // contrast peak threshold
  g_filter->edge_thresh = 10.0f;  // edge threshold
  g_filter->norm_thresh = 0.0f;   // L2 norm threshold
  g_filter->magnif = 3.0f;        // magnification factor
  g_filter->windowSize =
      SPATIAL_BINS / 2;  // 2, width (radius) of Gaussian window in spatial bins
  g_filter->keypoint_min = 2;   // minimum keypoints number
  g_filter->keypoint_max = -1;  // maximum keypoints number
  // user defined
  if (feature->params.peak_thresh >= 0)
    g_filter->peak_thresh = feature->params.peak_thresh;
  if (feature->params.edge_thresh >= 0)
    g_filter->edge_thresh = feature->params.edge_thresh;
  if (feature->params.norm_thresh >= 0)
    g_filter->norm_thresh = feature->params.norm_thresh;
  if (feature->params.magnif >= 0) g_filter->magnif = feature->params.magnif;
  if (feature->params.window_size >= 0)
    g_filter->windowSize = feature->params.window_size;
  if (feature->params.keypoint_min > 0)
    g_filter->keypoint_min = feature->params.keypoint_min;
  if (feature->params.keypoint_max > 0)
    g_filter->keypoint_max = feature->params.keypoint_max;

  if (g_filter->grad)  // current GSS gradient
    sift_aligned_free((void **)(g_filter->grad));
  g_filter->grad_o = g_filter->o_min - 1;  // current GSS gradient octave
  g_filter->grad =
      (float *)sift_aligned_malloc(g_filter->octave_size * sizeof(float) *
                                       (g_filter->s_max - g_filter->s_min) * 2,
                                   32);  // 2*(S+2)

  // sift optimization type
  g_filter->optimize_type = cpu_simd_instruction();

  // initialize fast_expn stuff
  g_filter->fast_expn_tab = fast_expn_init();

  return ERR_OK;
}

void uninit_sift_filter_f(SIFT_FILTER *g_filter) {
  if (g_filter->fimg)  // image data
    sift_aligned_free((void **)(&(g_filter->fimg)));
  g_filter->width = 0;   // image width
  g_filter->height = 0;  // image height
  g_filter->fimg_size = 0;

  g_filter->o_cur = g_filter->o_min;  // current octave is the 1st octave
  g_filter->xper = powf(2.0, g_filter->o_cur);  // current octave scale

  if (g_filter->temp)  // temporary buffer
    sift_aligned_free((void **)(&(g_filter->temp)));
  if (g_filter->octave)  // current GSS
    sift_aligned_free((void **)(&(g_filter->octave)));
  if (g_filter->dog)  // current DoG
    sift_aligned_free((void **)(&(g_filter->dog)));
  g_filter->octave_width = 0;   // current octave width
  g_filter->octave_height = 0;  // current octave height
  g_filter->octave_size = 0;

  if (g_filter->gaussFilter)  // current Gaussian filter
    sift_aligned_free((void **)(&(g_filter->gaussFilter)));
  g_filter->gaussFilterSigma =
      0.0f;                        // current Gaussian filter standard deviation
  g_filter->gaussFilterWidth = 0;  // current Gaussian filter width (radius)

  if (g_filter->keys)  // detected keypoints
  {
    free(g_filter->keys);
    g_filter->keys = NULL;
  }

  g_filter->nkeys = 0;     // number of detected keypoints
  g_filter->keys_res = 0;  // capacity of the keypoints buffer

  if (g_filter->grad)  // current GSS gradient
    sift_aligned_free((void **)(&(g_filter->grad)));
  g_filter->grad_o = g_filter->o_min - 1;  // current GSS gradient octave
}

// SIFT_EXPORT int sift_feature_extraction_f (SIFT_FEATURE *feature, SIFT_FILTER
// *g_filter)
SIFT_EXPORT int sift_feature_extraction_f(SIFT_FEATURE *feature) {
  int i, j, first;
  SIFT_FILTER g_filter;
  g_filter.fimg = NULL;
  g_filter.temp = NULL;
  g_filter.octave = NULL;
  g_filter.dog = NULL;
  g_filter.gaussFilter = NULL;
  g_filter.grad = NULL;
  g_filter.keys = NULL;

  if (!feature) return ERR_NULL_PTR;

  init_sift_filter_f(feature, &g_filter);

  // process each octave
  first = 1;
  while (1) {
    int err;
    const SIFT_KEYPOINT *keys;
    int nkeys = 0;

    if (first) {
      if (g_filter.optimize_type & OPTIMIZE_AVX)
        err = sift_process_first_octave_f_avx(&g_filter);
      else if (g_filter.optimize_type & OPTIMIZE_SSE)
        err = sift_process_first_octave_f_sse(&g_filter);
      else
        err = sift_process_first_octave_f(&g_filter);
      first = 0;
    } else {
      if (g_filter.optimize_type & OPTIMIZE_AVX)
        err = sift_process_next_octave_f_avx(&g_filter);
      else if (g_filter.optimize_type & OPTIMIZE_SSE)
        err = sift_process_next_octave_f_sse(&g_filter);
      else
        err = sift_process_next_octave_f(&g_filter);
    }

    if (err) break;

    if (g_filter.optimize_type & OPTIMIZE_AVX)
      sift_detect_f_avx(&g_filter);
    else if (g_filter.optimize_type & OPTIMIZE_SSE2)
      sift_detect_f_sse2(&g_filter);
    else
      sift_detect_f(&g_filter);

    if (g_filter.optimize_type & OPTIMIZE_AVX)
      sift_update_gradient_f_avx(&g_filter);
    else if (g_filter.optimize_type & OPTIMIZE_SSE2)
      sift_update_gradient_f_sse2(&g_filter);
    else
      sift_update_gradient_f(&g_filter);

    keys = g_filter.keys;
    nkeys = g_filter.nkeys;
    for (i = 0; i < nkeys; ++i) {
      float angles[4];
      int nangles;
      int q;
      const SIFT_KEYPOINT *k;
      k = keys + i;
      nangles = sift_calc_keypoint_orientations_f(&g_filter, angles,
                                                  k);  // compute orientation

      for (q = 0; q < nangles; ++q) {
        float descr[DESCRIPTOR_SIZE] __attribute__((aligned(32)));
        float rdescr[DESCRIPTOR_SIZE] __attribute__((aligned(32)));

        if (g_filter.optimize_type & OPTIMIZE_AVX)
          sift_calc_keypoint_descriptor_f_avx(&g_filter, descr, k, angles[q]);
        else if (g_filter.optimize_type & OPTIMIZE_SSE4_1)
          sift_calc_keypoint_descriptor_f_sse4_1(&g_filter, descr, k,
                                                 angles[q]);
        else  // pure c
          sift_calc_keypoint_descriptor_f(&g_filter, descr, k, angles[q]);

        // transpose
        for (j = 0; j < SPATIAL_BINS; j++) {
          int k, t;
          int jp = SPATIAL_BINS - 1 - j;
          for (k = 0; k < SPATIAL_BINS; k++) {
            int o = ORIENTATION_BINS * k + SPATIAL_BINS * ORIENTATION_BINS * j;
            int op =
                ORIENTATION_BINS * k + SPATIAL_BINS * ORIENTATION_BINS * jp;
            rdescr[op] = descr[o];  // vertical
            for (t = 1; t < ORIENTATION_BINS; t++)
              rdescr[ORIENTATION_BINS - t + op] =
                  descr[t + o];  // transpose orientation
          }
        }

        // make room for more keypoints and descriptor
        if (feature->keypoint_num >= feature->keypoint_res) {
          feature->keypoint_res += 1000;
          feature->keypoint_ptr = (OUTPUT_KEYPOINT_TYPE *)sift_aligned_realloc(
              (void **)(&feature->keypoint_ptr),
              KEYPOINT_SIZE * feature->keypoint_res *
                  sizeof(OUTPUT_KEYPOINT_TYPE),
              32);  // x, y, scale, orientation
          feature->descriptor_ptr = (int *)sift_aligned_realloc(
              (void **)(&feature->descriptor_ptr),
              DESCRIPTOR_SIZE * feature->keypoint_res * sizeof(int),
              32);  // 128 dimensional descriptor
        }

        // transpose
        feature->keypoint_ptr[KEYPOINT_SIZE * feature->keypoint_num + 0] = k->y;
        feature->keypoint_ptr[KEYPOINT_SIZE * feature->keypoint_num + 1] = k->x;
        feature->keypoint_ptr[KEYPOINT_SIZE * feature->keypoint_num + 2] =
            k->sigma;
        feature->keypoint_ptr[KEYPOINT_SIZE * feature->keypoint_num + 3] =
            PI / 2 - angles[q];

        for (j = 0; j < DESCRIPTOR_SIZE; ++j) {
          float x = 512.0F * rdescr[j];
          x = (x < 255.0F) ? x : 255.0F;
          feature->descriptor_ptr[DESCRIPTOR_SIZE * feature->keypoint_num + j] =
              (int)x;
        }
        feature->keypoint_num++;
      }  // for (q = 0 ; q < nangles ; ++q)
    }    // for (i = 0; i < nkeys; ++i)

    if (nkeys <= g_filter.keypoint_min) break;
    if (nkeys >= g_filter.keypoint_max &&
        g_filter.keypoint_max >
            0)  // stop if the current ovtave already have enough feature point
      break;
  }  // while (1)

  // trim
  feature->keypoint_ptr = (OUTPUT_KEYPOINT_TYPE *)sift_aligned_realloc(
      (void **)(&feature->keypoint_ptr),
      KEYPOINT_SIZE * feature->keypoint_num * sizeof(OUTPUT_KEYPOINT_TYPE),
      32);  // x, y, scale, orientation
  feature->descriptor_ptr = (int *)sift_aligned_realloc(
      (void **)(&feature->descriptor_ptr),
      DESCRIPTOR_SIZE * feature->keypoint_num * sizeof(int),
      32);  // 128 dimensional descriptor
  feature->keypoint_res = feature->keypoint_num;

  uninit_sift_filter_f(&g_filter);

  return ERR_OK;
}
