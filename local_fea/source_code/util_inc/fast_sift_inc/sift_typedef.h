/*
 * sift_typedef.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_TYPEDEF_H_
#define FAST_SIFT_SIFT_TYPEDEF_H_

#ifdef __cplusplus
#define SIFT_EXPORT extern "C"
#else
#define SIFT_EXPORT
#endif

#define OUTPUT_KEYPOINT_TYPE double

// error code
enum ERR_CODE {
  ERR_OK = 0,    // 0, no error
  ERR_NULL_PTR,  // 1, NULL pointer
  ERR_OVERFLOW,  // 2, buffer overflow
  ERR_ALLOC,     // 3, resource allocation error
  ERR_BAD_ARG,   // 4, bad argument or illegal data error
  ERR_IO,        // 5, input/output error
  ERR_EOF,       // 6, end-of-file or end-of-sequence error
  ERR_CODE_MAX   // 7
};

// Gaussian filter control flag
enum FILTER_FLAG {
  FILTER_PAD_BY_ZERO = (0x0 << 0),        // 0, pad with zeros
  FILTER_PAD_BY_CONTINUITY = (0x1 << 0),  // 1, pad by continuity
  FILTER_PAD_MASK = (0x3),                // 3, padding field selector
  FILTER_TRANSPOSE_OUT = (0x1 << 2),      // 4, transpose result
  FILTER_FLAG_MAX                         // 5
};

// optimization control flag
enum OPTIMIZE_FLAG {
  NO_OPTIMIZE = (0x00),           // 0
  OPTIMIZE_SSE = (0x01 << 0),     // 1
  OPTIMIZE_SSE2 = (0x01 << 1),    // 2
  OPTIMIZE_SSE3 = (0x01 << 2),    // 4
  OPTIMIZE_SSSE3 = (0x01 << 3),   // 8
  OPTIMIZE_SSE4_1 = (0x01 << 4),  // 16
  OPTIMIZE_SSE4_2 = (0x01 << 5),  // 32
  OPTIMIZE_AVX = (0x01 << 6),     // 64
  OPTIMIZE_AVX2 = (0x01 << 7),    // 128
  OPTIMIZE_MAX                    // 9
};

#define ORIENTATION_BINS 8  // 8 bins for the orientation dimension
#define SPATIAL_BINS 4      // 4 bins for each spatial dimension
#define SPATIAL_BINS2 (SPATIAL_BINS / 2)
#define KEYPOINT_SIZE 4  // 4
#define DESCRIPTOR_SIZE \
  (SPATIAL_BINS *SPATIAL_BINS *ORIENTATION_BINS)  // 4*4*8=128

typedef struct _SIFT_PARAMS {
  const unsigned char *img;
  int width;
  int stride;
  int height;

  int octaves;        // number of octaves
  int levels;         // number of levels per octave
  int o_min;          // minimum octave index
  float peak_thresh;  // contrast peak threshold
  float edge_thresh;  // edge threshold
  float norm_thresh;  // norm threshold
  float magnif;       // magnification factor
  float window_size;  // size of Gaussian window (in spatial bins)
  int keypoint_min;   // lower limit of keypoint number
  int keypoint_max;   // upper limit of keypoint number
} SIFT_PARAMS;

typedef struct _SIFT_FEATURE {
  SIFT_PARAMS params;
  int keypoint_res;                    // maximum keypoint number
  int keypoint_num;                    // current keypoint number
  OUTPUT_KEYPOINT_TYPE *keypoint_ptr;  // x, y, scale, orientation
  int *descriptor_ptr;                 // 128 dimensional descriptor
} SIFT_FEATURE;

// SIFT keypoint
typedef struct _SIFT_KEYPOINT {
  int o;  // o coordinate (octave)

  int ix;  // integer unnormalized x coordinate
  int iy;  // integer unnormalized y coordinate
  int is;  // integer s coordinate

  float x;      // x coordinate
  float y;      // y coordinate
  float s;      // s coordinate
  float sigma;  // scale
} SIFT_KEYPOINT;

// SIFT filter, this filter implements the SIFT detector and descriptor
typedef struct _SIFT_FILTER {
  float *fimg;    // image data, no stride
  int width;      // image width
  int height;     // image height
  int fimg_size;  // image size

  float sigman;   // nominal image smoothing
  float sigmak;   // k-smoothing
  float sigma0;   // smoothing of pyramid base
  float dsigma0;  // delta-smoothing

  int O;      // number of octaves
  int S;      // number of levels per octave
  int o_min;  // minimum octave index
  int s_min;  // minimum level index
  int s_max;  // maximum level index
  int o_cur;  // current octave
  float xper;

  float *temp;        // temporary pixel buffer
  float *octave;      // current GSS data
  float *dog;         // current DoG data
  int octave_width;   // current octave width, no stride
  int octave_height;  // current octave height
  int octave_size;

  float *gaussFilter;      // current Gaussian filter
  float gaussFilterSigma;  // current Gaussian filter std
  int gaussFilterWidth;    // current Gaussian filter width

  SIFT_KEYPOINT *keys;  // detected keypoints
  int nkeys;            // number of detected keypoints
  int keys_res;         // size of the keys buffer

  float peak_thresh;  // contrast peak threshold
  float edge_thresh;  // edge threshold
  float norm_thresh;  // norm threshold
  float magnif;       // magnification factor
  float windowSize;   // size of Gaussian window (in spatial bins)
  int keypoint_min;   // lower limit of keypoint number
  int keypoint_max;   // upper limit of keypoint number

  float *grad;  // GSS gradient data
  int grad_o;   // GSS gradient data octave

  int optimize_type;  // optimization type

  float *fast_expn_tab;  // fast_expn table, 0~25.0
} SIFT_FILTER;

#endif /* FAST_SIFT_SIFT_TYPEDEF_H_ */
