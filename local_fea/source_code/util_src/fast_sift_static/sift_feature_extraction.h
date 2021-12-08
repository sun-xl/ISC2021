/*
 * sift_feature_extraction.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_FEATURE_EXTRACTION_H_
#define FAST_SIFT_SIFT_FEATURE_EXTRACTION_H_

#include "sift_typedef.h"

// const int octaves=5;
// const int levels=7;
// const int o_min=0;
// const double peak_thresh=0.5;
// const double edge_thresh=10.0;
// const double norm_thresh=-1;
// const double magnif=-1;
// const double window_size=-1;
// const int keypoint_min=2;
// const int keypoint_max=-1;
SIFT_EXPORT int new_sift_feature_f(
    SIFT_FEATURE *feature, const unsigned char *img, const int width,
    const int stride, const int height, const int octaves, const int level,
    const int o_min, const float peak_thresh, const float edge_thresh,
    const float norm_thresh, const float magnif, const float window_size,
    const int keypoint_min, const int keypoint_max);

SIFT_EXPORT int new_sift_feature_container_f(
    SIFT_FEATURE *feature, int initFeatNum, const int width, const int stride,
    const int height, const int octaves, const int levels, const int o_min,
    const float peak_thresh, const float edge_thresh, const float norm_thresh,
    const float magnif, const float window_size, const int keypoint_min,
    const int keypoint_max);

SIFT_EXPORT void delete_sift_feature_f(SIFT_FEATURE *feature);

// SIFT_EXPORT int sift_feature_extraction_f (SIFT_FEATURE *feature, SIFT_FILTER
// *g_filter);
SIFT_EXPORT int sift_feature_extraction_f(SIFT_FEATURE *feature);

#endif /* FAST_SIFT_SIFT_FEATURE_EXTRACTION_H_ */
