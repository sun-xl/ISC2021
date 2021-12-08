/*
 * sift_keypoint_detection.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_KEYPOINT_DETECTION_H_
#define FAST_SIFT_SIFT_KEYPOINT_DETECTION_H_

#include "sift_typedef.h"

void sift_detect_f(SIFT_FILTER* f);

void sift_detect_f_sse2(SIFT_FILTER* f);

void sift_detect_f_avx(SIFT_FILTER* f);

#endif /* FAST_SIFT_SIFT_KEYPOINT_DETECTION_H_ */
