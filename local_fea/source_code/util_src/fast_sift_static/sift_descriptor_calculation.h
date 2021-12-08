/*
 * sift_descriptor_calculation.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_DESCRIPTOR_CALCULATION_H_
#define FAST_SIFT_SIFT_DESCRIPTOR_CALCULATION_H_

#include "sift_typedef.h"

void sift_calc_keypoint_descriptor_f(const SIFT_FILTER *f, float *descr,
                                     const SIFT_KEYPOINT *k,
                                     const float angle0);

void sift_calc_keypoint_descriptor_f_sse4_1(const SIFT_FILTER *f, float *descr,
                                            const SIFT_KEYPOINT *k,
                                            const float angle0);

void sift_calc_keypoint_descriptor_f_avx(const SIFT_FILTER *f, float *descr,
                                         const SIFT_KEYPOINT *k,
                                         const float angle0);

#endif /* FAST_SIFT_SIFT_DESCRIPTOR_CALCULATION_H_ */
