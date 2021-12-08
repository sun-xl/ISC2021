/*
 * sift_prefilter.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_PREFILTER_H_
#define FAST_SIFT_SIFT_PREFILTER_H_

#include "sift_typedef.h"

SIFT_EXPORT void sharp_filter_sse4_1(const unsigned char* const input,
                                     unsigned char* const output, const int w,
                                     const int h);

SIFT_EXPORT void add_weight_sse4_1(unsigned char* in1, const int width,
                                   const int stride1, const int height,
                                   const int factor1, unsigned char* in2out,
                                   const int stride2, const int factor2);

#endif /* FAST_SIFT_SIFT_PREFILTER_H_ */
