/*
 * sift_scale_space_calculation.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_SCALE_SPACE_CALCULATION_H_
#define FAST_SIFT_SIFT_SCALE_SPACE_CALCULATION_H_

#include "sift_typedef.h"

int sift_process_first_octave_f(SIFT_FILTER *f);
int sift_process_next_octave_f(SIFT_FILTER *f);

int sift_process_first_octave_f_sse(SIFT_FILTER *f);
int sift_process_next_octave_f_sse(SIFT_FILTER *f);

int sift_process_first_octave_f_avx(SIFT_FILTER *f);
int sift_process_next_octave_f_avx(SIFT_FILTER *f);

#endif /* FAST_SIFT_SIFT_SCALE_SPACE_CALCULATION_H_ */
