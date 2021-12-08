/*
 * sift_gradient_calculation.h
 *
 *  Created on: Nov 3, 2014
 *      Author: peter
 */

#ifndef SIFT_GRADIENT_CALCULATION_H_
#define SIFT_GRADIENT_CALCULATION_H_

#include "sift_typedef.h"

void sift_update_gradient_f(SIFT_FILTER *f);

void sift_update_gradient_f_sse2(SIFT_FILTER *f);

void sift_update_gradient_f_avx(SIFT_FILTER *f);

#endif /* SIFT_GRADIENT_CALCULATION_H_ */
