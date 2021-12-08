/*
 * sift_orientation_assignment.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_ORIENTATION_ASSIGNMENT_H_
#define FAST_SIFT_SIFT_ORIENTATION_ASSIGNMENT_H_

#include "sift_typedef.h"

int sift_calc_keypoint_orientations_f(const SIFT_FILTER *const f,
                                      float angles[4],
                                      const SIFT_KEYPOINT *const k);

#endif /* FAST_SIFT_SIFT_ORIENTATION_ASSIGNMENT_H_ */
