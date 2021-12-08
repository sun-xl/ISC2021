/*
 * sift_feature_match.h
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#ifndef FAST_SIFT_SIFT_FEATURE_MATCH_H_
#define FAST_SIFT_SIFT_FEATURE_MATCH_H_

#include "sift_typedef.h"

typedef struct _PAIR {
  int k1;
  int k2;
  int score;
} PAIR;

SIFT_EXPORT void match_keypoint_d_sse2(const int *L1_pt, const int *L2_pt,
                                       const double *pos_1, const double *pos_2,
                                       const int K1, const int K2,
                                       const float thresh, int *nMatches,
                                       int *Matches);

#endif /* FAST_SIFT_SIFT_FEATURE_MATCH_H_ */
