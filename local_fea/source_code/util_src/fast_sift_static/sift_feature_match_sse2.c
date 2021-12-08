/*
 * sift_feature_match.c
 *
 *  Created on: Oct 31, 2014
 *      Author: peter
 */

#include "sift_feature_match.h"
#include "sift_sse4_2_typedef.h"
#include "sift_util.h"
#include <math.h>

SIFT_EXPORT void match_keypoint_d_sse2(const int *L1_pt, const int *L2_pt,
                                       const double *pos_1, const double *pos_2,
                                       const int K1, const int K2,
                                       const float thresh, int *nMatches,
                                       int *Matches) {
  const int ND = 128;
  PAIR *pairs_begin = (PAIR *)calloc(sizeof(PAIR) * K1, 1);
  PAIR *pairs_end = pairs_begin;
  PAIR *search;

  int k1, k2;
  short *L1_pt_aligned =
      (short *)sift_aligned_malloc(sizeof(short) * ND * K1, 32);
  short *L2_pt_aligned =
      (short *)sift_aligned_malloc(sizeof(short) * ND * K2, 32);

  for (k1 = 0; k1 < ND * K1; k1++) L1_pt_aligned[k1] = L1_pt[k1];
  for (k2 = 0; k2 < ND * K2; k2++) L2_pt_aligned[k2] = L2_pt[k2];

  (*nMatches) = 0;  // set to 0
  for (k1 = 0; k1 < K1; k1++, L1_pt_aligned += ND) {
    int best = 0x7fffffff;
    int second_best = 0x7fffffff;
    int bestk = -1;  // initial

    for (k2 = 0; k2 < K2; k2++, L2_pt_aligned += ND) {
      int bin;
      union {
        __m128i acc128;
        int acc32[4];
      } accumulaton;

      accumulaton.acc128 = _mm_setzero_si128();

      for (bin = 0; bin < 128; bin += 8) {
        __m128i diff =
            _mm_sub_epi16(_mm_load_si128((__m128i *)(L1_pt_aligned + bin)),
                          _mm_load_si128((__m128i *)(L2_pt_aligned + bin)));
        accumulaton.acc128 =
            _mm_add_epi32(accumulaton.acc128, _mm_madd_epi16(diff, diff));
      }
      accumulaton.acc128 = _mm_add_epi32(accumulaton.acc128,
                                         _mm_srli_si128(accumulaton.acc128, 8));
      accumulaton.acc128 = _mm_add_epi32(accumulaton.acc128,
                                         _mm_srli_si128(accumulaton.acc128, 4));

      // best & 2nd best
      if (accumulaton.acc32[0] < best) {
        second_best = best;
        best = accumulaton.acc32[0];
        bestk = k2;
      } else if (accumulaton.acc32[0] < second_best) {
        second_best = accumulaton.acc32[0];
      }
    }  // for (j2=0; j2<K2; j2++, pt2_aligned+=ND)
    L2_pt_aligned -= ND * K2;  // reset the pointer

    // Lowe's method: accept the match only if unique
    if (thresh * (float)best < (float)second_best && bestk != -1) {
      int check = -1;
      for (search = pairs_begin; search < pairs_end; search++) {
        if ((int)rintf(pos_1[search->k1 * 4]) == (int)rintf(pos_1[k1 * 4]) &&
            (int)rintf(pos_1[1 + search->k1 * 4]) ==
                (int)rintf(pos_1[1 + k1 * 4])) {
          check = 1;
          if (best <= search->score) {
            search->k1 = k1;  // updating
            search->k2 = bestk;
            search->score = best;
          }
        }

        if ((int)rintf(pos_2[search->k2 * 4]) == (int)rintf(pos_2[bestk * 4]) &&
            (int)rintf(pos_2[1 + search->k2 * 4]) ==
                (int)rintf(pos_2[1 + bestk * 4])) {
          check = 1;
          if (best <= search->score) {
            search->k1 = k1;
            search->k2 = bestk;  // updating
            search->score = best;
          }
        }
      }

      if (check == -1) {
        pairs_end->k1 = k1;
        pairs_end->k2 = bestk;
        pairs_end->score = best;
        pairs_end++;
        (*nMatches)++;
      }
    }  //  if(thresh * (float)best < (float) second_best && bestk != -1)
  }    // for (j1=0; j1<K1; j1++, L1_pt+=ND)
  L1_pt_aligned -= ND * K1;  // reset the pointer

  int index = 0;
  for (search = pairs_begin; search < pairs_end; search++) {
    Matches[index++] = search->k1;
    Matches[index++] = search->k2;
  }

  free(pairs_begin);
  sift_aligned_free((void **)(&L1_pt_aligned));
  sift_aligned_free((void **)(&L2_pt_aligned));
  return;
}
