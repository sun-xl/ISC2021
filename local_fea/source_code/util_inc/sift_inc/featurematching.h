#ifndef FEATUREMATCHING_H
#define FEATUREMATCHING_H

#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <assert.h>
typedef unsigned char uchar;
typedef struct
{
  int k1 ;
  int k2 ;
  int score ;
} Pair ;

typedef struct
{
	int img_width;
	int img_height;
	int img_linesize;
	int octaves, levels, o_min;
	double edge_thresh, peak_thresh;
	int featureLimit;
	int ftureN;
//	int iniFT;
	int* fture;
	double* position;
	int frameInd;
}Feature;

typedef struct
{
	int matchN;
	int* matchPair;
}Match;

void initMatch(Match* match,int size);;
void resizeMatch(Match* match,int size);
void releaseMatch(Match* match);

#endif
