#include <stdio.h>
extern "C" {
#include<libavcodec/avcodec.h>
#include<libavutil/avutil.h>
#include<libavformat/avformat.h>
}
#include "sift_typedef.h"
#include "featurematching.h"
#ifdef __MINGW32__
#undef WIN32
#endif

void imageFeatureExtractionOne(AVFrame *dec_frame, SIFT_FEATURE* feature, int numFeatute,int isRootSift=0);

int imgCompareTwo_xl(SIFT_FEATURE* feat1, SIFT_FEATURE* feat2,int isStrit=0,int flip=1);