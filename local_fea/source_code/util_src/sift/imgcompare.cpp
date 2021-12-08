#define __STDC_CONSTANT_MACROS

#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <algorithm>
#include <ctime>
#include <sstream>
#include "sift_feature_match.h"
#include "sift_prefilter.h"
#include "sift_feature_extraction.h"
#include "imgcompare.h"
extern "C"
{
	#include "util.h"
	#include "avf_conf.h"
}

using namespace std;

pthread_mutex_t proc_lock = PTHREAD_MUTEX_INITIALIZER;

void imageFeatureExtractionOne(AVFrame *dec_frame, SIFT_FEATURE* feature, int numFeatute,int isRootSift)
{
		double thres=10.0f;
		pthread_mutex_t write_locker = PTHREAD_MUTEX_INITIALIZER;
		
		pthread_mutex_lock(&write_locker);
		new_sift_feature_f (feature,
				dec_frame->data[0], dec_frame->width, dec_frame->linesize[0], dec_frame->height,5, 7, 0, 0.5f, thres, -1.0f, -1.0f, -1.0f, 2, numFeatute);
		pthread_mutex_unlock(&write_locker);
			      //5, 7, 0, 0.5f, thres, -1.0f, -1.0f, -1.0f, 2, numFeatute);
		sift_feature_extraction_f (feature);
		
		if(feature->params.keypoint_max > 0)
		{
			while(feature->keypoint_num > feature->params.keypoint_max){ //if it is more than limit
				thres -= (0.1 * feature->keypoint_num/numFeatute);
				pthread_mutex_lock(&write_locker);
				delete_sift_feature_f(feature);
				new_sift_feature_f (feature,
						    		dec_frame->data[0], 
									dec_frame->width, dec_frame->linesize[0], dec_frame->height,
									5, 7, 0, 0.4f, thres, -1.0f, -1.0f, -1.0f, 2, numFeatute);
				pthread_mutex_unlock(&write_locker);
				sift_feature_extraction_f (feature);
//printf("%d e%lf p%lf\n", feature->keypoint_num, feature->params.edge_thresh, feature->params.peak_thresh);
				if( (feature->keypoint_num > feature->params.keypoint_max*0.8) 
				  && (feature->keypoint_num < feature->params.keypoint_max) ) break;
			}
		}
////change to root sift
		if(isRootSift==1){
			float tmp;
			float feattmp;
			for(int i=0;i<feature->keypoint_num;i++){
				tmp=0;
				for(int j=0;j<128;j++){
					tmp+=abs(feature->descriptor_ptr[i*128+j]);
				}
				if(tmp<0.001)
					continue;
				for(int j=0;j<128;j++){
					feattmp = sqrt(fabs(1.0*feature->descriptor_ptr[i*128+j]/tmp))*255;
					feature->descriptor_ptr[i*128+j]= int(feattmp+0.5);
				}
			}
		}
}


int imageFeatureMatchingTwo(SIFT_FEATURE* feature1, SIFT_FEATURE* feature2,int strict)
{
	double hMatrix[9];
	Match match;
	initMatch(&match,1);
	int matchKeyP=0;
	resizeMatch(&match,feature1->keypoint_num);

	match.matchN=0;

	match_keypoint_d_sse2 (feature1->descriptor_ptr,feature2->descriptor_ptr,
					feature1->keypoint_ptr,feature2->keypoint_ptr,
					feature1->keypoint_num, feature2->keypoint_num,
					1.8, &match.matchN, match.matchPair);

	matchKeyP = match.matchN;
	releaseMatch(&match);
	return matchKeyP;
        
}



int imgCompareTwo_xl(SIFT_FEATURE* feat1, SIFT_FEATURE* feat2,int isStrit,int flip){ // sunxl
	
	int matchKeyP = imageFeatureMatchingTwo(feat1,feat2,isStrit);
	
	if(flip==1){

		int binorder[16] = {12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3};
		int bininorder[8] = {0,7,6,5,4,3,2,1};
		float tmpdata[128];	
		int indexsrc,indexref;	
		for(int i=0;i<feat1->keypoint_num;i++){
			//memcpy (tmpdata,feat1->descriptor_ptr+i*128, sizeof(float)*128 );
			for(int j=0;j<128;j++){
				tmpdata[j] = feat1->descriptor_ptr[i*128+j];
			}
			for(int n=0;n<16;n++){
				for(int m=0;m<8;m++){
					indexsrc = n*8+m;
					indexref = binorder[n]*8+bininorder[m];
//cout<<"src: "<<n<<","<<m<<":"<<indexsrc<<"\tref: "<<binorder[n]<<","<<bininorder[m]<<": "<<indexref<<endl;
					feat1->descriptor_ptr[i*128+indexsrc] = tmpdata[indexref];
				}
			}
		}
		int matchKeyP2 = imageFeatureMatchingTwo(feat1, feat2,isStrit);
		
		matchKeyP = std::max(matchKeyP, matchKeyP2);
	}	
	//pthread_mutex_unlock(&proc_lock);
	
	return matchKeyP;
}
