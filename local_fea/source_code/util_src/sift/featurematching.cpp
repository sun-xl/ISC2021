#include "featurematching.h"
#include <iostream>
#define verbose 0

void initMatch(Match* match,int size){
	match->matchN=0;
	match->matchPair = new int[2*size];
	//match->matchPair = (int*)calloc(2*size,sizeof(int));
}

void resizeMatch(Match* match,int size){
	//match->matchPair = (int*)realloc(match->matchPair,2*sizeof(int)*size); // = Y X Scale Angle
	delete [] match->matchPair;
	match->matchPair = new int[2*size];
}

void releaseMatch(Match* match){
	delete [] match->matchPair;
	//free(match->matchPair);
}

