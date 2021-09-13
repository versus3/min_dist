// vc.h
// Voronoi Cell Method
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 10/06/2010
//		Last Modified: 10/06/2010

#ifndef __VORONOI
#define __VORONOI

#include <queue>
#include <map>
#include "../point/point.h"
#include "../voronoi/voronoiCell.h"
#include "../simpleblockfile/sblockfile.h"

using namespace std;

typedef struct QfNode{
	int nBlock;
	float fMinDist;
}QfNode;

struct cmp{
    bool operator() (QfNode* pNode1, QfNode* pNode2){
		return pNode1->fMinDist > pNode2->fMinDist;
	}
};

float dd(float* pfdnn, PotentialLocation* p, VoronoiCell* pVC, RTNode* pNodeC, float* fCMbr);

int findNNinF(PotentialLocation* p, RTNode* pNodeF, float pNNinF[4][2], int& nSpaceOverhead);

void voronoiCell(float* pfdnn, SBlockFile* psbf_p, RTree * rC, RTree * rP, RTree * rF);

//////////////////////////////////////////////////////////////////////////////////////

void dd1(VoronoiCell** pVC, /*int nVcStart,*/ int nVcEnd, RTNode* pNodeC, float* fCMbr, float* pfdnn, PotentialLocation& r);

void voronoiCell1(float* pfdnn, SBlockFile* psbf_p, RTree * rC, RTree * rP, RTree * rF);

#endif
