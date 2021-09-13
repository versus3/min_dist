// bb.h
// Branch and Bound Method
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 08/06/2010
//		Last Modified: 24/09/2012

#ifndef __BRANCH_BOUND_H
#define __BRANCH_BOUND_H

#include <queue>
#include <vector>
#include <map>
#include "../point/point.h"

#include "../bbrtree/bbrtree.h"
#include "../bbrtree/bbentry.h"
#include "../bbrtree/bbrtnode.h"

using namespace std;

typedef struct QpNode{
	Entry* pP;
	int nLevel;
	//vector<RTNode> Qc;
	double fMinDd;
	double fmaxDr;
	struct QpNode* pNext;
}QpNode;

typedef struct CNode{
	float fMBR/*[2]*/[4];
	//RTNode* pRRNode;
	double fMaxFDist;
}CNode;

typedef struct PQpNode{
	Entry* pP;
	vector<CBBEntry*> Qc;
	double fMinDd;
	double fmaxDr;
}PQpNode;

struct bbcmp{
    bool operator() (PQpNode* pNode1, PQpNode* pNode2){
		return pNode1->fmaxDr < pNode2->fmaxDr;
	}
};

typedef map<int, CNode*> ExtRC;

// Distance metrics
/* 
	Compute the min-dist between two MBRs
	Input: fMbrC, fMbrP (two MBRs)
	output: min-dist(fMbrC, fMbrP)
*/
double minDist(float* fMbrC, float* fMbrP);

/* 
	Compute the min-max-dist between a point and an MBR
	Input: x, y (a point), fMbrP (an MBR)
	output: min-max-dist((x, y), fMbrP)
*/
double minMaxDist(float x, float y, float* fMbrP);

/* 
	Compute the min-exist-dnn between two MBRs
	Input: fMbrC, fMbrP (two MBRs)
	output: min-exist-dnn(fMbrC, fMbrP)
*/
double minExistDNN(float* fMbrP, float* fMbrC);

void pis(QpNode* pQpNode, float* pfdnn, RTNode* pNodeC, ExtRC::iterator itr, double fMinDist, ExtRC* pExtRc);

double pisEntry(Entry* pCanLoc, float* pfdnn, RTNode* pNodeC, ExtRC* pExtRc);

/* 
	Branch and bound processing for a pair of nodes
	Input: pPNode (R-tree node of P), 
		fPMbr (MBR of pPNode),
		pCNode (BB-Rtree node of C), 
		fCMbr (MBR of pCNode),
		fMaxFDist (maxFDist of pCNode), 
		fpd (Pruning distance)
		r (Optimal location)
		pfdr (An array to hold all the dr values for the potential locations)
	output: min-dist location
*/
void bbJoin(RTNode* pPNode, float* fPMbr, CBBRTNode* pCNode, float* fCMbr, float fMaxFDist, float& fpd, PotentialLocation* r, float* pfdr);

/* 
	Branch and bound method
	Input: rC (bb-rtree on C), rp (r-tree on P)
	output: min-dist location
*/
void branchBound(CBBRTree* rC, RTree* rP);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An experimental version, implementing the bb-rtree with an r-tree and the extra information stored in a separated data structure
void branchBound_Ex(float* pfdnn, RTree* rC, /*RTree* rR, */RTree* rP, ExtRC* pExtRc);

void clearExtRc(ExtRC* pExtRc, RTNode* pNodeC);
double initExtRc(ExtRC* pExtRc, float* pfdnn, RTNode* pNodeC, float* fCMbr/*, RTNode* pNodeR, float* fRMbr*/);


#endif
