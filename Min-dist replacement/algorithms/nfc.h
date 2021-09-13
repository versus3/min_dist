// nfc.h
// Nearest Facility Circle Algorithm
//
// Yuan (Andy), Xue && Jianzhong Qi
// Date Created : 09/04/2010
// Last Modified: 16/06/2010

#ifndef __NFC
#define __NFC

#include <map>
#include <list>
#include "../simpleblockfile/sblockfile.h"

#include "../mndrtree/mndrtree.h"
#include "../mndrtree/mndrtnode.h"
#include "../mndrtree/mndentry.h"

using namespace std;

typedef struct ISNode{
	float fmaxDr;
	list<int> iS;
}ISNode;

typedef map<int, ISNode*> IS;

typedef map<int, float> MAP_DD;

//void join(RTNode* pNode1, float* fMbr1, RTNode* pNode2, float* fMbr2, 
//		  PotentialLocation* p, float* pfdnn, IS* pIS);  

void pisJoin(RTNode* pNode1, float* fMbr1, RTNode* pNode2, float* fMbr2, float* pfdnn, IS* pIS, int& nSpaceOverhead);  

void pisJoin1(RTNode* pNode1, float* fMbr1, RTNode* pNode2, float* fMbr2, float* pfdnn, MAP_DD* pDd, PotentialLocation* r);

void pisJoin2(RTNode* pNode1, float* fMbr1, RTNode* pNode2, float* fMbr2, float* pfdnn, float* pfDd, PotentialLocation* r);

void pisJoin3(RTNode* pNode1, float* fMbr1, CMNDRTNode* pNode2, float* fMbr2, float* pfdnn, float* pfDd, PotentialLocation* r);

//void pis(RTNode* pNodeP, RTNode* pNodeR, PotentialLocation* p, float* pfdnn, IS* pIS); // find the PIS of p[i]
int cmp(const void *a1, const void *a2); // used in qsort()

void printRTree(RTNode *n); // print the info in an rtree

void nfc(float* pfdnn, SBlockFile* psbf_p, RTree * rR, RTree * rP);
void nfc1(float* pfdnn, RTree * rR, RTree * rP);
void nfc2(float* pfdnn, RTree * rR, RTree * rP);		// Join Method

void nfc3(float* pfdnn, CMNDRTree * rR, RTree * rP);		// Join Method on Expended Client Tree

#endif
