// mnd.h
// Maximum NFC Distance Algorithm
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 08/06/2010
//		Last Modified: 10/06/2010

#ifndef __MND_H
#define __MND_H

#include <vector>
#include <map>
#include "../point/point.h"
#include "../mndrtree/mndrtree.h"
#include "../mndrtree/mndrtnode.h"
#include "../mndrtree/mndentry.h"

#include "../rtree/rtree.h"
#include "../rtree/rtnode.h"
#include "../rtree/entry.h"

using namespace std;

void mdrJoin(RTNode* pNode1, float* fMbr1, CMNDRTNode* pNode2, float* fMbr2, float* pfDd, PotentialLocation* r);

void mdr(CMNDRTree * rC, RTree * rP);		// Join Method on Expended Client Tree

#endif
