// mnd_rep_query.h
// Maximum NFC Distance Algorithm for the replacement query
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 01/10/2012
//		Last Modified: 01/10/2012

#ifndef __MND_REPLACEMENT_QUERY_H
#define __MND_REPLACEMENT_QUERY_H

#include <vector>
#include <map>
#include "../point/point.h"
#include "../simpleblockfile/sblockfile.h"
#include "../mndrepqueryrtree/mndrepqueryrtree.h"
#include "../mndrepqueryrtree/mndrepqueryrtnode.h"
#include "../mndrepqueryrtree/mndrepqueryentry.h"

#include "../rtree/rtree.h"
#include "../rtree/rtnode.h"
#include "../rtree/entry.h"

using namespace std;

/*
	Maximum NFD distance method for the replacement query, interface
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry)
		rF (R-tree on F), 
		rP (R-tree on P)
	Output: Min-dist replacement pair (facility, potential location)
*/
void mnd_repquery(CMNDRepqueryRTree* rC, RTree* rF, RTree* rP);

/*
	Maximum NFD distance method for the replacement query, actuall processing
	Input: pNodeC (CMNDRepqueryRTNode in rC, with mnd and 2mnd on each entry)
		fMbrC (MBR of pNodeC),
		pNodeF (RTNode in rF),
		fMbrF (MBR of pNodeF),
		pNodeP (RTNode in rP),
		fMbrP (MBR of pNodeP),
		pfDr (Distance reduction),
		optPair (optimal facility location replacement pair)
	Output: optPair
*/
void mndJoin_repquery(CMNDRepqueryRTNode* pNodeC, float* fMbrC, float fmnd, RTNode* pNodeF, float* fMbrF, RTNode* pNodeP, float* fMbrP, 
					  float* pfDr, opt_facility_location_pair_t& optPair);


/*
	Maximum NFD distance method (with pruning on F&P) for the replacement query, interface
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry)
		rF (R-tree on F), 
		rP (R-tree on P)
	Output: Min-dist replacement pair (facility, potential location)
*/
void mnd_repquery_with_pruning(CMNDRepqueryRTree* rC, RTree* rF, RTree* rP);


/*
	Pruning unpromising F & P
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry)
		fMbrC (MBR of pNodeC),
		pNodePrune (RTNode in rF or rP),
		fMbrPrune (MBR of pNodePrune),
	Output: Min-dist replacement pair (facility, potential location)
*/
int pruning_dateset(CMNDRepqueryRTNode* pNodeC, float* fMbrC, float fmnd, RTNode* pNodePrune, float* fMbrPrune, int& nPruneCounter);

/*
	Maximum NFD distance method (with pruning on F&P) for the replacement query, actuall processing
	Input: pNodeC (CMNDRepqueryRTNode in rC, with mnd and 2mnd on each entry)
		fMbrC (MBR of pNodeC),
		pNodeF (RTNode in rF),
		fMbrF (MBR of pNodeF),
		pNodeP (RTNode in rP),
		fMbrP (MBR of pNodeP),
		pfDr (Distance reduction),
		optPair (optimal facility location replacement pair)
	Output: optPair
*/
void mndJoin_repquery_with_pruning(CMNDRepqueryRTNode* pNodeC, float* fMbrC, float fmnd, RTNode* pNodeF, float* fMbrF, RTNode* pNodeP, float* fMbrP, 
					  float* pfDr, opt_facility_location_pair_t& optPair);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SS + MND
/*
	SS (for f-p pair generation) + MND (for influenced client identification) method
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry),
		psbf_f (facility set), psbf_p (Potential location set)
	Output: Min-dist replacement pair
*/
void ss_mnd_repquery(CMNDRepqueryRTree* rC, SBlockFile* psbf_f, SBlockFile* psbf_p);

/*
	Maximum NFD distance method for the replacement query, actuall processing
	Input: pNodeC (CMNDRepqueryRTNode in rC, with mnd and 2mnd on each entry)
		fMbrC (MBR of pNodeC),
		pNodeF (A Block of F),
		nFNum (number of entries in pNodeF)
		fMbrF (MBR of pNodeF),
		pNodeP (a Block of P),
		nPNum (number of entries in pNodeP)
		fMbrP (MBR of pNodeP),
		pfDr (Distance reduction),
		optPair (optimal facility location replacement pair)
	Output: optPair
*/
void ss_mndJoin_repquery(CMNDRepqueryRTNode* pNodeC, float* fMbrC, float fmnd, Facility* pNodeF, int nFNum, float* fMbrF, 
						 PotentialLocation* pNodeP, int nPNum, float* fMbrP, float* pfDr, opt_facility_location_pair_t& optPair);

#endif
