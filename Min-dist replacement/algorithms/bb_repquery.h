// bb_repquery.h
// Branch and Bound Algorithm for the replacement query
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 25/09/2012
//		Last Modified: 25/09/2012

#ifndef __BRANCH_BOUND_REPQUERY_H
#define __BRANCH_BOUND_REPQUERY_H

#include <queue>
#include <vector>
#include "../point/point.h"

#include "../mndrepqueryrtree/mndrepqueryrtree.h"
#include "../mndrepqueryrtree/mndrepqueryrtnode.h"
#include "../mndrepqueryrtree/mndrepqueryentry.h"

using namespace std;

// Branch and bound for the replacement query
/*
	Branch and bound method
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry),
		rF_BB_repquery (CMNDRepqueryRTree on F, with mnd and maxdr), 
		rP_BB_repquery (CMNDRepqueryRTree on P, with mnd and maxdr)
	Output: Min-dist replacement pair
*/
void bb_repquery(CMNDRepqueryRTree* rC_MND_repquery, CMNDRepqueryRTree* rF_BB_repquery, CMNDRepqueryRTree* rP_BB_repquery);

void bbjoin_repquery(CMNDRepqueryRTNode* pNodeF, CMNDRepqueryEntry* pEntryF, CMNDRepqueryRTNode* pNodeP, CMNDRepqueryEntry* pEntryP, 
					 CMNDRepqueryRTNode* pNodeC, opt_facility_location_pair_t& optPair);

void bbjoin_repquery_entry_processing(CMNDRepqueryEntry* pEntryF, CMNDRepqueryEntry* pEntryP, CMNDRepqueryRTNode* pNodeC, float& fMaxDr);

#endif
