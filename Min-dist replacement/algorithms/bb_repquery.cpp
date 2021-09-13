#pragma warning(disable:4786)
// bb_repquery.cpp
// Branch and Bound Algorithm for the replacement query
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 25/09/2012
//		Last Modified: 25/09/2012

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../point/point.h"
#include "../rtree/rtree.h"
#include "../rtree/rtnode.h"
#include "../rtree/entry.h"

#include "defs.h"
#include "bb_repquery.h"

// Branch and bound for the replacement query
/*
	Branch and bound method
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry),
		rF_BB_repquery (CMNDRepqueryRTree on F, with mnd and maxdr), 
		rP_BB_repquery (CMNDRepqueryRTree on P, with mnd and maxdr)
	Output: Min-dist replacement pair
*/
void bb_repquery(CMNDRepqueryRTree* rC_MND_repquery, CMNDRepqueryRTree* rF_BB_repquery, CMNDRepqueryRTree* rP_BB_repquery) {
	
	// store the result - optimal facility-location pair
	opt_facility_location_pair_t optPair;

	/******* Stage 2: Finding Optiomal Facility-Location Pair*******/
	rC_MND_repquery->load_root();
	rF_BB_repquery->load_root();
	rP_BB_repquery->load_root();

	optPair.p.maxDr = 0;
		
	bbjoin_repquery(rF_BB_repquery->root_ptr, NULL, rP_BB_repquery->root_ptr, NULL, rC_MND_repquery->root_ptr, optPair);

	/******* Stage 3: Output *******/
	printf("dr{f[%d]=(%.2f, %.2f), p[%d]=(%.2f, %.2f)} = %.2f\n", 
		optPair.f.id, optPair.f.x, optPair.f.y, optPair.p.id, optPair.p.x, optPair.p.y, optPair.p.maxDr);

	rC_MND_repquery->del_root();
	rF_BB_repquery->del_root();
	rP_BB_repquery->del_root();
}

void bbjoin_repquery(CMNDRepqueryRTNode* pNodeF, CMNDRepqueryEntry* pEntryF, CMNDRepqueryRTNode* pNodeP, CMNDRepqueryEntry* pEntryP, 
					 CMNDRepqueryRTNode* pNodeC, opt_facility_location_pair_t& optPair) {
	
	int i, j;
	float fMbrF[4];
	float* pfMbrF, * pfMbrP;
	int nMaxID;
	float fMaxDr;

	if (pNodeF->level == 0 && pNodeP->level == 0) {
		for (i = 0; i < pNodeF->num_entries; i++) {
			fMbrF[0] = pNodeF->entries[i].bounces[0] - pNodeF->entries[i].f2mnd;
			fMbrF[1] = pNodeF->entries[i].bounces[1] + pNodeF->entries[i].f2mnd;
			fMbrF[2] = pNodeF->entries[i].bounces[2] - pNodeF->entries[i].f2mnd;
			fMbrF[3] = pNodeF->entries[i].bounces[3] + pNodeF->entries[i].f2mnd;

			for (j = 0; j < pNodeP->num_entries; j++) {
				if (pNodeP->entries[j].fmnd - pNodeF->entries[i].fmnd >= optPair.p.maxDr || 
					(pNodeP->entries[j].fmnd >= optPair.p.maxDr && mbrIntersect(fMbrF, pNodeP->entries[j].bounces)) ){
					
					fMaxDr = 0;
					bbjoin_repquery_entry_processing(&(pNodeF->entries[i]), &(pNodeP->entries[j]), pNodeC, fMaxDr);
					if (optPair.p.maxDr <= fMaxDr) {
						optPair.f.id = pNodeF->entries[i].son;
						optPair.f.x = pNodeF->entries[i].bounces[0];
						optPair.f.y = pNodeF->entries[i].bounces[2];

						optPair.p.id = pNodeP->entries[j].son;
						optPair.p.x = pNodeP->entries[j].bounces[0];
						optPair.p.y = pNodeP->entries[j].bounces[2];
						optPair.p.maxDr = fMaxDr;
					}
				}
			}
		}
	} else if (pNodeF->level == 0) {
		if (!pEntryF) {
			pfMbrF = pNodeF->get_mbr();

			fMbrF[0] = pfMbrF[0] - DATA_DOMAIN;
			fMbrF[1] = pfMbrF[1] + DATA_DOMAIN;
			fMbrF[2] = pfMbrF[2] - DATA_DOMAIN;
			fMbrF[3] = pfMbrF[3] + DATA_DOMAIN;

			nMaxID = -1;
		
			for (i = 0; i < pNodeP->num_entries; i++) {	

				if (pNodeP->entries[i].fmnd >= optPair.p.maxDr && mbrIntersect(pfMbrF, pNodeP->entries[i].bounces)) {
					bbjoin_repquery(pNodeF, pEntryF, pNodeP->entries[i].get_son(), &(pNodeP->entries[i]), pNodeC, optPair);
					pNodeP->entries[i].del_son();
				} else if (nMaxID == -1 || pNodeP->entries[i].fmnd > fMaxDr){
					fMaxDr = pNodeP->entries[i].fmnd;
					nMaxID = i;
				}
			}

			if (nMaxID != -1 && fMaxDr >= optPair.p.maxDr) {
				if (fMaxDr > optPair.p.maxDr) {
					optPair.p.maxDr = fMaxDr;
				}
				bbjoin_repquery(pNodeF, pEntryF, pNodeP->entries[nMaxID].get_son(), &(pNodeP->entries[nMaxID]), pNodeC, optPair);
				pNodeP->entries[nMaxID].del_son();
			}
			delete[] pfMbrF;
		} else {
			nMaxID = -1;
			
			fMbrF[0] = pEntryF->bounces[0] - pEntryF->f2mnd;
			fMbrF[1] = pEntryF->bounces[1] + pEntryF->f2mnd;
			fMbrF[2] = pEntryF->bounces[2] - pEntryF->f2mnd;
			fMbrF[3] = pEntryF->bounces[3] + pEntryF->f2mnd;

			for (i = 0; i < pNodeP->num_entries; i++) {	

				if (pNodeP->entries[i].fmnd >= optPair.p.maxDr && mbrIntersect(fMbrF, pNodeP->entries[i].bounces)) {
					bbjoin_repquery(pNodeF, pEntryF, pNodeP->entries[i].get_son(), &(pNodeP->entries[i]), pNodeC, optPair);
					pNodeP->entries[i].del_son();
				} else if (nMaxID == -1 || pNodeP->entries[i].fmnd - pEntryF->f2mnd >= fMaxDr){
					fMaxDr = pNodeP->entries[i].fmnd;
					nMaxID = i;
				}
			}

			if (nMaxID != -1 && fMaxDr >= optPair.p.maxDr) {
				if (fMaxDr > optPair.p.maxDr) {
					optPair.p.maxDr = fMaxDr;
				}				
				bbjoin_repquery(pNodeF, pEntryF, pNodeP->entries[nMaxID].get_son(), &(pNodeP->entries[nMaxID]), pNodeC, optPair);
				pNodeP->entries[nMaxID].del_son();
			}
		}
	} else if (pNodeP->level == 0) {
		if (!pEntryP) {
			pfMbrP = pNodeP->get_mbr();
			nMaxID = -1;
		
			for (i = 0; i < pNodeF->num_entries; i++) {	
				fMbrF[0] = pNodeF->entries[i].bounces[0] - pNodeF->entries[i].f2mnd;
				fMbrF[1] = pNodeF->entries[i].bounces[1] + pNodeF->entries[i].f2mnd;
				fMbrF[2] = pNodeF->entries[i].bounces[2] - pNodeF->entries[i].f2mnd;
				fMbrF[3] = pNodeF->entries[i].bounces[3] + pNodeF->entries[i].f2mnd;

				if (mbrIntersect(fMbrF, pfMbrP)) {
					bbjoin_repquery(pNodeF->entries[i].get_son(), &(pNodeF->entries[i]), pNodeP, pEntryP, pNodeC, optPair);
					pNodeF->entries[i].del_son();
				} else if (nMaxID == -1 || pNodeF->entries[i].fmnd < fMaxDr){
					fMaxDr = pNodeF->entries[i].fmnd;
					nMaxID = i;
				}
			}

			if (nMaxID != -1) {
				bbjoin_repquery(pNodeF->entries[nMaxID].get_son(), &(pNodeF->entries[nMaxID]), pNodeP, pEntryP, pNodeC, optPair);
				pNodeF->entries[i].del_son();
			}
			delete[] pfMbrP;
		} else {
			nMaxID = -1;
			
			for (i = 0; i < pNodeF->num_entries; i++) {	
				fMbrF[0] = pNodeF->entries[i].bounces[0] - pNodeF->entries[i].f2mnd;
				fMbrF[1] = pNodeF->entries[i].bounces[1] + pNodeF->entries[i].f2mnd;
				fMbrF[2] = pNodeF->entries[i].bounces[2] - pNodeF->entries[i].f2mnd;
				fMbrF[3] = pNodeF->entries[i].bounces[3] + pNodeF->entries[i].f2mnd;

				if (mbrIntersect(fMbrF, pEntryP->bounces)) {
					bbjoin_repquery(pNodeF->entries[i].get_son(), &(pNodeF->entries[i]), pNodeP, pEntryP, pNodeC, optPair);
					pNodeF->entries[i].del_son();
				} else if (nMaxID == -1 || pEntryP->fmnd - pNodeF->entries[i].fmnd < fMaxDr){
					fMaxDr = pEntryP->fmnd - pNodeF->entries[i].fmnd;
					nMaxID = i;
				}
			}
			if (nMaxID != -1 && fMaxDr >= optPair.p.maxDr) {
				if (fMaxDr > optPair.p.maxDr) {
					optPair.p.maxDr = fMaxDr;
				}
				bbjoin_repquery(pNodeF->entries[nMaxID].get_son(), &(pNodeF->entries[nMaxID]), pNodeP, pEntryP, pNodeC, optPair);
				pNodeF->entries[i].del_son();
			}
		}
	
	} else {
		for (i = 0; i < pNodeF->num_entries; i++) {
			fMbrF[0] = pNodeF->entries[i].bounces[0] - pNodeF->entries[i].f2mnd;
			fMbrF[1] = pNodeF->entries[i].bounces[1] + pNodeF->entries[i].f2mnd;
			fMbrF[2] = pNodeF->entries[i].bounces[2] - pNodeF->entries[i].f2mnd;
			fMbrF[3] = pNodeF->entries[i].bounces[3] + pNodeF->entries[i].f2mnd;
			for (j = 0; j < pNodeP->num_entries; j++) {
				if (pNodeP->entries[j].fmnd - pNodeF->entries[i].fmnd >= optPair.p.maxDr || 
					(pNodeP->entries[j].fmnd >= optPair.p.maxDr && mbrIntersect(fMbrF, pNodeP->entries[j].bounces)) ) {
					if (pNodeP->entries[j].fmnd - pNodeF->entries[i].fmnd >= optPair.p.maxDr) {
						optPair.p.maxDr = pNodeP->entries[j].fmnd - pNodeF->entries[i].fmnd;
					}
					bbjoin_repquery(pNodeF->entries[i].get_son(), &(pNodeF->entries[i]), 
						pNodeP->entries[j].get_son(), &(pNodeP->entries[j]), pNodeC, optPair);				
					pNodeP->entries[j].del_son();
				}
			}
			pNodeF->entries[i].del_son();
		}
	}
}

void bbjoin_repquery_entry_processing(CMNDRepqueryEntry* pEntryF, CMNDRepqueryEntry* pEntryP, CMNDRepqueryRTNode* pNodeC, float& fMaxDr) {
	
	int i;
	float fMbr[4];
	float fDistCP, fDistCF;

	if (pNodeC->level == 0) {
		for (i = 0; i < pNodeC->num_entries; i++) {
			fDistCP = pointDist(pNodeC->entries[i].bounces[0], pNodeC->entries[i].bounces[2], 
				pEntryP->bounces[0], pEntryP->bounces[2]);
			fDistCF = pointDist(pNodeC->entries[i].bounces[0], pNodeC->entries[i].bounces[2], 
				pEntryF->bounces[0], pEntryF->bounces[2]);

			if (fDistCP <= pNodeC->entries[i].fmnd) {
				fMaxDr += (pNodeC->entries[i].fmnd - fDistCP);
			} else if (fDistCP >= pNodeC->entries[i].f2mnd) {
				if (fDistCF > pNodeC->entries[i].fmnd) {
					fMaxDr += 0;
				} else if (fDistCF < pNodeC->entries[i].f2mnd) {
					fMaxDr += (pNodeC->entries[i].fmnd - pNodeC->entries[i].f2mnd);
				}
			} else {
				if (fDistCF > pNodeC->entries[i].fmnd) {
					fMaxDr += 0;
				} else if (fDistCF < pNodeC->entries[i].f2mnd) {
					fMaxDr += (pNodeC->entries[i].fmnd - fDistCP);
				}
			}
		}
	} else {
		for (i = 0; i < pNodeC->num_entries; i++) {
			fMbr[0] = pNodeC->entries[i].bounces[0] - pNodeC->entries[i].fmnd;
			fMbr[1] = pNodeC->entries[i].bounces[1] + pNodeC->entries[i].fmnd;
			fMbr[2] = pNodeC->entries[i].bounces[2] - pNodeC->entries[i].fmnd;
			fMbr[3] = pNodeC->entries[i].bounces[3] + pNodeC->entries[i].fmnd;
			if (ptInMbr(pEntryF->bounces[0], pEntryF->bounces[2], fMbr) || ptInMbr(pEntryP->bounces[0], pEntryP->bounces[2], fMbr)){
				bbjoin_repquery_entry_processing(pEntryF, pEntryP, pNodeC->entries[i].get_son(), fMaxDr);
			}
			pNodeC->entries[i].del_son();
		}
	}
}