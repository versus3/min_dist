#pragma warning(disable:4786)
// mnd_repquery.cpp
// Maximum NFC Distance Algorithm for the replacement query
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 01/10/2012
//		Last Modified: 01/10/2012

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"
#include "mnd_repquery.h"
#include "../stopwatch/stopwatch.h"	// stopwatch to record the time

/*
	Maximum NFD distance method for the replacement query, interface
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry)
		rF (R-tree on F), 
		rP (R-tree on P)
	Output: Min-dist replacement pair (facility, potential location)
*/
void mnd_repquery(CMNDRepqueryRTree* rC, RTree* rF, RTree* rP) {
	// store the result - optimal facility-location pair
	opt_facility_location_pair_t optPair;
	
	float* fMbr1, * fMbr2, * fMbr3;
	/******* Stage 2: Finding Optiomal Candidate Location ******/
	rC->load_root();
	rF->load_root();
	rP->load_root();
	
	optPair.p.maxDr = 0;
	fMbr1 = rC->root_ptr->get_mbr();
	fMbr2 = rC->root_ptr->get_mbr();
	fMbr3 = rP->root_ptr->get_mbr();
	float* pfDr = new float[rF->root_ptr->capacity * rP->root_ptr->capacity];

	mndJoin_repquery(rC->root_ptr, fMbr1, DATA_DOMAIN, rF->root_ptr, fMbr2, rP->root_ptr, fMbr3, pfDr, optPair);
	
	/******* Stage 3: Output ******/
	printf("dr{f[%d]=(%.2f, %.2f), p[%d]=(%.2f, %.2f)} = %.2f\n", 
		optPair.f.id, optPair.f.x, optPair.f.y, optPair.p.id, optPair.p.x, optPair.p.y, optPair.p.maxDr);
	
	rC->del_root();
	rP->del_root();
	rF->del_root();
	delete[] pfDr;
	delete[] fMbr1;
	delete[] fMbr2;
	delete[] fMbr3;
}


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
					 float* pfDr, opt_facility_location_pair_t& optPair) {
	int i, j, k;
	RTNode* pNodeF1;
	float fMbr[4];

	if (pNodeF->level == 0 && pNodeP->level == 0) {
		if (pNodeC == pNodeC->my_tree->root_ptr){
			memset(pfDr, 0, sizeof(float) * pNodeF->capacity * pNodeP->capacity);
		}
		if (pNodeC->level == 0) {
			float fDistCP, fDistCF;

			fMbr[0] = fMbrC[0] - fmnd;
			fMbr[1] = fMbrC[1] + fmnd;
			fMbr[2] = fMbrC[2] - fmnd;
			fMbr[3] = fMbrC[3] + fmnd;

			for (i = 0; i < pNodeF->num_entries; i++){
				for (j = 0; j < pNodeP->num_entries; j++) {
					if (ptInMbr(pNodeF->entries[i].bounces[0], pNodeF->entries[i].bounces[2], fMbr) || 
						ptInMbr(pNodeP->entries[j].bounces[0], pNodeP->entries[j].bounces[2], fMbr) ) {
						for (k = 0; k < pNodeC->num_entries; k++) {
							fDistCP = pointDist(pNodeC->entries[k].bounces[0], pNodeC->entries[k].bounces[2], 
								pNodeP->entries[j].bounces[0], pNodeP->entries[j].bounces[2]);
							fDistCF = pointDist(pNodeC->entries[k].bounces[0], pNodeC->entries[k].bounces[2], 
								pNodeF->entries[i].bounces[0], pNodeF->entries[i].bounces[2]);

							if (fDistCP <= pNodeC->entries[k].fmnd) {
								pfDr[i*pNodeP->num_entries+j] += (pNodeC->entries[k].fmnd - fDistCP);
							} else if (fDistCP >= pNodeC->entries[k].f2mnd) {
								if (fDistCF > pNodeC->entries[k].fmnd) {
									pfDr[i*pNodeP->num_entries+j] += 0;
								} else if (fDistCF < pNodeC->entries[k].f2mnd) {
									pfDr[i*pNodeP->num_entries+j] += (pNodeC->entries[k].fmnd - pNodeC->entries[k].f2mnd);
								}
							} else {
								if (fDistCF > pNodeC->entries[k].fmnd) {
									pfDr[i*pNodeP->num_entries+j] += 0;
								} else if (fDistCF < pNodeC->entries[k].f2mnd) {
									pfDr[i*pNodeP->num_entries+j] += (pNodeC->entries[k].fmnd - fDistCP);
								}
							}
						}
						if (pfDr[i*pNodeP->num_entries+j] > optPair.p.maxDr) {
							optPair.f.id = pNodeF->entries[i].son;
							optPair.f.x = pNodeF->entries[i].bounces[0];
							optPair.f.y = pNodeF->entries[i].bounces[2];

							optPair.p.id = pNodeP->entries[j].son;
							optPair.p.x = pNodeP->entries[j].bounces[0];
							optPair.p.y = pNodeP->entries[j].bounces[2];

							optPair.p.maxDr = pfDr[i*pNodeP->num_entries+j];
						}
					}
				}
			}
		} else {
			for(i = 0; i < pNodeC->num_entries; i++){
				fMbr[0] = pNodeC->entries[i].bounces[0] - pNodeC->entries[i].fmnd;
				fMbr[1] = pNodeC->entries[i].bounces[1] + pNodeC->entries[i].fmnd;
				fMbr[2] = pNodeC->entries[i].bounces[2] - pNodeC->entries[i].fmnd;
				fMbr[3] = pNodeC->entries[i].bounces[3] + pNodeC->entries[i].fmnd;
				
				if (mbrIntersect(fMbr, fMbrF) || mbrIntersect(fMbr, fMbrP)) {
					mndJoin_repquery(pNodeC->entries[i].get_son(), pNodeC->entries[i].bounces, pNodeC->entries[i].fmnd, pNodeF, fMbrF, 
										pNodeP, fMbrP, pfDr, optPair);
					pNodeC->entries[i].del_son();				
				}
			}	
		}

	} else if (pNodeF->level == 0 && pNodeP->level != 0) {
		for(i = 0; i < pNodeP->num_entries; i++){
			mndJoin_repquery(pNodeC, fMbrC, 0, pNodeF, fMbrF, 
					pNodeP->entries[i].get_son(), pNodeP->entries[i].bounces, pfDr, optPair);
			pNodeP->entries[i].del_son();
		}
	} else if (pNodeF->level != 0 && pNodeP->level == 0) {
		for(i = 0; i < pNodeF->num_entries; i++){
			mndJoin_repquery(pNodeC, fMbrC, 0, pNodeF->entries[i].get_son(), pNodeF->entries[i].bounces, 
					pNodeP, fMbrP, pfDr, optPair);
			pNodeF->entries[i].del_son();
		}
	} else {
		for(i = 0; i < pNodeF->num_entries; i++){
			pNodeF1 = pNodeF->entries[i].get_son();
			for(j = 0; j < pNodeP->num_entries; j++){
				mndJoin_repquery(pNodeC, fMbrC, 0, pNodeF1, pNodeF->entries[i].bounces, 
					pNodeP->entries[j].get_son(), pNodeP->entries[j].bounces, pfDr, optPair);
				pNodeP->entries[j].del_son();
			}
			pNodeF->entries[i].del_son();
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// with pruning
/*
	Maximum NFD distance method (with pruning on F&P) for the replacement query, interface
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry)
		rF (R-tree on F), 
		rP (R-tree on P)
	Output: Min-dist replacement pair (facility, potential location)
*/
void mnd_repquery_with_pruning(CMNDRepqueryRTree* rC, RTree* rF, RTree* rP) {
	// store the result - optimal facility-location pair
	opt_facility_location_pair_t optPair;
	
	float* fMbr1, * fMbr2, * fMbr3;
	/******* Stage 2: Finding Optiomal Candidate Location ******/
	rC->load_root();
	rF->load_root();
	rP->load_root();
	
	
	optPair.p.maxDr = 0;
	fMbr1 = rC->root_ptr->get_mbr();
	fMbr2 = rC->root_ptr->get_mbr();
	fMbr3 = rP->root_ptr->get_mbr();
	float* pfDr = new float[rF->root_ptr->capacity * rP->root_ptr->capacity];

	int nPruneCounter = 0;
	start(); 
	pruning_dateset(rC->root_ptr, fMbr1, DATA_DOMAIN, rF->root_ptr, fMbr2, nPruneCounter);
	printf("%d\t%d\n", rF->num_of_dnodes, nPruneCounter); 
	stop(); 
	start(); 
	nPruneCounter = 0;
	pruning_dateset(rC->root_ptr, fMbr1, DATA_DOMAIN, rP->root_ptr, fMbr3, nPruneCounter);
	printf("%d\t%d\n", rP->num_of_dnodes, nPruneCounter);
	stop(); 

	// return;
	start();
	mndJoin_repquery_with_pruning(rC->root_ptr, fMbr1, 0, rF->root_ptr, fMbr2, rP->root_ptr, fMbr3, pfDr, optPair);
	stop(); 

	/******* Stage 3: Output ******/
	printf("dr{f[%d]=(%.2f, %.2f), p[%d]=(%.2f, %.2f)} = %.2f\n", 
		optPair.f.id, optPair.f.x, optPair.f.y, optPair.p.id, optPair.p.x, optPair.p.y, optPair.p.maxDr);
	
	rC->del_root();
	rP->del_root();
	rF->del_root();
	delete[] pfDr;
	delete[] fMbr1;
	delete[] fMbr2;
	delete[] fMbr3;
}

/*
	Pruning unpromising F & P
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry)
		fMbrC (MBR of pNodeC),
		pNodePrune (RTNode in rF or rP),
		fMbrPrune (MBR of pNodePrune),
	Output: Min-dist replacement pair (facility, potential location)
*/
int pruning_dateset(CMNDRepqueryRTNode* pNodeC, float* fMbrC, float fmnd, RTNode* pNodePrune, float* fMbrPrune, int& nPruneCounter) {
	int i, j;
	float fMbr[4], fMbr1[4];
	int flag;
	int ret = 0;

	if (pNodeC->level == 0 && pNodePrune->level == 0) {
		for (i = 0; i < pNodePrune->num_entries; i++) {
			fMbr1[0] = fMbrC[0] - fmnd;
			fMbr1[1] = fMbrC[1] + fmnd;
			fMbr1[2] = fMbrC[2] - fmnd;
			fMbr1[3] = fMbrC[3] + fmnd;
						
			if (mbrIntersect(fMbr1, pNodePrune->entries[i].bounces)) {
				
				for (j = 0; j < pNodeC->num_entries; j++) {
					fMbr[0] = pNodeC->entries[j].bounces[0] - pNodeC->entries[j].fmnd;
					fMbr[1] = pNodeC->entries[j].bounces[1] + pNodeC->entries[j].fmnd;
					fMbr[2] = pNodeC->entries[j].bounces[2] - pNodeC->entries[j].fmnd;
					fMbr[3] = pNodeC->entries[j].bounces[3] + pNodeC->entries[j].fmnd;

					if (ptInMbr(pNodePrune->entries[i].bounces[0], pNodePrune->entries[i].bounces[2], fMbr)) {
						pNodePrune->entries[i].nn_num = 1;	
						nPruneCounter++;
						ret = 1;
						break;
					}
				}
			}
		}
	} else if (pNodeC->level == 0) {
		fMbr[0] = fMbrC[0] - fmnd;
		fMbr[1] = fMbrC[1] + fmnd;
		fMbr[2] = fMbrC[2] - fmnd;
		fMbr[3] = fMbrC[3] + fmnd;

		for (i = 0; i < pNodePrune->num_entries; i++) {	
			if (mbrIntersect(fMbr, pNodePrune->entries[i].bounces)) {
				flag = pruning_dateset (pNodeC, fMbrC, fmnd, pNodePrune->entries[i].get_son(), pNodePrune->entries[i].bounces, nPruneCounter);
				pNodePrune->entries[i].del_son();
				if (flag == 1) {
					pNodePrune->entries[i].nn_num = 1;
					ret = 1;
				}
			}
		}
		return ret;
	} else if (pNodePrune->level == 0) {
		for (i = 0; i < pNodeC->num_entries; i++) {	
			fMbr[0] = pNodeC->entries[i].bounces[0] - pNodeC->entries[i].fmnd;
			fMbr[1] = pNodeC->entries[i].bounces[1] + pNodeC->entries[i].fmnd;
			fMbr[2] = pNodeC->entries[i].bounces[2] - pNodeC->entries[i].fmnd;
			fMbr[3] = pNodeC->entries[i].bounces[3] + pNodeC->entries[i].fmnd;

			if (mbrIntersect(fMbr, fMbrPrune)) {
				flag = pruning_dateset (pNodeC->entries[i].get_son(), pNodeC->entries[i].bounces, pNodeC->entries[i].fmnd, pNodePrune, fMbrPrune, nPruneCounter);
				pNodeC->entries[i].del_son();
				if (flag == 1) {
					return 1;
				}				
			}
		}

		return 0;
	} else {
		
		for (i = 0; i < pNodePrune->num_entries; i++) {
			fMbr1[0] = fMbrC[0] - fmnd;
			fMbr1[1] = fMbrC[1] + fmnd;
			fMbr1[2] = fMbrC[2] - fmnd;
			fMbr1[3] = fMbrC[3] + fmnd;
						
			if (mbrIntersect(fMbr1, pNodePrune->entries[i].bounces)) {
				
				for (j = 0; j < pNodeC->num_entries; j++) {
					fMbr[0] = pNodeC->entries[j].bounces[0] - pNodeC->entries[j].fmnd;
					fMbr[1] = pNodeC->entries[j].bounces[1] + pNodeC->entries[j].fmnd;
					fMbr[2] = pNodeC->entries[j].bounces[2] - pNodeC->entries[j].fmnd;
					fMbr[3] = pNodeC->entries[j].bounces[3] + pNodeC->entries[j].fmnd;

					if (mbrIntersect(fMbr, pNodePrune->entries[i].bounces)) {

						flag = pruning_dateset (pNodeC->entries[j].get_son(), pNodeC->entries[j].bounces, pNodeC->entries[j].fmnd, 
							pNodePrune->entries[i].get_son(), pNodePrune->entries[i].bounces, nPruneCounter);
						pNodeC->entries[j].del_son();
						if (flag == 1) {
							pNodePrune->entries[i].nn_num = 1;	
						}
						break;
					}
				}
				pNodePrune->entries[i].del_son();
			}
		}
	}
	return ret;
}

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
								 float* pfDr, opt_facility_location_pair_t& optPair) {
	int i, j, k;
	RTNode* pNodeF1;
	float fMbr[4];

	if (pNodeF->level == 0 && pNodeP->level == 0) {
		if (pNodeC == pNodeC->my_tree->root_ptr){
			memset(pfDr, 0, sizeof(float) * pNodeF->capacity * pNodeP->capacity);
		}
		if (pNodeC->level == 0) {
			float fDistCP, fDistCF;

			fMbr[0] = fMbrC[0] - fmnd;
			fMbr[1] = fMbrC[1] + fmnd;
			fMbr[2] = fMbrC[2] - fmnd;
			fMbr[3] = fMbrC[3] + fmnd;

			for (i = 0; i < pNodeF->num_entries; i++){
				if (pNodeF->entries[i].nn_num == 1) {
					continue;
				}
				for (j = 0; j < pNodeP->num_entries; j++) {
					if (pNodeP->entries[j].nn_num == 1) {
						continue;
					}
					if (ptInMbr(pNodeF->entries[i].bounces[0], pNodeF->entries[i].bounces[2], fMbr) || 
						ptInMbr(pNodeP->entries[j].bounces[0], pNodeP->entries[j].bounces[2], fMbr) ) {
						for (k = 0; k < pNodeC->num_entries; k++) {
							fDistCP = pointDist(pNodeC->entries[k].bounces[0], pNodeC->entries[k].bounces[2], 
								pNodeP->entries[j].bounces[0], pNodeP->entries[j].bounces[2]);
							fDistCF = pointDist(pNodeC->entries[k].bounces[0], pNodeC->entries[k].bounces[2], 
								pNodeF->entries[i].bounces[0], pNodeF->entries[i].bounces[2]);

							if (fDistCP <= pNodeC->entries[k].fmnd) {
								pfDr[i*pNodeP->num_entries+j] += (pNodeC->entries[k].fmnd - fDistCP);
							} else if (fDistCP >= pNodeC->entries[k].f2mnd) {
								if (fDistCF > pNodeC->entries[k].fmnd) {
									pfDr[i*pNodeP->num_entries+j] += 0;
								} else if (fDistCF < pNodeC->entries[k].f2mnd) {
									pfDr[i*pNodeP->num_entries+j] += (pNodeC->entries[k].fmnd - pNodeC->entries[k].f2mnd);
								}
							} else {
								if (fDistCF > pNodeC->entries[k].fmnd) {
									pfDr[i*pNodeP->num_entries+j] += 0;
								} else if (fDistCF < pNodeC->entries[k].f2mnd) {
									pfDr[i*pNodeP->num_entries+j] += (pNodeC->entries[k].fmnd - fDistCP);
								}
							}
						}
						if (pfDr[i*pNodeP->num_entries+j] > optPair.p.maxDr) {
							optPair.f.id = pNodeF->entries[i].son;
							optPair.f.x = pNodeF->entries[i].bounces[0];
							optPair.f.y = pNodeF->entries[i].bounces[2];

							optPair.p.id = pNodeP->entries[j].son;
							optPair.p.x = pNodeP->entries[j].bounces[0];
							optPair.p.y = pNodeP->entries[j].bounces[2];

							optPair.p.maxDr = pfDr[i*pNodeP->num_entries+j];
						}
					}
				}
			}
		} else {
			for(i = 0; i < pNodeC->num_entries; i++){
				fMbr[0] = pNodeC->entries[i].bounces[0] - pNodeC->entries[i].fmnd;
				fMbr[1] = pNodeC->entries[i].bounces[1] + pNodeC->entries[i].fmnd;
				fMbr[2] = pNodeC->entries[i].bounces[2] - pNodeC->entries[i].fmnd;
				fMbr[3] = pNodeC->entries[i].bounces[3] + pNodeC->entries[i].fmnd;
				
				if (mbrIntersect(fMbr, fMbrF) || mbrIntersect(fMbr, fMbrP)) {
					mndJoin_repquery_with_pruning(pNodeC->entries[i].get_son(), pNodeC->entries[i].bounces, pNodeC->entries[i].fmnd, pNodeF, fMbrF, 
										pNodeP, fMbrP, pfDr, optPair);
					pNodeC->entries[i].del_son();				
				}
			}	
		}

	} else if (pNodeF->level == 0 && pNodeP->level != 0) {
		for(i = 0; i < pNodeP->num_entries; i++){
			if (pNodeP->entries[i].nn_num == 1) {
				mndJoin_repquery_with_pruning(pNodeC, fMbrC, 0, pNodeF, fMbrF, 
						pNodeP->entries[i].get_son(), pNodeP->entries[i].bounces, pfDr, optPair);
				pNodeP->entries[i].del_son();
			}
		}
	} else if (pNodeF->level != 0 && pNodeP->level == 0) {
		for(i = 0; i < pNodeF->num_entries; i++){
			if (pNodeF->entries[i].nn_num == 1) {
				mndJoin_repquery_with_pruning(pNodeC, fMbrC, 0, pNodeF->entries[i].get_son(), pNodeF->entries[i].bounces, 
						pNodeP, fMbrP, pfDr, optPair);
				pNodeF->entries[i].del_son();
			}
		}
	} else {
		for(i = 0; i < pNodeF->num_entries; i++){
			if (pNodeF->entries[i].nn_num == 1) {
				pNodeF1 = pNodeF->entries[i].get_son();
				for(j = 0; j < pNodeP->num_entries; j++){
					if (pNodeP->entries[j].nn_num == 1) {
						mndJoin_repquery_with_pruning(pNodeC, fMbrC, 0, pNodeF1, pNodeF->entries[i].bounces, 
								pNodeP->entries[j].get_son(), pNodeP->entries[j].bounces, pfDr, optPair);
						pNodeP->entries[j].del_son();					
					}
				}
				pNodeF->entries[i].del_son();
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SS + MND
/*
	SS (for f-p pair generation) + MND (for influenced client identification) method
	Input: rC (CMNDRepqueryRTree on C, with mnd and 2mnd on each entry),
		psbf_f (facility set), psbf_p (Potential location set)
	Output: Min-dist replacement pair
*/
void ss_mnd_repquery(CMNDRepqueryRTree* rC, SBlockFile* psbf_f, SBlockFile* psbf_p) {
	// store the result - optimal replacement pair (facility, potential location)
	opt_facility_location_pair_t optPair;
	optPair.p.maxDr = 0;

	// Temporary variables
	int i, j;
	float* fMbr1;
	float fMbr2[4], fMbr3[4];
	float* fMaxDr = new float[psbf_f->m_nItemPerPage*psbf_p->m_nItemPerPage];

	// find dr(p_i) for each p_i by finding clients who consider p_i as their new NN
	Facility* f = new Facility[psbf_f->m_nItemPerPage];						// An array to hold a block of facilities read from the Facility set data file
	PotentialLocation* p = new PotentialLocation[psbf_p->m_nItemPerPage];	// An array to hold a block of potential locations read from the Potential Location set data file
	int nFRecord, nPRecord;													// Number of records read from the data file

	rC->load_root();
	fMbr1 = rC->root_ptr->get_mbr();

	nFRecord = nPRecord = 0;
	psbf_f->rewind();
	for (i = 0; i < psbf_f->m_nItemTotal; i+= nFRecord) {
		nFRecord = psbf_f->readBlock(f);
	
		psbf_p->rewind();
		for (j = 0; j < psbf_p->m_nItemTotal; j += nPRecord) {
			nPRecord = psbf_p->readBlock(p);
			
			memset(fMaxDr, 0, sizeof(float)*psbf_f->m_nItemPerPage*psbf_p->m_nItemPerPage);
			getMbr(f, nFRecord, fMbr2);
			getMbr(p, nPRecord, fMbr3);

			ss_mndJoin_repquery(rC->root_ptr, fMbr1, DATA_DOMAIN, f, nFRecord, fMbr2, p, nPRecord, fMbr3, fMaxDr, optPair);
		}
	}

	// Output
	printf("dr{f[%d]=(%.2f, %.2f), p[%d]=(%.2f, %.2f)} = %.2f\n", 
		optPair.f.id, optPair.f.x, optPair.f.y, optPair.p.id, optPair.p.x, optPair.p.y, optPair.p.maxDr);

	delete[] f;
	delete[] p;

	delete[] fMaxDr;
	delete[] fMbr1;
}

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
						PotentialLocation* pNodeP, int nPNum, float* fMbrP, float* pfDr, opt_facility_location_pair_t& optPair) {
	
	int i, j, k;
	float fMbr[4];						
	float fDistCP, fDistCF;

	if (pNodeC->level == 0) {
		
		fMbr[0] = fMbrC[0] - fmnd;
		fMbr[1] = fMbrC[1] + fmnd;
		fMbr[2] = fMbrC[2] - fmnd;
		fMbr[3] = fMbrC[3] + fmnd;

		for (i = 0; i < nFNum; i++){
			for (j = 0; j < nPNum; j++) {
				if (ptInMbr(pNodeF[i].x, pNodeF[i].y, fMbr) || 
					ptInMbr(pNodeP[j].x, pNodeP[j].y, fMbr) ) {
					for (k = 0; k < pNodeC->num_entries; k++) {
						fDistCP = pointDist(pNodeC->entries[k].bounces[0], pNodeC->entries[k].bounces[2], 
							pNodeP[j].x, pNodeP[j].y);
						fDistCF = pointDist(pNodeC->entries[k].bounces[0], pNodeC->entries[k].bounces[2], 
							pNodeF[i].x, pNodeF[i].y);

						if (fDistCP <= pNodeC->entries[k].fmnd) {
							pfDr[i*nPNum+j] += (pNodeC->entries[k].fmnd - fDistCP);
						} else if (fDistCP >= pNodeC->entries[k].f2mnd) {
							if (fDistCF > pNodeC->entries[k].fmnd) {
								pfDr[i*nPNum+j] += 0;
							} else if (fDistCF < pNodeC->entries[k].f2mnd) {
								pfDr[i*nPNum+j] += (pNodeC->entries[k].fmnd - pNodeC->entries[k].f2mnd);
							}
						} else {
							if (fDistCF > pNodeC->entries[k].fmnd) {
								pfDr[i*nPNum+j] += 0;
							} else if (fDistCF < pNodeC->entries[k].f2mnd) {
								pfDr[i*nPNum+j] += (pNodeC->entries[k].fmnd - fDistCP);
							}
						}
					}
					if (pfDr[i*nPNum+j] > optPair.p.maxDr) {
						optPair.f = pNodeF[i];
						optPair.p = pNodeP[j];

						optPair.p.maxDr = pfDr[i*nPNum+j];
					}
				}
			}
		}
	} else {
		for(i = 0; i < pNodeC->num_entries; i++){
			fMbr[0] = pNodeC->entries[i].bounces[0] - pNodeC->entries[i].fmnd;
			fMbr[1] = pNodeC->entries[i].bounces[1] + pNodeC->entries[i].fmnd;
			fMbr[2] = pNodeC->entries[i].bounces[2] - pNodeC->entries[i].fmnd;
			fMbr[3] = pNodeC->entries[i].bounces[3] + pNodeC->entries[i].fmnd;

			if (mbrIntersect(fMbr, fMbrF) || mbrIntersect(fMbr, fMbrP)) {
				ss_mndJoin_repquery(pNodeC->entries[i].get_son(), pNodeC->entries[i].bounces, pNodeC->entries[i].fmnd, pNodeF, nFNum, fMbrF, 
					pNodeP, nPNum, fMbrP, pfDr, optPair);
				pNodeC->entries[i].del_son();				
			}
		}	
	}
}
