#pragma warning(disable:4786)
// mnd.cpp
// Maximum NFC Distance Algorithm
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 08/06/2010
//		Last Modified: 10/06/2010

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mnd.h"


void mdr(CMNDRTree * rC, RTree * rP)
{
	// store the result - optimal candidate location
	PotentialLocation r;
	float* fMbr1, * fMbr2, fMbr3[4];
	float* fmnd = new float[rP->num_of_data];
	memset(fmnd, 0, sizeof(float)*rP->num_of_data);
	/******* Stage 2: Finding Optiomal Candidate Location *******/
	rC->load_root();
	rP->load_root();

	r.maxDr = 0;
	
	fMbr1 = rP->root_ptr->get_mbr();
	fMbr2 = rC->root_ptr->get_mbr();
	
	fMbr3[0] = fMbr2[0] - 1000;
	fMbr3[1] = fMbr2[1] + 1000;
	fMbr3[2] = fMbr2[2] - 1000;
	fMbr3[3] = fMbr2[3] + 1000;
	if(mbrIntersect(fMbr1, fMbr2)){
		mdrJoin(rP->root_ptr, fMbr1, rC->root_ptr, fMbr3, fmnd, &r);
	}
	
	/******* Stage 3: Output *******/
	if(r.maxDr > 0){
		printf("dd{p[%d]=(%.2f, %.2f)} = %.2f\n", r.id, r.x, r.y, r.maxDr);
	}

	delete[] fmnd;
	rC->del_root();
	rP->del_root();
}

void mdrJoin(RTNode* pNode1, float* fMbr1, CMNDRTNode* pNode2, float* fMbr2, float* pfDd, PotentialLocation* r)
{
	int i, j;
	RTNode* rtn;
	float fmnd;
	float fMbr[4];

	if(pNode1->level == 0){
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], fMbr2)) {
					for(j = 0; j < pNode2->num_entries; j++){
						fMbr[0] = pNode2->entries[j].bounces[0] - pNode2->entries[j].fmnd;
						fMbr[1] = pNode2->entries[j].bounces[1] + pNode2->entries[j].fmnd;
						fMbr[2] = pNode2->entries[j].bounces[2] - pNode2->entries[j].fmnd;
						fMbr[3] = pNode2->entries[j].bounces[3] + pNode2->entries[j].fmnd;
						if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], fMbr)) {
							fmnd = pNode2->entries[j].fmnd - pointDist(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], 
									pNode2->entries[j].bounces[0], pNode2->entries[j].bounces[2]);
							if(fmnd > 0){
								pfDd[pNode1->entries[i].son-1] += fmnd;
							}
						}
					}
					if(pfDd[pNode1->entries[i].son-1] > r->maxDr){
						r->maxDr = pfDd[pNode1->entries[i].son-1];
						r->id = pNode1->entries[i].son;
						r->x = pNode1->entries[i].bounces[0];
						r->y = pNode1->entries[i].bounces[2];
					}
				}
			}
		}else{
			for(i = 0; i < pNode2->num_entries; i++){
				fMbr[0] = pNode2->entries[i].bounces[0] - pNode2->entries[i].fmnd;
				fMbr[1] = pNode2->entries[i].bounces[1] + pNode2->entries[i].fmnd;
				fMbr[2] = pNode2->entries[i].bounces[2] - pNode2->entries[i].fmnd;
				fMbr[3] = pNode2->entries[i].bounces[3] + pNode2->entries[i].fmnd;
				if(mbrIntersect(fMbr1, fMbr)){
					mdrJoin(pNode1, fMbr1, pNode2->entries[i].get_son(), fMbr, pfDd, r);
					pNode2->entries[i].del_son();
				}
			}
		}
	}else{
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(mbrIntersect(pNode1->entries[i].bounces, fMbr2)){
					mdrJoin(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, pfDd, r);
					pNode1->entries[i].del_son();
				}
			}
		}else{
			for(i = 0; i < pNode1->num_entries; i++){
				rtn = pNode1->entries[i].get_son();
				for(j = 0; j < pNode2->num_entries; j++){
					fMbr[0] = pNode2->entries[j].bounces[0] - pNode2->entries[j].fmnd;
					fMbr[1] = pNode2->entries[j].bounces[1] + pNode2->entries[j].fmnd;
					fMbr[2] = pNode2->entries[j].bounces[2] - pNode2->entries[j].fmnd;
					fMbr[3] = pNode2->entries[j].bounces[3] + pNode2->entries[j].fmnd;
					if(mbrIntersect(pNode1->entries[i].bounces, fMbr)){
						mdrJoin(rtn, pNode1->entries[i].bounces, pNode2->entries[j].get_son(), fMbr, pfDd, r);
						pNode2->entries[j].del_son();
					}
				}
				pNode1->entries[i].del_son();
			}
		}
	}
}