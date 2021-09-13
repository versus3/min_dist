// nfc.cpp
// Nearest Facility Circle Algorithm
//
// Yuan (Andy), Xue && Jianzhong Qi
// Date Created : 01/03/2010
// Last Modified: 16/06/2010

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../point/point.h"
#include "../rtree/rtree.h"
#include "../rtree/rtnode.h"
#include "../rtree/entry.h"

#include "nfc.h"
void nfc(float* pfdnn, SBlockFile* psbf_p, RTree * rR, RTree * rP) 
{
	// store the result - optimal candidate location
	PotentialLocation r;
	PotentialLocation* p = new PotentialLocation[psbf_p->m_nItemPerPage];
	
	IS* pIS;
	IS::iterator itrISMap;

	int nPRecord;
	int nSpaceOverhead = 0;
	int i, j, k; // iterators

	int counter = 0;
	float temp; // temporary variable
	RTNode* pNode;
	list<int>::iterator itrISList;

	/******* Stage 2: Finding Optiomal Candidate Location *******/
	// Stage 2.1: Find PIS(p)
	/*
	for(i = 0; i < psbf_p->m_nItemTotal; i++){
		p[i].id = -1;
		p[i].maxDr = 0;
		//p[i].iS = new SortedLinList;
	}
	*/
	rR->load_root();
	rP->load_root();

	r.maxDr = 0;
	pIS = new IS;
	pisJoin(rP->root_ptr, rP->root_ptr->get_mbr(), rR->root_ptr, rR->root_ptr->get_mbr(), pfdnn, pIS, nSpaceOverhead);
	//pis(rP->root_ptr, rR->root_ptr, p, pfdnn, pIS);

	// Stage 2.2: Find dd(p)
	// sort the array of candiate locations

	nPRecord = 0;
	psbf_p->rewind();
	for (i = 0; i < psbf_p->m_nItemTotal; i += nPRecord) {
		nPRecord = psbf_p->readBlock(p);
		
		for(j = 0; j < nPRecord; j++){
			itrISMap = pIS->find(p[j].id);
			if(itrISMap != pIS->end()){
				p[j].maxDr = itrISMap->second->fmaxDr;
			}else{
				p[j].maxDr = 0;
			}
		}
		
		qsort(p, nPRecord, sizeof(PotentialLocation), cmp);	

		j = 0;
		while (p[j].maxDr > r.maxDr) {
			counter++;
			p[j].maxDr = 0;
			// calculate the actual distance deduction of p[i]
			
			itrISMap = pIS->find(p[j].id);
			itrISList = itrISMap->second->iS.begin();
			//nSpaceOverhead += itrISMap->second->iS.size() + 1;
			while(itrISList != itrISMap->second->iS.end()){
				pNode = new RTNode(rR, *itrISList);

				for(k = 0; k < pNode->num_entries; k++){
					temp = pfdnn[pNode->entries[k].son-1] - 
						p[j].dist(pNode->entries[k].bounces[0] + pfdnn[pNode->entries[k].son-1], pNode->entries[k].bounces[2] + pfdnn[pNode->entries[k].son-1]);
					if (temp > 0) {
						p[j].maxDr += temp;
					}	
				}

				delete pNode;
				itrISList++;
			}
			delete itrISMap->second;
			
			if (p[j].maxDr > r.maxDr) { // update the optimal location
				r.id = p[j].id;
				r.x = p[j].x;
				r.y = p[j].y;
				r.maxDr = p[j].maxDr;
			}
			j++;
		}
	}
	
	
/*
	//qsort(p, psbf_p->m_nItemTotal, sizeof(PotentialLocation), cmp);

	// find dd(p_i) in IS(p_i)
	int counter = 0; // record the number of items in the p[] calculated
	float temp; // temporary variable
	//Linkable *l;
	RTNode* pNode;
	list<int>::iterator itrIS;

	i = 0;
	while (p[i].maxDr > r.maxDr) {
		counter++;
		p[i].maxDr = 0;
		// calculate the actual distance deduction of p[i]
		/*
		itrIS = p[i].iS.begin();
		while(itrIS != p[i].iS.end()){
			pNode = new RTNode(rR, *itrIS);
			
			for(j = 0; j < pNode->num_entries; j++){
				temp = pfdnn[pNode->entries[j].son-1] - 
					p[i].dist(pNode->entries[j].bounces[0] + pfdnn[pNode->entries[j].son-1], pNode->entries[j].bounces[2] + pfdnn[pNode->entries[j].son-1]);
				if (temp > 0) {
					p[i].maxDr += temp;
				}	
			}

			delete pNode;
			itrIS++;
		}
		*/
		/*
		for (j = 0; j < p[i].iS->get_num(); j++) {
			if (j == 0) {
				l = p[i].iS->get_first();
			} else {
				l = p[i].iS->get_next();
			}
			temp = pfdnn[l->son - 1] - p[i].dist(Point (l->bounces[0] + pfdnn[l->son - 1], l->bounces[2] + pfdnn[l->son - 1]));
			if (temp > 0) {
				p[i].maxDr += temp;
			}
		}
		*/
		/*
		if (p[i].maxDr > r.maxDr) { // update the optimal location
			r.id = p[i].id;
			r.x = p[i].x;
			r.y = p[i].y;
			r.maxDr = p[i].maxDr;
		}
		i++;
	}
*/
	/******* Stage 3: Output *******/
	printf("dd{p[%d]=(%.2f, %.2f)} = %.2f\n", r.id, r.x, r.y, r.maxDr);
	printf("Counter: %d out of %d\n", counter, psbf_p->m_nItemTotal);
	printf("Space overhead: %d\n", nSpaceOverhead);

	pIS->clear();
	delete pIS;
	delete[] p;
	rR->del_root();
	rP->del_root();
}

void nfc1(float* pfdnn, SBlockFile* psbf_p, RTree * rR, RTree * rP) 
{
	// store the result - optimal candidate location
	PotentialLocation r;
	
	MAP_DD* pDd = new map<int, float>;
	
	int nSpaceOverhead = 0;

	int counter = 0;

	/******* Stage 2: Finding Optiomal Candidate Location *******/
	rR->load_root();
	rP->load_root();

	r.maxDr = 0;
	pisJoin1(rP->root_ptr, rP->root_ptr->get_mbr(), rR->root_ptr, rR->root_ptr->get_mbr(), pfdnn, pDd, &r);
	nSpaceOverhead = pDd->size() * 2; 

	/******* Stage 3: Output *******/
	printf("dd{p[%d]=(%.2f, %.2f)} = %.2f\n", r.id, r.x, r.y, r.maxDr);
	printf("Counter: %d out of %d\n", counter, psbf_p->m_nItemTotal);
	printf("Space overhead: %d\n", nSpaceOverhead);

	pDd->clear();
	delete pDd;
	rR->del_root();
	rP->del_root();
}

void nfc2(float* pfdnn, RTree * rR, RTree * rP)
{
		// store the result - optimal candidate location
	PotentialLocation r;
	float* fMbr1, * fMbr2;
	float* fmnd = new float[rP->num_of_data];
	memset(fmnd, 0, sizeof(float)*rP->num_of_data);
	/******* Stage 2: Finding Optiomal Candidate Location *******/
	rR->load_root();
	rP->load_root();

	r.maxDr = 0;
	
	fMbr1 = rP->root_ptr->get_mbr();
	fMbr2 = rR->root_ptr->get_mbr();
	if(mbrIntersect(fMbr1, fMbr2)){
		pisJoin2(rP->root_ptr, fMbr1, rR->root_ptr, fMbr2, pfdnn, fmnd, &r);
	}
	
	/******* Stage 3: Output *******/
	if(r.maxDr > 0){
		printf("dd{p[%d]=(%.2f, %.2f)} = %.2f\n", r.id, r.x, r.y, r.maxDr);
	}

	delete[] fmnd;
	rR->del_root();
	rP->del_root();
}

void nfc3(float* pfdnn, CMNDRTree * rR, RTree * rP)
{
		// store the result - optimal candidate location
	PotentialLocation r;
	float* fMbr1, * fMbr2;
	float* fmnd = new float[rP->num_of_data];
	memset(fmnd, 0, sizeof(float)*rP->num_of_data);
	/******* Stage 2: Finding Optiomal Candidate Location *******/
	rR->load_root();
	rP->load_root();

	r.maxDr = 0;
	
	fMbr1 = rP->root_ptr->get_mbr();
	fMbr2 = rR->root_ptr->get_mbr();
	if(mbrIntersect(fMbr1, fMbr2)){
		pisJoin3(rP->root_ptr, fMbr1, rR->root_ptr, fMbr2, pfdnn, fmnd, &r);
	}
	
	/******* Stage 3: Output *******/
	if(r.maxDr > 0){
		printf("dd{p[%d]=(%.2f, %.2f)} = %.2f\n", r.id, r.x, r.y, r.maxDr);
	}

	delete[] fmnd;
	rR->del_root();
	rP->del_root();
}

// used in qsort()
int cmp(const void *a1, const void *a2) {
	PotentialLocation *p1, *p2;
	p1 = (PotentialLocation *) a1;
	p2 = (PotentialLocation *) a2;

	if (p1->maxDr < p2->maxDr) {
		return 1;
	} else if (p1->maxDr > p2->maxDr) {
		return -1;
	} else {
		return 0;
	}
}

// print the info in an rtree (For DEBUG ONLY)
void printRTree(RTNode *n) {
	if (n->level != 0) {
		printf("Level=%d, numEntry=%d\n", n->level, n->num_entries);
		for (int i = 0; i < n->num_entries; i++) {
			printf("%d.%d: [%f - %f] [%f - %f]\n", n->level, i,
					n->entries[i].bounces[0], n->entries[i].bounces[1],
					n->entries[i].bounces[2], n->entries[i].bounces[3]);
			printRTree(n->entries[i].get_son());
		}
	} else { // leaf node
		printf("Level=%d, numEntry=%d\n", n->level, n->num_entries);
		for (int i = 0; i < min(5, n->num_entries); i++) {
			printf("%d.%d: [%f - %f] [%f - %f]\n", n->level, i,
					n->entries[i].bounces[0], n->entries[i].bounces[1],
					n->entries[i].bounces[2], n->entries[i].bounces[3]);
		}
		printf("\n");
	}
}

//join two nodes, find the PIS
//void join(RTNode* pNode1, float* fMbr1, RTNode* pNode2, float* fMbr2, PotentialLocation* p, float* pfdnn, IS* pIS)  
void pisJoin(RTNode* pNode1, float* fMbr1, RTNode* pNode2, float* fMbr2, float* pfdnn, IS* pIS, int& nSpaceOverhead)
{ 
	int i, j;
	//Linkable *copy;
	RTNode* rtn;
	int nFlag;
	
	ISNode* pNodeIS;
	IS::iterator itrIS;

	if(pNode1->level == 0){
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], fMbr2)) {
					nFlag = 0;
					itrIS = pIS->find(pNode1->entries[i].son);
					for(j = 0; j < pNode2->num_entries; j++){
						if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], pNode2->entries[j].bounces)) {
							if(itrIS == pIS->end()){
								pNodeIS = new ISNode;
								pNodeIS->fmaxDr = pfdnn[pNode2->entries[j].son-1];

								pIS->insert(pair<int, ISNode*>(pNode1->entries[i].son, pNodeIS));
								itrIS = pIS->find(pNode1->entries[i].son);
								nSpaceOverhead++;
							}else{
								itrIS->second->fmaxDr += pfdnn[pNode2->entries[j].son-1];
							}
							/*
							if(p[pNode1->entries[i].son-1].id == -1){

								p[pNode1->entries[i].son-1].id = pNode1->entries[i].son;
								p[pNode1->entries[i].son-1].x = pNode1->entries[i].bounces[0];
								p[pNode1->entries[i].son-1].y = pNode1->entries[i].bounces[2];
							}
							p[pNode1->entries[i].son-1].maxDr += pfdnn[pNode2->entries[j].son-1];
							//copy = pNode2->entries[j].gen_Linkable();
							//p[pNode1->entries[i].son-1].iS->insert(copy);
							*/
							nFlag = 1;
						}
					}
					
					if(nFlag){
						//p[pNode1->entries[i].son-1].iS.push_back(pNode2->block);
						itrIS->second->iS.push_back(pNode2->block );
						nSpaceOverhead++;
					}
					
				}
			}
		}else{
			for(i = 0; i < pNode2->num_entries; i++){
				if(mbrIntersect(fMbr1, pNode2->entries[i].bounces)){
					//join(pNode1, fMbr1, pNode2->entries[i].get_son(), pNode2->entries[i].bounces, p, pfdnn, pIS);
					pisJoin(pNode1, fMbr1, pNode2->entries[i].get_son(), pNode2->entries[i].bounces, pfdnn, pIS, nSpaceOverhead);
					pNode2->entries[i].del_son();
				}
			}
		}
	}else{
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(mbrIntersect(pNode1->entries[i].bounces, fMbr2)){
					//join(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, p, pfdnn, pIS);
					pisJoin(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, pfdnn, pIS, nSpaceOverhead);
					pNode1->entries[i].del_son();
				}
			}
		}else{
			for(i = 0; i < pNode1->num_entries; i++){
				rtn = pNode1->entries[i].get_son();
				for(j = 0; j < pNode2->num_entries; j++){
					if(mbrIntersect(pNode1->entries[i].bounces, pNode2->entries[j].bounces)){
						//join(rtn, pNode1->entries[i].bounces, 
						//	pNode2->entries[j].get_son(), pNode2->entries[j].bounces, p, pfdnn, pIS);
						pisJoin(rtn, pNode1->entries[i].bounces, 
							pNode2->entries[j].get_son(), pNode2->entries[j].bounces, pfdnn, pIS, nSpaceOverhead);
						pNode2->entries[j].del_son();
					}
				}
				pNode1->entries[i].del_son();
			}
		}
	}

}

void pisJoin1(RTNode* pNode1, float* fMbr1, RTNode* pNode2, float* fMbr2, float* pfdnn, MAP_DD* pDd, PotentialLocation* r)
{
	int i, j;
	RTNode* rtn;
	float fmnd;

	MAP_DD::iterator itrDd;

	if(pNode1->level == 0){
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], fMbr2)) {
					itrDd = pDd->find(pNode1->entries[i].son);
					for(j = 0; j < pNode2->num_entries; j++){
						if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], pNode2->entries[j].bounces)) {
							fmnd = pfdnn[pNode2->entries[j].son-1] - pointDist(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], 
									pNode2->entries[j].bounces[0] + pfdnn[pNode2->entries[j].son-1], pNode2->entries[j].bounces[2] + pfdnn[pNode2->entries[j].son-1]);
							if(fmnd > 0){
								if(itrDd == pDd->end()){

									pDd->insert(pair<int, float>(pNode1->entries[i].son, fmnd));
									itrDd = pDd->find(pNode1->entries[i].son);
								}else{
									itrDd->second += fmnd;
								}
								if(itrDd->second > r->maxDr){
									r->maxDr = itrDd->second;
									r->id = itrDd->first;
									r->x = pNode1->entries[i].bounces[0];
									r->y = pNode1->entries[i].bounces[2];
								}
							}
						}
					}
				}
			}
		}else{
			for(i = 0; i < pNode2->num_entries; i++){
				if(mbrIntersect(fMbr1, pNode2->entries[i].bounces)){
					//join(pNode1, fMbr1, pNode2->entries[i].get_son(), pNode2->entries[i].bounces, p, pfdnn, pIS);
					pisJoin1(pNode1, fMbr1, pNode2->entries[i].get_son(), pNode2->entries[i].bounces, pfdnn, pDd, r);
					pNode2->entries[i].del_son();
				}
			}
		}
	}else{
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(mbrIntersect(pNode1->entries[i].bounces, fMbr2)){
					//join(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, p, pfdnn, pIS);
					pisJoin1(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, pfdnn, pDd, r);
					pNode1->entries[i].del_son();
				}
			}
		}else{
			for(i = 0; i < pNode1->num_entries; i++){
				rtn = pNode1->entries[i].get_son();
				for(j = 0; j < pNode2->num_entries; j++){
					if(mbrIntersect(pNode1->entries[i].bounces, pNode2->entries[j].bounces)){
						//join(rtn, pNode1->entries[i].bounces, 
						//	pNode2->entries[j].get_son(), pNode2->entries[j].bounces, p, pfdnn, pIS);
						pisJoin1(rtn, pNode1->entries[i].bounces, 
							pNode2->entries[j].get_son(), pNode2->entries[j].bounces, pfdnn, pDd, r);
						pNode2->entries[j].del_son();
					}
				}
				pNode1->entries[i].del_son();
			}
		}
	}
}

void pisJoin2(RTNode* pNode1, float* fMbr1, RTNode* pNode2, float* fMbr2, float* pfdnn, float* pfDd, PotentialLocation* r)
{
	int i, j;
	RTNode* rtn;
	float fmnd;

	if(pNode1->level == 0){
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], fMbr2)) {
					for(j = 0; j < pNode2->num_entries; j++){
						if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], pNode2->entries[j].bounces)) {
							fmnd = pfdnn[pNode2->entries[j].son-1] - pointDist(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], 
									pNode2->entries[j].bounces[0] + pfdnn[pNode2->entries[j].son-1], pNode2->entries[j].bounces[2] + pfdnn[pNode2->entries[j].son-1]);
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
				if(mbrIntersect(fMbr1, pNode2->entries[i].bounces)){
					//join(pNode1, fMbr1, pNode2->entries[i].get_son(), pNode2->entries[i].bounces, p, pfdnn, pIS);
					pisJoin2(pNode1, fMbr1, pNode2->entries[i].get_son(), pNode2->entries[i].bounces, pfdnn, pfDd, r);
					pNode2->entries[i].del_son();
				}
			}
		}
	}else{
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(mbrIntersect(pNode1->entries[i].bounces, fMbr2)){
					//join(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, p, pfdnn, pIS);
					pisJoin2(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, pfdnn, pfDd, r);
					pNode1->entries[i].del_son();
				}
			}
		}else{
			for(i = 0; i < pNode1->num_entries; i++){
				rtn = pNode1->entries[i].get_son();
				for(j = 0; j < pNode2->num_entries; j++){
					if(mbrIntersect(pNode1->entries[i].bounces, pNode2->entries[j].bounces)){
						//join(rtn, pNode1->entries[i].bounces, 
						//	pNode2->entries[j].get_son(), pNode2->entries[j].bounces, p, pfdnn, pIS);
						pisJoin2(rtn, pNode1->entries[i].bounces, 
							pNode2->entries[j].get_son(), pNode2->entries[j].bounces, pfdnn, pfDd, r);
						pNode2->entries[j].del_son();
					}
				}
				pNode1->entries[i].del_son();
			}
		}
	}
}

void pisJoin3(RTNode* pNode1, float* fMbr1, CMNDRTNode* pNode2, float* fMbr2, float* pfdnn, float* pfDd, PotentialLocation* r)
{
	int i, j;
	RTNode* rtn;
	float fmnd;

	if(pNode1->level == 0){
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], fMbr2)) {
					for(j = 0; j < pNode2->num_entries; j++){
						if(ptInMbr(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], pNode2->entries[j].bounces)) {
							fmnd = pfdnn[pNode2->entries[j].son-1] - pointDist(pNode1->entries[i].bounces[0], pNode1->entries[i].bounces[2], 
									pNode2->entries[j].bounces[0] + pfdnn[pNode2->entries[j].son-1], pNode2->entries[j].bounces[2] + pfdnn[pNode2->entries[j].son-1]);
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
				if(mbrIntersect(fMbr1, pNode2->entries[i].bounces)){
					//join(pNode1, fMbr1, pNode2->entries[i].get_son(), pNode2->entries[i].bounces, p, pfdnn, pIS);
					pisJoin3(pNode1, fMbr1, pNode2->entries[i].get_son(), pNode2->entries[i].bounces, pfdnn, pfDd, r);
					pNode2->entries[i].del_son();
				}
			}
		}
	}else{
		if(pNode2->level == 0){
			for(i = 0; i < pNode1->num_entries; i++){
				if(mbrIntersect(pNode1->entries[i].bounces, fMbr2)){
					//join(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, p, pfdnn, pIS);
					pisJoin3(pNode1->entries[i].get_son(), pNode1->entries[i].bounces, pNode2, fMbr2, pfdnn, pfDd, r);
					pNode1->entries[i].del_son();
				}
			}
		}else{
			for(i = 0; i < pNode1->num_entries; i++){
				rtn = pNode1->entries[i].get_son();
				for(j = 0; j < pNode2->num_entries; j++){
					if(mbrIntersect(pNode1->entries[i].bounces, pNode2->entries[j].bounces)){
						//join(rtn, pNode1->entries[i].bounces, 
						//	pNode2->entries[j].get_son(), pNode2->entries[j].bounces, p, pfdnn, pIS);
						pisJoin3(rtn, pNode1->entries[i].bounces, 
							pNode2->entries[j].get_son(), pNode2->entries[j].bounces, pfdnn, pfDd, r);
						pNode2->entries[j].del_son();
					}
				}
				pNode1->entries[i].del_son();
			}
		}
	}
}