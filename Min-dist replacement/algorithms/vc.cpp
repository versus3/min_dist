// vc.cpp
// Voronoi Cell Algorithm
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 10/06/2010
//		Last Modified: 10/06/2010

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../point/point.h"
#include "../rtree/rtree.h"
#include "../rtree/rtnode.h"
#include "../rtree/entry.h"

#include "vc.h"

// calculate distance deduction for p
float dd(float* pfdnn, PotentialLocation* p, VoronoiCell* pVC, RTNode* pNodeC, float* fCMbr)
{
	int i;
	float fmnd = 0, fDist;
	queue<Entry*> Q_CEntry;
	Entry* pEntryHead;
	RTNode* pNode;
	vector<Entry*> vNodeList;

	if(!mbrIntersect(pVC->m_fMbr, fCMbr)){
		return fmnd;
	}
	
	if(pNodeC->level == 0){
		for(i = 0; i < pNodeC->num_entries; i++){
			fDist = p->dist(Point(pNodeC->entries[i].bounces[0], pNodeC->entries[i].bounces[2]));
			if(fDist < pfdnn[pNodeC->entries[i].son-1]){
				fmnd += pfdnn[pNodeC->entries[i].son-1] - fDist;
			}
		}
	}else{
		for(i = 0; i < pNodeC->num_entries; i++){
			Q_CEntry.push(&pNodeC->entries[i]);
		}
	}

	while(!Q_CEntry.empty()){
		pEntryHead = Q_CEntry.front();
		Q_CEntry.pop();

		if(!mbrIntersect(pVC->m_fMbr, pEntryHead->bounces)){
			continue;
		}

		pNode = pEntryHead->get_son();
		vNodeList.push_back(pEntryHead);
		if(pNode->level == 0){
			for(i = 0; i < pNode->num_entries; i++){
				
				fDist = p->dist(Point(pNode->entries[i].bounces[0], pNode->entries[i].bounces[2]));
				if(fDist < pfdnn[pNode->entries[i].son-1]){
					fmnd += pfdnn[pNode->entries[i].son-1] - fDist;
				}
			}
		}else{
			for(i = 0; i < pNode->num_entries; i++){
				Q_CEntry.push(&pNode->entries[i]);
			}
		}
	}
	
	
	while(vNodeList.size()){
		pEntryHead = vNodeList.back();
		vNodeList.pop_back();
		pEntryHead->del_son();
	}
	
	return fmnd;
}

// find 4 entries in F that are the nearest neighbors of p in the 4 quadrants
int findNNinF(PotentialLocation* p, RTNode* pNodeF, float pNNinF[4][2], int& nSpaceOverhead)
{
	int i;
	float fx, fy, fDist;
	double fMinDist[4] = { 2147483647.0, 2147483647.0, 2147483647.0, 2147483647.0 };

	priority_queue<QfNode*, vector<QfNode*>, cmp> PQ_FNodes;
	
	QfNode* pHead, * pNew;
	RTNode* pNode;

	if(pNodeF->level == 0){
		pNode = pNodeF;
		for(i = 0; i < pNode->num_entries; i++){
			fx = pNode->entries[i].bounces[0];
			fy = pNode->entries[i].bounces[2];
			if(fx == p->x && fy == p->y){
				pNNinF[0][0] = pNNinF[1][0] = pNNinF[2][0] = pNNinF[3][0] = -1;
				return 0;
			}else{
				fDist = p->dist(Point(fx, fy));
				// Entry in quadrant 1
				if(p->x < fx && p->y <= fy && fDist < fMinDist[0]){
					fMinDist[0] = fDist;
					pNNinF[0][0] = (pNode->entries[i]).bounces[0];
					pNNinF[0][1] = (pNode->entries[i]).bounces[2];
				}
				// Entry in quadrant 2
				if(p->x >= fx && p->y < fy && fDist < fMinDist[1]){
					fMinDist[1] = fDist;
					pNNinF[1][0] = (pNode->entries[i]).bounces[0];
					pNNinF[1][1] = (pNode->entries[i]).bounces[2];
				}
				// Entry in quadrant 3
				if(p->x > fx && p->y >= fy && fDist < fMinDist[2]){
					fMinDist[2] = fDist;
					pNNinF[2][0] = (pNode->entries[i]).bounces[0];
					pNNinF[2][1] = (pNode->entries[i]).bounces[2];
				}
				// Entry in quadrant 4
				if(p->x <= fx && p->y > fy && fDist < fMinDist[3]){
					fMinDist[3] = fDist;
					pNNinF[3][0] = (pNode->entries[i]).bounces[0];
					pNNinF[3][1] = (pNode->entries[i]).bounces[2];
				}
			}
		}
		return 1;
	}

	for(i = 0; i < pNodeF->num_entries; i++){
		pNew = new QfNode;
		pNew->nBlock = pNodeF->entries[i].son;
		pNew->fMinDist = ptMbrDist(p, pNodeF->entries[i].bounces);
		PQ_FNodes.push(pNew);
	}
	
	while(!PQ_FNodes.empty()){
		if(PQ_FNodes.size() > (unsigned int) nSpaceOverhead){
			nSpaceOverhead = PQ_FNodes.size();
		}

		pHead = PQ_FNodes.top();
		PQ_FNodes.pop();

		if(pHead->fMinDist >= fMinDist[0] && pHead->fMinDist >= fMinDist[1] && pHead->fMinDist >= fMinDist[2] && pHead->fMinDist >= fMinDist[3]){
			delete pHead;
			while(!PQ_FNodes.empty()){
				pHead = PQ_FNodes.top();
				PQ_FNodes.pop();
				delete pHead;
			}
			break;
		}

		pNode = new RTNode(pNodeF->my_tree, pHead->nBlock);
		if(pNode->level == 0){
			for(i = 0; i < pNode->num_entries; i++){
				fx = pNode->entries[i].bounces[0];
				fy = pNode->entries[i].bounces[2];
				if(fx == p->x && fy == p->y){
					pNNinF[0][0] = pNNinF[1][0] = pNNinF[2][0] = pNNinF[3][0] = -1;
					
					delete pHead;
					while(!PQ_FNodes.empty()){
						pHead = PQ_FNodes.top();
						PQ_FNodes.pop();
						delete pHead;
					}
					return 0;
				}else{
					fDist = p->dist(Point(fx, fy));
					// Entry in quadrant 1
					if(p->x < fx && p->y <= fy && fDist < fMinDist[0]){
						fMinDist[0] = fDist;
						pNNinF[0][0] = (pNode->entries[i]).bounces[0];
						pNNinF[0][1] = (pNode->entries[i]).bounces[2];
					}
					// Entry in quadrant 2
					if(p->x >= fx && p->y < fy && fDist < fMinDist[1]){
						fMinDist[1] = fDist;
						pNNinF[1][0] = (pNode->entries[i]).bounces[0];
						pNNinF[1][1] = (pNode->entries[i]).bounces[2];
					}
					// Entry in quadrant 3
					if(p->x > fx && p->y >= fy && fDist < fMinDist[2]){
						fMinDist[2] = fDist;
						pNNinF[2][0] = (pNode->entries[i]).bounces[0];
						pNNinF[2][1] = (pNode->entries[i]).bounces[2];
					}
					// Entry in quadrant 4
					if(p->x <= fx && p->y > fy && fDist < fMinDist[3]){
						fMinDist[3] = fDist;
						pNNinF[3][0] = (pNode->entries[i]).bounces[0];
						pNNinF[3][1] = (pNode->entries[i]).bounces[2];
					}
				}
			}
		}else{
			for(i = 0; i < pNode->num_entries; i++){
				fDist = ptMbrDist(p, pNode->entries[i].bounces);

				if(fDist < fMinDist[0] || fDist < fMinDist[1] || fDist < fMinDist[2] || fDist < fMinDist[3]){
					pNew = new QfNode;
					pNew->nBlock = pNode->entries[i].son;
					pNew->fMinDist = fDist;
				
					PQ_FNodes.push(pNew);
				}
			}
		}
		delete pNode;
		delete pHead;
	}

	return 1;
}

// the voronoi cell method for finding the candidate location with max dd
void voronoiCell(float* pfdnn, SBlockFile* psbf_p, RTree * rC, RTree * rP, RTree * rF)
{
	int i, j, k, m, nBestP = -1; // iterators
	float pNNinF[4][2];
	float fmaxDr = 0, fmnd;
	int nSpaceOverhead = 0;

	VoronoiCell* pVC;
	Point* pt[5];
	float* fCMbr;
	PotentialLocation* p = new PotentialLocation[psbf_p->m_nItemPerPage];
	PotentialLocation r;
	r.maxDr = 0;
	int nPRecord = 0;
	
	for(i = 0; i < 5; i++){
		pt[i] = new Point();
	}

	rC->load_root();
	rF->load_root();

	fCMbr = rC->root_ptr->get_mbr();
	// Construct VC(p), find IS(p) and then calculate dd(p)
	psbf_p->rewind();
	for(m = 0; m < psbf_p->m_nItemTotal; m += nPRecord){
		nPRecord = psbf_p->readBlock(p);
		for (i = 0; i < nPRecord; i++) {
			for(j = 0; j < 4; j++){
				pNNinF[j][0] = -1;
			}
			// Find nearest neighbors in F in 4 quadrants for p
			if(findNNinF(&p[i], rF->root_ptr, pNNinF, nSpaceOverhead)){
				pt[0]->x = p[i].x;
				pt[0]->y = p[i].y;
				k = 1;
				for(j = 0; j < 4; j++){
					if(pNNinF[j][0] != -1){
						pt[k]->x = pNNinF[j][0];
						pt[k]->y = pNNinF[j][1];
						k++;
					}
				}

				// construct Voronoi Cell on p
				pVC = new VoronoiCell(pt, k, fCMbr);

				// calculate distance diduction
				fmnd = dd(pfdnn, &p[i], pVC, rC->root_ptr, fCMbr);

				if(fmnd > r.maxDr){
					r.id = p[i].id;
					r.x = p[i].x;
					r.y = p[i].y;
					r.maxDr = fmnd;
				}

				delete pVC;
			}
		}
	}

	/******* Stage 3: Output *******/
	if(r.maxDr > 0){
		printf("dd{p[%d]=(%.2f, %.2f)} = %.2f\nMax SpaceOverhead: %d\n", r.id, r.x, r.y, r.maxDr, nSpaceOverhead);
	}

	for(i = 0; i < 5; i++){
		delete pt[i];
	}
	delete[] p;
	rC->del_root();
	rF->del_root();
}

void dd1(VoronoiCell** pVC, /*int nVcStart, */int nVcEnd, RTNode* pNodeC, float* fCMbr, float* pfdnn, PotentialLocation& r)
{
	int i, j, /*nStart1,*/ nEnd1;
	float fmnd;
	VoronoiCell** pNewVC;

	if(pNodeC->level == 0){
		for(i = /*nVcStart*/0; i < nVcEnd; i++){
			if(mbrIntersect(pVC[i]->m_fMbr, fCMbr)){
				for(j = 0; j < pNodeC->num_entries; j++){
					fmnd = pfdnn[pNodeC->entries[j].son-1] - pVC[i]->m_points[0].dist(pNodeC->entries[j].bounces[0], pNodeC->entries[j].bounces[2]);
					if(fmnd > 0 ){
						pVC[i]->m_fDd += fmnd;
					}
				}
			}
			if(pVC[i]->m_fDd > r.maxDr){
				r.maxDr = pVC[i]->m_fDd;
				r.id = pVC[i]->m_points[0].id;
				r.x = pVC[i]->m_points[0].x;
				r.y = pVC[i]->m_points[0].y;
			}
		}
	}else{
		pNewVC = new VoronoiCell*[nVcEnd];

		for(i = 0; i < pNodeC->num_entries; i++){
			/*
			j = nVcStart;
			while(j < nVcEnd){
				if(j < nVcEnd && !mbrIntersect(pVC[j]->m_fMbr, pNodeC->entries[i].bounces) ){
					j++;
					continue;
				}
				if(j < nVcEnd){
					nStart1 = j;
					j++;
					while(j < nVcEnd && mbrIntersect(pVC[j]->m_fMbr, pNodeC->entries[i].bounces)){
						j++;
					}
					nEnd1 = j;
					
					dd1(pVC, nStart1, nEnd1, pNodeC->entries[i].get_son(), pNodeC->entries[i].bounces, pfdnn, r);
					//dd1(pVC, nStart1, nVcEnd, pNodeC->entries[i].get_son(), pNodeC->entries[i].bounces, pfdnn, r);
					//break;
					j++;
				}
			}
			*/
			nEnd1 = 0; 
			j = 0;
			while(j < nVcEnd){
				if(mbrIntersect(pVC[j]->m_fMbr, pNodeC->entries[i].bounces) ){
					pNewVC[nEnd1++] = pVC[j];
				}
				j++;
			}
			if(nEnd1 > 0){
				dd1(pNewVC, nEnd1, pNodeC->entries[i].get_son(), pNodeC->entries[i].bounces, pfdnn, r);
			}
			pNodeC->entries[i].del_son();
		}

		delete[] pNewVC;
	}
}

// the voronoi cell method for finding the candidate location with max dd
void voronoiCell1(float* pfdnn, SBlockFile* psbf_p, RTree * rC, RTree * rP, RTree * rF)
{
	int i, j, k, m, nVcCounter, nBestP = -1; // iterators
	float pNNinF[4][2];
	int nSpaceOverhead = 0;

	VoronoiCell** pVC = new VoronoiCell*[psbf_p->m_nItemPerPage];

	Point* pt[5];
	float* fCMbr;
	PotentialLocation* p = new PotentialLocation[psbf_p->m_nItemPerPage];
	PotentialLocation r;
	r.maxDr = 0;
	int nPRecord = 0;

	for(i = 0; i < 5; i++){
		pt[i] = new Point();
	}

	rC->load_root();
	rF->load_root();

	fCMbr = rC->root_ptr->get_mbr();
	// Construct VC(p), find IS(p) and then calculate dd(p)
	psbf_p->rewind();
	for(m = 0; m < psbf_p->m_nItemTotal; m += nPRecord){
		nPRecord = psbf_p->readBlock(p);
		nVcCounter = 0;
		for (i = 0; i < nPRecord; i++) {
			for(j = 0; j < 4; j++){
				pNNinF[j][0] = -1;
			}
			// Find nearest neighbors in F in 4 quadrants for p
			if(findNNinF(&p[i], rF->root_ptr, pNNinF, nSpaceOverhead)){
				pt[0]->x = p[i].x;
				pt[0]->y = p[i].y;
				k = 1;
				for(j = 0; j < 4; j++){
					if(pNNinF[j][0] != -1){
						pt[k]->x = pNNinF[j][0];
						pt[k]->y = pNNinF[j][1];
						k++;
					}
				}

				// construct Voronoi Cell on p
				pVC[nVcCounter] = new VoronoiCell(pt, k, fCMbr);
				pVC[nVcCounter++]->m_points[0].id = p[i].id;
			}
		}

		// calculate distance diduction
		if(nVcCounter > 0){
			//dd1(pVC, 0, nVcCounter, rC->root_ptr, fCMbr, pfdnn, r);
			dd1(pVC, nVcCounter, rC->root_ptr, fCMbr, pfdnn, r);
			
			for(j = 0; j < nVcCounter; j++){
				delete pVC[j];
			}
		}
	}


	/******* Stage 3: Output *******/
	if(r.maxDr > 0){
		printf("dd{p[%d]=(%.2f, %.2f)} = %.2f\nMax SpaceOverhead: %d\n", r.id, r.x, r.y, r.maxDr, nSpaceOverhead);
	}

	for(i = 0; i < 5; i++){
		delete pt[i];
	}

	delete[] pVC;
	delete[] p;
	rC->del_root();
	rF->del_root();
}