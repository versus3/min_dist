#pragma warning(disable:4786)
// bb.cpp
// Branch and Bound Algorithm
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 08/06/2010
//		Last Modified: 24/09/2012

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../point/point.h"
#include "../rtree/rtree.h"
#include "../rtree/rtnode.h"
#include "../rtree/entry.h"

#include "bb.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance metrics
/* 
	Compute the min-dist between two MBRs
	Input: fMbrC, fMbrP (two MBRs)
	output: min-dist(fMbrC, fMbrP)
*/
double minDist(float* fMbrC, float* fMbrP)
{	
	double fMinDist;

	if(fMbrC[0] > fMbrP[1]){
		if(fMbrC[3] < fMbrP[2]){
			fMinDist = pointDist(fMbrC[0], fMbrC[3], fMbrP[1], fMbrP[2]);
		}else if(fMbrC[2] > fMbrP[3]){
			fMinDist = pointDist(fMbrC[0], fMbrC[2], fMbrP[1], fMbrP[3]);
		}else{
			fMinDist = 	fMbrC[0] - fMbrP[1];
		}
	}else if(fMbrC[1] < fMbrP[0]){
		if(fMbrC[3] < fMbrP[2]){
			fMinDist = pointDist(fMbrC[1], fMbrC[3], fMbrP[0], fMbrP[2]);
		}else if(fMbrC[2] > fMbrP[3]){
			fMinDist = pointDist(fMbrC[1], fMbrC[2], fMbrP[0], fMbrP[3]);
		}else{
			fMinDist = 	fMbrP[0] - fMbrC[1];
		}
	}else if(fMbrC[3] < fMbrP[2]){
		fMinDist = 	fMbrP[2] - fMbrC[3];
	}else if(fMbrC[2] > fMbrP[3]){
		fMinDist = 	fMbrC[2] - fMbrP[3];
	}else{
		fMinDist = 0;
	} 

	return fMinDist;
}

/* 
	Compute the min-max-dist between a point and an MBR
	Input: x, y (a point), fMbrP (an MBR)
	output: min-max-dist((x, y), fMbrP)
*/
double minMaxDist(float x, float y, float* fMbrP)
{
	double fMin1, fMin2, fMin3;
	fMin1 = sqrt(pow(x - fMbrP[0], 2) + pow(y - fMbrP[2], 2));
	fMin2 = sqrt(pow(x - fMbrP[0], 2) + pow(y - fMbrP[3], 2));
	
	if(fMin1 > fMin2){
		fMin3 = fMin1;
		fMin1 = fMin2;
		fMin2 = fMin3;
	}

	fMin3 = sqrt(pow(x - fMbrP[1], 2) + pow(y - fMbrP[2], 2));

	if(fMin3 <= fMin1){
		fMin2 = fMin1;
		fMin1 = fMin3;
	}else if(fMin3 < fMin2){
		fMin2 = fMin3;
	}

	fMin3 = sqrt(pow(x - fMbrP[1], 2) + pow(y - fMbrP[3], 2));

	if(fMin3 <= fMin1){
		return fMin1;
	}else if(fMin3 < fMin2){
		return fMin3;
	}

	return fMin2;
}

/* 
	Compute the min-exist-dnn between two MBRs
	Input: fMbrC, fMbrP (two MBRs)
	output: min-exist-dnn(fMbrC, fMbrP)
*/
double minExistDNN(float* fMbrP, float* fMbrC)
{
	double fMinExistDNN = 0, fMinExistDNN1;

	fMinExistDNN1 = minMaxDist(fMbrC/*[0]*/[0], fMbrC/*[0]*/[2], fMbrP); 
	if(fMinExistDNN1 > fMinExistDNN){
		fMinExistDNN = fMinExistDNN1;
	}

	fMinExistDNN1 = minMaxDist(fMbrC/*[0]*/[0], fMbrC/*[0]*/[3], fMbrP); 
	if(fMinExistDNN1 > fMinExistDNN){
		fMinExistDNN = fMinExistDNN1;
	}

	fMinExistDNN1 = minMaxDist(fMbrC/*[0]*/[1], fMbrC/*[0]*/[2], fMbrP); 
	if(fMinExistDNN1 > fMinExistDNN){
		fMinExistDNN = fMinExistDNN1;
	}

	fMinExistDNN1 = minMaxDist(fMbrC/*[0]*/[1], fMbrC/*[0]*/[3], fMbrP); 
	if(fMinExistDNN1 > fMinExistDNN){
		fMinExistDNN = fMinExistDNN1;
	}
	return fMinExistDNN;
}


// Find potential influent set for pi and calculate the lower/upper bound of its dd
double pisEntry(Entry* pCanLoc, float* pfdnn, RTNode* pNodeC, ExtRC* pExtRc)
{
	double fdr = 0;
	double fDist;
	int i;

	for(i = 0; i < pNodeC->num_entries; i++){
		if(pNodeC->level == 0){
			fDist = sqrt(pow(pCanLoc->bounces[0] - pNodeC->entries[i].bounces[0], 2) + pow(pCanLoc->bounces[2] - pNodeC->entries[i].bounces[2], 2));

			fDist -= pfdnn[pNodeC->entries[i].son-1];
			if(fDist < 0){
				fdr += (-fDist);
			}
		}else{
			ExtRC::iterator itr = pExtRc->find(pNodeC->entries[i].son);

			if(/*pCanLoc->bounces[0] < itr->second->fMBR[1][1] && pCanLoc->bounces[0] > itr->second->fMBR[1][0] && 
				pCanLoc->bounces[2] < itr->second->fMBR[1][3] && pCanLoc->bounces[2] > itr->second->fMBR[1][2] &&*/
				minDist(itr->second->fMBR/*[0]*/, pCanLoc->bounces) < itr->second->fMaxFDist){

					fdr += pisEntry(pCanLoc, pfdnn, pNodeC->entries[i].get_son(), pExtRc);
					//pNodeC->entries[i].del_son();
			}
		}
	}

	return fdr;
}

// Find potential influent set for Pi and calculate the lower/upper bound of its dd
void pis(QpNode* pQpNode, float* pfdnn, RTNode* pNodeC, ExtRC::iterator itr, double fMinDist, ExtRC* pExtRc)
{
	int i;
	double fMinDd;
	Entry* pCanLoc = pQpNode->pP;

	if(pNodeC->level == 0){
		pQpNode->fmaxDr += pNodeC->num_entries * (itr->second->fMaxFDist - fMinDist);

		fMinDd = itr->second->fMaxFDist - minExistDNN(pCanLoc->bounces, itr->second->fMBR/*[0]*/);
		if(fMinDd > pQpNode->fMinDd){
			pQpNode->fMinDd = fMinDd;
		}
	}else{
		for(i = 0; i < pNodeC->num_entries; i++){
			itr = pExtRc->find(pNodeC->entries[i].son);
			fMinDist = minDist(itr->second->fMBR/*[0]*/, pCanLoc->bounces );

			if(/*pCanLoc->bounces[0] < itr->second->fMBR[1][1] && pCanLoc->bounces[0] > itr->second->fMBR[1][0] && 
				pCanLoc->bounces[2] < itr->second->fMBR[1][3] && pCanLoc->bounces[2] > itr->second->fMBR[1][2] && */
				fMinDist < itr->second->fMaxFDist){

				pis(pQpNode, pfdnn, pNodeC->entries[i].get_son(), itr, fMinDist, pExtRc); 
				//pNodeC->entries[i].del_son();
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void bbEntryJoin(Entry* pEntry, CBBRTNode* pNode2, float* fMbr2, float fMaxFDist, /*float& fpd, */PotentialLocation* r, float* pfdr)
{
	int i, j;
	float fdr, fMinDist/*, fMinDd*/;
	priority_queue<PQpNode*, vector<PQpNode*>, bbcmp> PQ_PNodes;
	PQpNode* pPQpNode;
	int nFlag;
	
	if(pNode2->level == 0){
		for(j = 0; j < pNode2->num_entries; j++){
			fdr = pNode2->entries[j].fMaxFDist - pointDist(pEntry->bounces[0], pEntry->bounces[2], 
				pNode2->entries[j].bounces[0], pNode2->entries[j].bounces[2]);
			if(fdr > 0){
				pfdr[pEntry->son-1] += fdr;
			}
		}
		if(pfdr[pEntry->son-1] > r->maxDr){
			r->maxDr = pfdr[pEntry->son-1];
			r->id = pEntry->son;
			r->x = pEntry->bounces[0];
			r->y = pEntry->bounces[2];
		}
	}else{
		nFlag = 0;
		for(j = 0; j < pNode2->num_entries; j++){
			fMinDist = ptMbrDist(pEntry->bounces[0], pEntry->bounces[2], pNode2->entries[j].bounces);
			if(fMinDist < pNode2->entries[j].fMaxFDist){
				if(nFlag == 0){
					pPQpNode = new PQpNode;
					pPQpNode->pP = pEntry;
					//pPQpNode->fMinDd = pNode2->entries[j].fMaxFDist - minExistDNN(pNode1->entries[i].bounces, pNode2->entries[j].bounces);
					pPQpNode->fmaxDr = (pNode2->entries[j].fMaxFDist - fMinDist) * pNode2->entries[j].nDataEntryNum;
					pPQpNode->Qc.push_back(&(pNode2->entries[j]));

					nFlag = 1;
				}else{
					/*
					fMinDd = pNode2->entries[j].fMaxFDist - minExistDNN(pNode1->entries[i].bounces, pNode2->entries[j].bounces);
					if(fMinDd > pPQpNode->fMinDd){
					pPQpNode->fMinDd = fMinDd;
					}
					*/
					pPQpNode->fmaxDr += (pNode2->entries[j].fMaxFDist - fMinDist) * pNode2->entries[j].nDataEntryNum;
					pPQpNode->Qc.push_back(&(pNode2->entries[j]));
				}
			}
		}
		if(nFlag){
			if(pPQpNode->fmaxDr > r->maxDr){
				PQ_PNodes.push(pPQpNode);
			}else{
				pPQpNode->Qc.clear();
				delete pPQpNode;
			}
		}

		while(!PQ_PNodes.empty()){
			pPQpNode = PQ_PNodes.top();
			PQ_PNodes.pop();

			if(pPQpNode->fmaxDr > r->maxDr){
				for(i = 0; (unsigned int)i < pPQpNode->Qc.size(); i++){
					bbEntryJoin(pPQpNode->pP, pPQpNode->Qc[i]->get_son(), pPQpNode->Qc[i]->bounces, pPQpNode->Qc[i]->fMaxFDist, /*fpd, */r, pfdr);
					pPQpNode->Qc[i]->del_son();
				}
			}

			pPQpNode->Qc.clear();

			delete pPQpNode;
		}
	}
}

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
void bbJoin(RTNode* pPNode, float* fPMbr, CBBRTNode* pCNode, float* fCMbr, float fMaxFDist, float& fpd, PotentialLocation* r, float* pfdr)
{
	int i, j;
	float fdr, fMinDist, fMinDd;
	//priority_queue<PQpNode*, vector<PQpNode*>, bbcmp> PQ_PNodes;
	queue<PQpNode*> PQ_PNodes;
	PQpNode* pPQpNode;
	int nFlag;
	
	if(pPNode->level == 0){
		if(pCNode->level == 0){
			// Both pPNode and pCNode are leaf nodes, update dr values for the entries in pPNode
			for(i = 0; i < pPNode->num_entries; i++){
				for(j = 0; j < pCNode->num_entries; j++){
					fdr = pCNode->entries[j].fMaxFDist - pointDist(pPNode->entries[i].bounces[0], pPNode->entries[i].bounces[2], 
						pCNode->entries[j].bounces[0], pCNode->entries[j].bounces[2]);
					if(fdr > 0){
						pfdr[pPNode->entries[i].son-1] += fdr;
					}
				}
				if(pfdr[pPNode->entries[i].son-1] > r->maxDr){
					r->maxDr = pfdr[pPNode->entries[i].son-1];
					r->id = pPNode->entries[i].son;
					r->x = pPNode->entries[i].bounces[0];
					r->y = pPNode->entries[i].bounces[2];
				}
			}
			if(r->maxDr > fpd){
				fpd = r->maxDr;
			}
		}else{
			// pPNode is a leaf node, pCNode is not, expand pCNode, join pPNode wtih pCNode's children
			for(i = 0; i < pCNode->num_entries; i++){
				fMinDist = minDist(fPMbr, pCNode->entries[i].bounces);
				if(fMinDist < pCNode->entries[i].fMaxFDist){
					bbJoin(pPNode, fPMbr, pCNode->entries[i].get_son(), pCNode->entries[i].bounces, pCNode->entries[i].fMaxFDist, fpd, r, pfdr);
					pCNode->entries[i].del_son();
				}
			}
		}
	}else{
		if(pCNode->level == 0){
			// pCNode is a leaf node, pPNode is not, expand pPNode, join pCNode wtih pPNode's children
			for(i = 0; i < pPNode->num_entries; i++){
				fMinDist = minDist(fCMbr, pPNode->entries[i].bounces);
				if(fMinDist < fMaxFDist && (fMaxFDist - fMinDist) * pCNode->num_entries > fpd){
					pPQpNode = new PQpNode;
					pPQpNode->pP = &(pPNode->entries[i]);
					pPQpNode->fMinDd = fMaxFDist - minExistDNN(pPNode->entries[i].bounces, fCMbr);
					pPQpNode->fmaxDr = (fMaxFDist - fMinDist) * pCNode->num_entries;
					
					if(pPQpNode->fMinDd > fpd){
						fpd = pPQpNode->fMinDd;
					}
					PQ_PNodes.push(pPQpNode);
				}
			}

			while(!PQ_PNodes.empty()){
				//pPQpNode = PQ_PNodes.top();
				pPQpNode = PQ_PNodes.front();
				PQ_PNodes.pop();

				if(pPQpNode->fmaxDr > fpd/*r->maxDr*/){
					if(pPQpNode->fMinDd > fpd){
						fpd = pPQpNode->fMinDd;
					}
					bbJoin(pPQpNode->pP->get_son(), pPQpNode->pP->bounces, pCNode, fCMbr, fMaxFDist, fpd, r, pfdr);

					pPQpNode->pP->del_son();
				}

				delete pPQpNode;
			}
		}else{
			for(i = 0; i < pPNode->num_entries; i++){
				nFlag = 0;
				for(j = 0; j < pCNode->num_entries; j++){
					fMinDist = minDist(pPNode->entries[i].bounces, pCNode->entries[j].bounces);
					if(fMinDist < pCNode->entries[j].fMaxFDist){
						if(nFlag == 0){
							pPQpNode = new PQpNode;
							pPQpNode->pP = &(pPNode->entries[i]);
							pPQpNode->fMinDd = pCNode->entries[j].fMaxFDist - minExistDNN(pPNode->entries[i].bounces, pCNode->entries[j].bounces);
							pPQpNode->fmaxDr = (pCNode->entries[j].fMaxFDist - fMinDist) * pCNode->entries[j].nDataEntryNum;
							pPQpNode->Qc.push_back(&(pCNode->entries[j]));

							nFlag = 1;
						}else{
							fMinDd = pCNode->entries[j].fMaxFDist - minExistDNN(pPNode->entries[i].bounces, pCNode->entries[j].bounces);
							if(fMinDd > pPQpNode->fMinDd){
								pPQpNode->fMinDd = fMinDd;
							}
							
							pPQpNode->fmaxDr += (pCNode->entries[j].fMaxFDist - fMinDist) * pCNode->entries[j].nDataEntryNum;
							pPQpNode->Qc.push_back(&(pCNode->entries[j]));
						}
					}
				}
				if(nFlag){
					if(pPQpNode->fmaxDr > fpd){
						if(pPQpNode->fMinDd > fpd){
							fpd = pPQpNode->fMinDd;
						}
						PQ_PNodes.push(pPQpNode);
					}else{
						pPQpNode->Qc.clear();
						delete pPQpNode;
					}
				}
			}
			while(!PQ_PNodes.empty()){
				//pPQpNode = PQ_PNodes.top();
				pPQpNode = PQ_PNodes.front();
				PQ_PNodes.pop();

				if(pPQpNode->fmaxDr > fpd/*r->maxDr*/){
					for(i = 0; (unsigned int)i < pPQpNode->Qc.size(); i++){
						bbJoin(pPQpNode->pP->get_son(), pPQpNode->pP->bounces, 
							pPQpNode->Qc[i]->get_son(), pPQpNode->Qc[i]->bounces, pPQpNode->Qc[i]->fMaxFDist, fpd, r, pfdr);
						pPQpNode->Qc[i]->del_son();
					}
					pPQpNode->pP->del_son();
				}

				pPQpNode->Qc.clear();

				delete pPQpNode;
			}
		}
	}
}

/* 
	Branch and bound method
	Input: rC (bb-rtree on C), rp (r-tree on P)
	output: min-dist location
*/
void branchBound(CBBRTree * rC, RTree * rP) {
	PotentialLocation r;
	float* pfdr = new float[rP->num_of_data];
	memset(pfdr, 0, sizeof(float)*rP->num_of_data);

	float* fMbr1, * fMbr2;
	float fpd = 0, fMaxFDist = 0;
	int i;
	/******* Stage 2: Finding the Optiomal Potential Location *******/
	r.maxDr = 0;
	
	rC->load_root();
	rP->load_root();

	fMbr1 = rP->root_ptr->get_mbr();
	fMbr2 = rC->root_ptr->get_mbr();

	float fMinDist = (float)minDist(fMbr1, fMbr2) ;
	
	for(i = 0; i < rC->root_ptr->num_entries; i++){
		if(rC->root_ptr->entries[i].fMaxFDist > fMaxFDist){
			fMaxFDist = rC->root_ptr->entries[i].fMaxFDist;
		}
	}
	if(fMinDist < fMaxFDist){ // If the two root nodes are close enough, then investigate the child nodes
		bbJoin(rP->root_ptr, fMbr1, rC->root_ptr, fMbr2, fMaxFDist, fpd, &r, pfdr);
	}
	
	/******* Stage 3: Output *******/
	if(r.maxDr > 0){
		printf("dr{p[%d]=(%.2f, %.2f)} = %.2f\n", r.id, r.x, r.y, r.maxDr);
	}else{
		printf("No valid solution found.\n");
	}

	delete[] pfdr;
	rC->del_root();
	rP->del_root();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An experimental version, implementing the bb-rtree with an r-tree and the extra information stored in a separated data structure
void branchBound_Ex(float* pfdnn, RTree * rC, /*RTree * rR, */RTree * rP, ExtRC* pExtRc) {
	QpNode* pQpHead = new QpNode, * pQpTail;
	QpNode* itrQp, * itrQp1;
	RTNode* pNode;

	int nPairTotal = 0;
	int nProcessedTotal = 0;

	double fpd = 0, fdr, fmaxDr = 0;
	int i;
	Entry* pBestCand = NULL;
	float* pMbr;

	ExtRC::iterator itr;
	double fMinDist;

	// Add root node of tree P to initiate the Branch and Bound method
	pQpHead->pP = new Entry(2, rP);
	pMbr = rP->root_ptr->get_mbr();
	for(i =0; i < 4; i++){
		pQpHead->pP->bounces[i] = pMbr[i];
	}
	pQpHead->nLevel = rP->root_ptr->level;
	pQpHead->fMinDd = pQpHead->fmaxDr = 0;
	pQpHead->pNext = NULL;

	pQpTail = pQpHead;
	itrQp = pQpHead;

	itr = pExtRc->find(rC->root_ptr->block);
	while(itrQp){
		// Find potential influent set for Pi and calculate the lower/upper bound of its dd
		
		fMinDist = minDist(itr->second->fMBR/*[0]*/, itrQp->pP->bounces);

		if(/*itrQp->pP->bounces[0] < itr->second->fMBR[1][1] && itrQp->pP->bounces[0] > itr->second->fMBR[1][0] && 
			itrQp->pP->bounces[2] < itr->second->fMBR[1][3] && itrQp->pP->bounces[2] > itr->second->fMBR[1][2] && */
			fMinDist < itr->second->fMaxFDist){

			pis(itrQp, pfdnn, rC->root_ptr, itr, fMinDist, pExtRc);
			nProcessedTotal++;

			// If Pi has an upper bound of its dd >= current minDdd, check its sons
			if(itrQp->fmaxDr >= fpd){
				if(itrQp->fMinDd > fpd){
					fpd = itrQp->fMinDd;
				}
				if(itrQp->nLevel == rP->root_ptr->level){
					pNode = rP->root_ptr;
				}else{
					pNode = itrQp->pP->get_son();
				}
				// Pi is a leaf node, decide if any of its son has a better dd than the current solution
				if(itrQp->nLevel == 0){
					for(i = 0; i < pNode->num_entries; i++){
						fMinDist = minDist(itr->second->fMBR/*[0]*/, pNode->entries[i].bounces);

						if(/*pNode->entries[i].bounces[0] < itr->second->fMBR[1][1] && pNode->entries[i].bounces[0] > itr->second->fMBR[1][0] && 
							pNode->entries[i].bounces[2] < itr->second->fMBR[1][3] && pNode->entries[i].bounces[2] > itr->second->fMBR[1][2] &&*/
							fMinDist < itr->second->fMaxFDist){

							fdr = pisEntry(&(pNode->entries[i]), pfdnn, rC->root_ptr, pExtRc);

							if(fdr >= fpd){
								fpd = fdr;
								if(fdr > fmaxDr){
									pBestCand = &(pNode->entries[i]);
									fmaxDr = fdr;
								}
							}
						}
					}
					// Pi is an inner node, add its sons for further consideration
				}else{
					for(i = 0; i < pNode->num_entries; i++){
						itrQp1 = new QpNode;
						itrQp1->pP = &(pNode->entries[i]);
						itrQp1->nLevel = itrQp->nLevel - 1;
						itrQp1->fMinDd = itrQp1->fmaxDr = 0;
						itrQp1->pNext = NULL;

						pQpTail->pNext = itrQp1; 
						pQpTail = itrQp1;

						nPairTotal++;
					}
				}
			}
		
		}
		itrQp1 = itrQp;
		itrQp = itrQp->pNext;
		delete itrQp1;
	}

	// Output result
	if(pBestCand){
		printf("dd{p[%d]=(%.2f, %.2f)} = %.2f\n", pBestCand->son, pBestCand->bounces[0], pBestCand->bounces[2], fmaxDr);
	}else{
		printf("No valid solution found.\n");
	}

	printf("Pair Generated Total: %d, Pair Processed Total: %d\n", nPairTotal, nProcessedTotal);
}

// Clear the fMaxFDist B-Tree
void clearExtRc(ExtRC* pExtRc, RTNode* pNodeC) {
	ExtRC::iterator itr;
	int i;

	itr = pExtRc->find(pNodeC->block);

	if(itr != pExtRc->end()){
		delete itr->second;
		pExtRc->erase(itr);
	}

	if(pNodeC->level > 0){
		for(i = 0; i < pNodeC->num_entries; i++){
			clearExtRc(pExtRc, pNodeC->entries[i].get_son());
			//pNodeC->entries[i].del_son();
		}
	}
}

// Create a B-Tree to store fMaxFDist for each node in Rc
double initExtRc(ExtRC* pExtRc, float* pfdnn, RTNode* pNodeC, float* fCMbr/*, RTNode* pNodeR, float* fRMbr*/) {
	CNode* pExtInfo = new CNode;
	int i;
	double fMaxFDist = 0, fMaxFDist1;

	//pExtInfo->pRRNode = pNodeR;
	for(i = 0; i < 4; i++){
		//pExtInfo->fMBR[0][i] = fCMbr[i];
		//pExtInfo->fMBR[1][i] = fRMbr[i];
		pExtInfo->fMBR[i] = fCMbr[i];
	}
	
	if(pNodeC->level == 0){
		for(i = 0; i < pNodeC->num_entries; i++){
			if(pfdnn[pNodeC->entries[i].son-1] > fMaxFDist){
				fMaxFDist = pfdnn[pNodeC->entries[i].son-1];
			}
		}
	}else{
		for(i = 0; i < pNodeC->num_entries; i++){
			fMaxFDist1 = initExtRc(pExtRc, pfdnn, pNodeC->entries[i].get_son(), pNodeC->entries[i].bounces/*, pNodeR->entries[i].get_son(), pNodeR->entries[i].bounces*/);

			//pNodeC->entries[i].del_son();
			//pNodeR->entries[i].del_son();

			if(fMaxFDist1 > fMaxFDist){
				fMaxFDist = fMaxFDist1;
			}
		}
	}
	pExtInfo->fMaxFDist = fMaxFDist;
	pExtRc->insert(pair<int, CNode*>(pNodeC->block, pExtInfo));

	return fMaxFDist;
}