#pragma warning(disable:4996)
// main.cpp
// main function of the min-dist location selection problem
//
// Yuan (Andy), Xue & Jianzhong Qi
// Date Created : 09/04/2010
// Last Modified: 22/06/2013

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "./blockfile/cache.h"
#include "./rtree/rtree.h"
#include "./rtree/rtnode.h"
#include "./rtree/entry.h"
#include "./simpleblockfile/sblockfile.h"

#include "./mndrtree/mndrtree.h"
#include "./mndrtree/mndrtnode.h"
#include "./mndrtree/mndentry.h"

#include "./bbrtree/bbrtree.h"
#include "./bbrtree/bbrtnode.h"
#include "./bbrtree/bbentry.h"

#include "stopwatch/stopwatch.h"	// stopwatch to record the time
#include "algorithms/defs.h"		// commonly used definitions

#include "algorithms/ss.h"			// Sequential Search Algorithm
#include "algorithms/nfc.h"			// Nearest Facility Circle Algorithm
#include "algorithms/bb.h"			// Branch and Bound Algorithm
#include "algorithms/vc.h"			// Voronoi Cell Algorithm
#include "algorithms/mnd.h"			// Maximum NFC Distance Algorithm

// For replacement query
#include "algorithms/ss_repquery.h"			// Sequential Search Algorithm
#include "algorithms/bb_repquery.h"			// Branch and Bound Algorithm
#include "algorithms/mnd_repquery.h"		// Maximum NFC Distance Algorithm

typedef struct{
	float fmnd;
	float f2mnd;
} influence_region_t;

// control whether to create or to restore rtree
#define CREATENEW 1

// control whether to perform location selection query or facility replacement query
#define ADDING_QUERY
// #define REPLACEMENT_QUERY

/*
	Reads the input from text files, and stores into array of client, potential location, and facility
	Input: argv (data file names), 
		n_c (|C|), n_cPerPage (number of clients per block), 
		n_p (|P|), n_pPerPage (number of potential locations per block), 
		n_f (|F|), n_fPerPage (number of facilities per block), 
		n_cPerPageRepquery (number of clients per block, for the replacement query)
	Output: threee sequential files containing C, F and P
*/
void createSimpleBlockFiles(char** argv, int& n_c, int& n_cPerPage, int& n_p, int& n_pPerPage, int& n_f, int& n_fPerPage, int& n_cPerPageRepquery);

/*
	Pre-Computations
	Input: cache (caches to be used by rtrees), 
		pfdnn (dnn array), 
		pfdnn (d2nn array), 
		pfIR (Influence region array), 
		sbf_c (simple block file for C), 
		sbf_f (simple block file for F),
		sbf_c (simple block file for C used in replacement query), 
		rC (R-tree for C), 
		rP (R-tree for P), 
		rF (R-tree for F), 
		rnnC (RNN-tree for C), 
		rC_MND (R-tree for C, with an extra field for mnd), 
		rnnC_DNN (Test implementation of RNN-tree for C, with an extra field for dnn), 
		rC_MND_capacity (Another R-tree variant, with two extra fields for each entry: dnn, number of entries), 
		rC_MND_repquery (R-tree for C, with two extra fields for mnd and 2mnd, used in the replacement query), 
		rP_BB_repquery (R-tree for P, with two extra fields for maxDr and 2mnd, used in the replacement query), 
		rF_BB_repquery (R-tree for F, with two extra fields for maxDr and 2mnd, used in the replacement query), 
		argc, argv (dataset file names given as command line parameters)
	Output: 
*/
void dnnPreCom(int n_c, int n_cPerPage, int n_f, int n_fPerPage, int n_cPerPageRepquery, int nBlockSize, 
			float* pfdnn, float* pfd2nn, influence_region_t* pfIR, 
			SBlockFile*& sbf_c, SBlockFile*& sbf_f, SBlockFile*& sbf_c_repquery);

void rtreePreCom(Cache** cache, float* pfdnn, float* pfd2nn, influence_region_t* pfIR, 
			SBlockFile* sbf_c, int nBlockSize, 
			RTree*& rC, RTree*& rP, RTree*& rF, RTree*& rnnC, 
			CMNDRTree*& rC_MND, CMNDRTree*& rnnC_DNN, CBBRTree*& rC_MND_capacity, 
			CMNDRepqueryRTree*& rC_MND_repquery, 
			CMNDRepqueryRTree*& rP_BB_repquery, 
			CMNDRepqueryRTree*& rF_BB_repquery, 
			int argc, char** argv);

/*
	Make a copy for a file
	Input: srcFileName (source file name), dstFileName (destination file name)
	Output: 
*/
void copyFile(const char* srcFileName, const char* dstFileName);

/*
	create an RNN-tree from updating an R-tree
*/
void createRNN(RTNode* n, SBlockFile* sbf_c); 

/*
	create an mnd-rtree
*/
float createMNDRtree(CMNDRTNode *n, float* fParentMbr, float* pfdnn); 

/*
	create an rnn-tree from an mnd-tree with dnn values, for test purpose only
*/
void createRNN_DNN(CMNDRTNode *n, SBlockFile* sbf_c); 

/*
	create a bb-rtree
*/
void createBB(CBBRTNode *n, float& fMaxFDist, int& nDataEntryNum, float* pfdnn); 

/*
	create an mnd-rtree (C) for the replacement query
*/
void createMNDRepqueryRtreeC(CMNDRepqueryRTNode *n, float* fParentMbr, float* pfdnn, float* pfd2nn, float& fmnd, float &f2mnd); 

/*
	create an mnd-rtree (P) for the replacement query
*/
void createBBRepqueryRtreeP(CMNDRepqueryRTNode* pNodeP, CMNDRepqueryRTNode* pNodeC);

/*
	initialise the data entries of an mnd-rtree (P) for the replacement query
*/
void initialLeafEntriesP(CMNDRepqueryRTNode* pNodeP, float* pMbrP, CMNDRepqueryRTNode* pNodeC, float* pMbrC);

/*
	initialise the inner node entries of an mnd-rtree (P) for the replacement query
*/
void initialInnerEntriesP(CMNDRepqueryRTNode* pNodeP, float* fParentMbr, float& fmnd, float& f2mnd);

/*
	create an mnd-rtree (F) for the replacement query
*/
void createBBRepqueryRtreeF(CMNDRepqueryRTNode* pNodeF, influence_region_t* pfIR);

/*
	initialise the data entries of an mnd-rtree (P) for the replacement query
*/
void initialLeafEntriesF(CMNDRepqueryRTNode* pNodeF, influence_region_t* pfIR);

/*
	initialise the inner node entries of an mnd-rtree (P) for the replacement query
*/
void initialInnerEntriesF(CMNDRepqueryRTNode* pNodeF, float* fParentMbr, float& fmnd, float& f2mnd);

/*****************
* Main Function *
*****************/
void main(int argc, char* argv[])
{
	int i;

	// cache of the r-trees, pagefault counter for the cache
	Cache** cache;
	int nPageCounter[PAGE_COUNTER_NUM];
	cache = new Cache*[CACHE_BLOCK_NUM];
	for (i = 0; i < CACHE_BLOCK_NUM; i++) {
		cache[i] = NULL;
	}
	// dataset arrays. Initialised in preCom()
	float* pfdnn = NULL;
	float* pfd2nn = NULL;
	influence_region_t* pfIR = NULL;

	SBlockFile* sbf_c = NULL;			// Client set
	SBlockFile* sbf_c_repquery = NULL;	// Client set for replacement query
	SBlockFile* sbf_p = NULL;			// Potential location set
	SBlockFile* sbf_f = NULL;			// Facility set

	// R tree of clients, candidate locations, facilities, and nearest facility circles
	RTree* rP = NULL, * rF = NULL, * rC = NULL;
	RTree* rnnC = NULL;								// RNN-tree, for nfc method
	CMNDRTree* rC_MND = NULL, *rnnC_DNN = NULL;		// Rdnn-tree, for mnd method
	CBBRTree* rC_MND_capacity = NULL;				// R-tree variant, for branch and bound method
	
	CMNDRepqueryRTree* rC_MND_repquery = NULL;
	CMNDRepqueryRTree* rF_BB_repquery = NULL;
	CMNDRepqueryRTree* rP_BB_repquery = NULL;

	// |C|, |P|, |F|
	int n_c, n_p, n_f;
	// Block capacity for C, P and F
	int n_cPerPage, n_pPerPage, n_fPerPage;
	
	int nBlockSize;

	// Block capacity for C in the replacement query
	int n_cPerPageRepquery;

	/**** Pre-Computation ****/
	// Initialization: read datasets from files and write them to sequential files
	if (argc != 6) {
		// tN: file which stores n_c (|C|), n_p (|P|), and n_f (|F|);
		// tC, tP, tF: file which stores points of clients, potential locations, and facilities
		// blocksize: the size of a block on hard disk
		fprintf(stderr, "Usage: rtree tN.txt tC.txt tP.txt tF.txt [blocksize]\n");
		exit(EXIT_FAILURE);
	}
	nBlockSize = atoi(argv[5]);
	createSimpleBlockFiles(argv, n_c, n_cPerPage, n_p, n_pPerPage, n_f, n_fPerPage, n_cPerPageRepquery);
	printf("n_c=%d, n_p=%d, n_f=%d\n", n_c, n_p, n_f);

	// Initialization: create simple block file pointers to point to the data files
	sbf_c = new SBlockFile(CLIENT_ARRAY_FILE, nBlockSize, sizeof(Client), n_cPerPage, n_c);
	sbf_p = new SBlockFile(POTENTIAL_LOCATION_ARRAY_FILE, nBlockSize, sizeof(PotentialLocation), n_pPerPage, n_p);
	sbf_f = new SBlockFile(FACILITY_ARRAY_FILE, nBlockSize, sizeof(Facility), n_fPerPage, n_f);

	pfdnn = new float[n_c];
	pfd2nn = new float[n_c];
	pfIR = new influence_region_t[n_f];

	dnnPreCom(n_c, n_cPerPage, n_f, n_fPerPage, n_cPerPageRepquery, nBlockSize, pfdnn, pfd2nn, pfIR, sbf_c, sbf_f, sbf_c_repquery);

	rtreePreCom(cache, pfdnn, pfd2nn, pfIR, sbf_c, nBlockSize, 
		rC, rP, rF, rnnC, rC_MND, rnnC_DNN, rC_MND_capacity, 
		rC_MND_repquery, rP_BB_repquery, rF_BB_repquery, 
		argc, argv);

	sbf_c->m_nPageAccesses = 0;
	sbf_p->m_nPageAccesses = 0;
	sbf_f->m_nPageAccesses = 0;

	sbf_c_repquery->m_nPageAccesses = 0;

	/**** Algorithms ****/
#ifdef ADDING_QUERY
/*	
#define SS
#define BB
#define VC
#define NFC
*/
#define MND
#endif

	printf("Query processing...\n");

#ifdef REPLACEMENT_QUERY
// #define SS_REPLACEMENT
// #define MND_REPLACEMENT
// #define MND_REPLACEMENT_WITH_PRUNING
#define BB_REPLACEMENT
#endif

#ifdef SS
	/* Sequential Scan */
	printf("\n|| SS ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	nPageCounter[0] += sbf_c->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	start(); // start stopwatch
	ss(sbf_c, sbf_p);
	stop();

	nPageCounter[1] += sbf_c->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);

#endif

#ifdef BB
	/* Branch and Bound */

	printf("\n|| Branch and Bound ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[0] += cache[i]->page_faults;
		}
	}
	nPageCounter[0] += sbf_c->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	start(); 
	branchBound(rC_MND_capacity, rP);
	stop();

	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[1] += cache[i]->page_faults;
		}
	}
	nPageCounter[1] += sbf_c->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);

#endif	

#ifdef VC
	/* Voronoi Cell */

	printf("\n|| VC ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[0] += cache[i]->page_faults;
		}
	}
	nPageCounter[0] += sbf_c->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	start();
	//voronoiCell(pfdnn, sbf_p, rC, rP, rF);
	voronoiCell1(pfdnn, sbf_p, rC, rP, rF);
	stop();

	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[1] += cache[i]->page_faults;
		}
	}
	nPageCounter[1] += sbf_c->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);
#endif

#ifdef NFC
	/* Nearest Facility Circle */

	printf("\n|| NFC ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[0] += cache[i]->page_faults;
		}
	}
	nPageCounter[0] += sbf_c->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	start(); 
	//nfc(pfdnn, sbf_p, rR, rP);
	// nfc2(pfdnn, rnnC, rP);  // Join Method
	nfc3(pfdnn, rnnC_DNN, rP);  // Join Method on Expanded Client Tree

	stop();

	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[1] += cache[i]->page_faults;
		}
	}
	nPageCounter[1] += sbf_c->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);
#endif

#ifdef MND
	/* Maximum NFC Distance */

	printf("\n|| MND ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[0] += cache[i]->page_faults;
		}
	}
	nPageCounter[0] += sbf_c->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	start(); 
	mdr(rC_MND, rP);  // Join Method on Expanded Client Tree

	stop();

	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[1] += cache[i]->page_faults;
		}
	}
	nPageCounter[1] += sbf_c->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);

#endif	

#ifdef SS_REPLACEMENT
	/* Sequential Scan */
	printf("\n|| SS_REPLACEMENT ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	nPageCounter[0] += sbf_c_repquery->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	start(); // start stopwatch
	ss_repquery(sbf_c_repquery, sbf_f, sbf_p);
	stop();

	nPageCounter[1] += sbf_c_repquery->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);
#endif

#ifdef MND_REPLACEMENT
	/* Maximum NFC Distance */

	printf("\n|| MND_REPLACEMENT ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[0] += cache[i]->page_faults;
		}
	}
	nPageCounter[0] += sbf_c->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	start(); 
	ss_mnd_repquery(rC_MND_repquery, sbf_f, sbf_p);  // Join Method on Expanded Client Tree

	stop();

	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[1] += cache[i]->page_faults;
		}
	}
	nPageCounter[1] += sbf_c->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);

#endif	

#ifdef MND_REPLACEMENT_WITH_PRUNING
	/* Maximum NFC Distance with Pruning */

	printf("\n|| MND_REPLACEMENT_WITH_PRUNING ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[0] += cache[i]->page_faults;
		}
	}
	nPageCounter[0] += sbf_c->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	// start(); 
	mnd_repquery_with_pruning(rC_MND_repquery, rF, rP);  // Join Method on Expanded Client Tree

	// stop();

	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[1] += cache[i]->page_faults;
		}
	}
	nPageCounter[1] += sbf_c->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);

#endif	

#ifdef BB_REPLACEMENT
	/* Branch & Bound for the replacement query */

	printf("\n|| BB_REPLACEMENT ||\n\n");

	nPageCounter[0] = 0;
	nPageCounter[1] = 0;
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[0] += cache[i]->page_faults;
		}
	}
	nPageCounter[0] += sbf_c->m_nPageAccesses;
	nPageCounter[0] += sbf_p->m_nPageAccesses;
	nPageCounter[0] += sbf_f->m_nPageAccesses;

	start(); 
	bb_repquery(rC_MND_repquery, rF_BB_repquery, rP_BB_repquery);  
	stop();

	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			nPageCounter[1] += cache[i]->page_faults;
		}
	}
	nPageCounter[1] += sbf_c->m_nPageAccesses;
	nPageCounter[1] += sbf_p->m_nPageAccesses;
	nPageCounter[1] += sbf_f->m_nPageAccesses;

	nPageCounter[0] = nPageCounter[1] - nPageCounter[0];
	printf("IO No.: %d\n\n", nPageCounter[0]);

#endif	

	printf("\n");

	delete sbf_c;
	delete sbf_p;
	delete sbf_f;
	delete pfdnn;

	delete rC;
	delete rP;
	delete rF;
	delete rnnC;
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		if (cache[i]) {
			delete cache[i];
		}
	}
	delete cache;
}

/*
	Pre-Computation
	Input: cache (caches to be used by rtrees), 
		pfdnn (dnn array), 
		sbf_c (simple block file for C), 
		sbf_f (simple block file for F),
		sbf_c (simple block file for C used in replacement query), 
		rC (R-tree for C), 
		rP (R-tree for P), 
		rF (R-tree for F), 
		rnnC (RNN-tree for C), 
		rC_MND (R-tree for C, with an extra field for mnd), 
		rnnC_DNN (Test implementation of RNN-tree for C, with an extra field for dnn), 
		rC_MND_capacity (Another R-tree variant, with two extra fields for each entry: dnn, number of entries), 
		rC_MND_repquery (R-tree for C, with two extra fields for mnd and 2mnd, used in the replacement query), 
		rP_BB_repquery (R-tree for P, with two extra fields for maxDr and 2mnd, used in the replacement query), 
		rF_BB_repquery (R-tree for F, with two extra fields for maxDr and 2mnd, used in the replacement query), 
		argc, argv (dataset file names given as command line parameters)
	Output: 
*/
void dnnPreCom(int n_c, int n_cPerPage, int n_f, int n_fPerPage, int n_cPerPageRepquery, int nBlockSize, 
			float* pfdnn, float* pfd2nn, influence_region_t* pfIR, 
			SBlockFile*& sbf_c, SBlockFile*& sbf_f, SBlockFile*& sbf_c_repquery) {
	
	// For dnn pre-computation
	Client* c;
	Facility* f;
	FILE* fpC;
	char* szPadding;			// for padding
	int nCRecord, nFRecord;

	int i, j, m, n;
	float temp;

	// For d2nn pre-computation
	int* pNFIDs;
	
	/*
	int* pDnnIDs = new int[n_c];
	memset(pDnnIDs, -1, n_c * sizeof(int));
	map<int, int> mFCCounter;
	*/
	printf("Precomputing dnn...\n");

	// Pre-compute dnn(c_i): nearest facility distance of each client, 
	szPadding = new char[nBlockSize];
	memset(szPadding, 0, nBlockSize);

	c = new Client[n_cPerPage];
	f = new Facility[n_fPerPage];

	if((fpC = fopen(CLIENT_ARRAY_FILE_WITH_DNN, "wb")) == NULL){
		fprintf(stderr, "Can not open %s", CLIENT_ARRAY_FILE_WITH_DNN);
		exit(0);
	}
	nCRecord = nFRecord = 0;
	for (i = 0; i < n_c; i += nCRecord) {
		nCRecord = sbf_c->readBlock(c);

		sbf_f->rewind();
		for(j = 0; j < n_f; j += nFRecord){
			nFRecord = sbf_f->readBlock(f);

			if(j == 0){
				for(m = 0; m < nCRecord; m++){
					c[m].dnn = c[m].dist(f[0]);
					pfdnn[i+m] = c[m].dnn;
				}
			}
			for(m = 0; m < nCRecord; m++){
				for(n = 0; n < nFRecord; n++){
					temp = c[m].dist(f[n]);
					if(temp < c[m].dnn){
						c[m].dnn = temp;
						pfdnn[i+m] = c[m].dnn;
						// pDnnIDs[i+m] = j + n;
					}
				}
			}
		}
		/*
		for(m = 0; m < nCRecord; m++){
			if (mFCCounter.find(pDnnIDs[i+m]) != mFCCounter.end()) {
				mFCCounter[pDnnIDs[i+m]]++;
			} else {
				mFCCounter.insert(make_pair<int, int>(pDnnIDs[i+m], 1));
			}
		}
		*/
		// Write a new sequential file for clients with dnn
		fwrite((void*)&nCRecord, sizeof(int), 1, fpC);
		fwrite((void*)c, sizeof(Client), nCRecord, fpC);
		fwrite((void*)szPadding, sizeof(char), (nBlockSize - sizeof(int) - sizeof(Client) * nCRecord), fpC);
	}
	fclose(fpC);
	delete[] c;

	// printf("Occupied facilities No.: %d\n", mFCCounter.size());

	// delete[] pDnnIDs;
	// Use the client sequential file with dnn in the following experiments
	delete sbf_c;
	sbf_c = new SBlockFile(CLIENT_ARRAY_FILE_WITH_DNN, nBlockSize, sizeof(Client), n_cPerPage, n_c);

	// Precompute dnn&&d2nn(c_i): second nearest facility distance of each client
	printf("Precomputing d2nn...\n");

	pNFIDs = new int[n_cPerPageRepquery];
	for (i = 0; i < n_f; i++) {
		pfIR[i].fmnd = pfIR[i].f2mnd = 0;
	}

	sbf_c_repquery = new SBlockFile(CLIENT_ARRAY_FILE_FOR_D2NN_COMPUTATION, nBlockSize, sizeof(ClientRepQuery), n_cPerPageRepquery, n_c);
	memset(szPadding, 0, nBlockSize);

	ClientRepQuery* c_repquery = new ClientRepQuery[n_cPerPageRepquery];

	if((fpC = fopen(CLIENT_ARRAY_FILE_WITH_D2NN, "wb")) == NULL){
		fprintf(stderr, "Can not open %s", CLIENT_ARRAY_FILE_WITH_D2NN);
		exit(0);
	}
	nCRecord = nFRecord = 0;
	for (i = 0; i < n_c; i += nCRecord) {
		nCRecord = sbf_c_repquery->readBlock(c_repquery);

		sbf_f->rewind();
		for(j = 0; j < n_f; j += nFRecord){
			nFRecord = sbf_f->readBlock(f);

			if(j == 0){
				for(m = 0; m < nCRecord; m++){
					c_repquery[m].dnn = c_repquery[m].dist(f[0]);
					c_repquery[m].d2nn = c_repquery[m].dist(f[1]);
					
					pNFIDs[m] = j+0;

					if (c_repquery[m].dnn > c_repquery[m].d2nn) {
						temp = c_repquery[m].dnn;
						c_repquery[m].dnn = c_repquery[m].d2nn;
						c_repquery[m].d2nn = temp;

						pNFIDs[m] = j+1;
					}
				}
				for(m = 0; m < nCRecord; m++){
					for(n = 2; n < nFRecord; n++){
						temp = c_repquery[m].dist(f[n]);
						if(temp < c_repquery[m].dnn){
							c_repquery[m].d2nn = c_repquery[m].dnn;
							c_repquery[m].dnn = temp;

							pfd2nn[i+m] = c_repquery[m].d2nn;

							pNFIDs[m] = j+n;
						} else if (temp < c_repquery[m].d2nn) {
							c_repquery[m].d2nn = temp;

							pfd2nn[i+m] = c_repquery[m].d2nn;
						}
					}
				}
			} else {
				for(m = 0; m < nCRecord; m++){
					for(n = 0; n < nFRecord; n++){
						temp = c_repquery[m].dist(f[n]);
						if(temp < c_repquery[m].dnn){
							c_repquery[m].d2nn = c_repquery[m].dnn;
							c_repquery[m].dnn = temp;

							pfd2nn[i+m] = c_repquery[m].d2nn;

							pNFIDs[m] = j+n;
						} else if (temp < c_repquery[m].d2nn) {
							c_repquery[m].d2nn = temp;

							pfd2nn[i+m] = c_repquery[m].d2nn;
						}
					}
				}
			}

			for(m = 0; m < nCRecord; m++){
				pfIR[pNFIDs[m]].fmnd += (c_repquery[m].d2nn - c_repquery[m].dnn);
				if (c_repquery[m].d2nn + c_repquery[m].dnn > pfIR[pNFIDs[m]].f2mnd) {
					pfIR[pNFIDs[m]].f2mnd = c_repquery[m].d2nn + c_repquery[m].dnn;
				}
			}
		}

		// Write a new sequential file for clients with dnn
		fwrite((void*)&nCRecord, sizeof(int), 1, fpC);
		fwrite((void*)c_repquery, sizeof(ClientRepQuery), nCRecord, fpC);
		fwrite((void*)szPadding, sizeof(char), (nBlockSize - sizeof(int) - sizeof(ClientRepQuery) * nCRecord), fpC);
	}
	fclose(fpC);

	if (szPadding) {
		delete[] szPadding;
	}
	if (pNFIDs) {
		delete[] pNFIDs;
		pNFIDs = NULL;
	}
	// Use the client sequential file with d2nn in the replacement query experiments
	delete sbf_c_repquery;
	sbf_c_repquery = new SBlockFile(CLIENT_ARRAY_FILE_WITH_D2NN, nBlockSize, sizeof(ClientRepQuery), n_cPerPageRepquery, n_c);

	delete[] f;
	
	delete[] c_repquery;
}

void rtreePreCom(Cache** cache, float* pfdnn, float* pfd2nn, influence_region_t* pfIR, 
			SBlockFile* sbf_c, int nBlockSize, 
			RTree*& rC, RTree*& rP, RTree*& rF, RTree*& rnnC, 
			CMNDRTree*& rC_MND, CMNDRTree*& rnnC_DNN, CBBRTree*& rC_MND_capacity, 
			CMNDRepqueryRTree*& rC_MND_repquery, 
			CMNDRepqueryRTree*& rP_BB_repquery, 
			CMNDRepqueryRTree*& rF_BB_repquery, 
			int argc, char** argv) {
	int i;
	float* pMbr;

	// Initial cach blocks for the R-trees
	for(i = 0; i < CACHE_BLOCK_NUM; i++){
		cache[i] = new Cache(0, nBlockSize);
	}
	// Create or Restore R trees

	printf("Creating/loading R-trees...\n");
	if (CREATENEW) {
#ifdef ADDING_QUERY
		// construct a new rtree using a text dataset file
		rC = new RTree(argv[2], CLIENT_RTREE, nBlockSize, cache[0], DIMENSION);
		rP = new RTree(argv[3], POTENTIAL_LOCATION_RTREE, nBlockSize, cache[1], DIMENSION);
		rF = new RTree(argv[4], FACILITY_RTREE, nBlockSize, cache[2], DIMENSION);

		// construct rnn-tree for C
		delete rC; // delete rC in order to copy the file
		copyFile(CLIENT_RTREE, CLIENT_RNN_TREE);
		rC = new RTree(CLIENT_RTREE, cache[0]);
		rnnC = new RTree(CLIENT_RNN_TREE, cache[3]);

		// The data stored in CLIENT_RNN_TREE is a copy of CLIENT_RTREE. We update its content to obtain CLIENT_RNN_TREE
		rnnC->load_root();
		createRNN(rnnC->root_ptr, sbf_c);
		rnnC->root = rnnC->root_ptr->block;
		rnnC->del_root();

		// construct mnd-rtree for C
		rC_MND = new CMNDRTree(argv[2], CLIENT_MND_TREE, nBlockSize, cache[4], DIMENSION);	// Load the points from argv[2]
		rC_MND->load_root();
		pMbr = rC_MND->root_ptr->get_mbr();
		createMNDRtree(rC_MND->root_ptr, pMbr, pfdnn);		// Update the mnd values
		delete[] pMbr;

		// construct rnn-tree for C (with dnn values), for test p
		delete rC_MND; 
		copyFile(CLIENT_MND_TREE, CLIENT_RNN_TREE_WITH_DNN);
		rC_MND = new CMNDRTree(CLIENT_MND_TREE, cache[4]);
		rnnC_DNN = new CMNDRTree(CLIENT_RNN_TREE_WITH_DNN, cache[5]);

		rnnC_DNN->load_root();
		createRNN_DNN(rnnC_DNN->root_ptr, sbf_c);
		rnnC_DNN->root = rnnC_DNN->root_ptr->block;
		rnnC_DNN->del_root();

		// construct bb-rtree for C
		float fMaxFDist;
		int nDataEntryNum;

		rC_MND_capacity = new CBBRTree(argv[2], CLIENT_BB_TREE, nBlockSize, cache[6], DIMENSION);	
		rC_MND_capacity->load_root();
		createBB(rC_MND_capacity->root_ptr, fMaxFDist, nDataEntryNum, pfdnn);					// Update the maxdnn values and capacity values
		delete rC_MND_capacity;
		rC_MND_capacity = new CBBRTree(CLIENT_BB_TREE, cache[6]);
#endif
		
#ifdef REPLACEMENT_QUERY
		rC = new RTree(argv[2], CLIENT_RTREE, nBlockSize, cache[0], DIMENSION);
		rP = new RTree(argv[3], POTENTIAL_LOCATION_RTREE, nBlockSize, cache[1], DIMENSION);
		rF = new RTree(argv[4], FACILITY_RTREE, nBlockSize, cache[2], DIMENSION);

		// construct mndrepquery-rtree for C
		float fmnd, f2mnd;
		rC_MND_repquery = new CMNDRepqueryRTree(argv[2], CLIENT_MND_REPQUERY_TREE, nBlockSize, cache[7], DIMENSION);				
		rC_MND_repquery->load_root();
		
		pMbr = rC_MND_repquery->root_ptr->get_mbr();
		createMNDRepqueryRtreeC(rC_MND_repquery->root_ptr, pMbr, pfdnn, pfd2nn, fmnd, f2mnd);	// Update the mnd values
		delete[] pMbr;

		delete rC_MND_repquery; 
		rC_MND_repquery = new CMNDRepqueryRTree(CLIENT_MND_REPQUERY_TREE, cache[7]);

		rC_MND_repquery->load_root();

		// construct mndrepquery-rtree for F
		rF_BB_repquery = new CMNDRepqueryRTree(argv[4], FACILITY_BB_REPQUERY_TREE, nBlockSize, cache[9], DIMENSION);				
		rF_BB_repquery->load_root();
		createBBRepqueryRtreeF(rF_BB_repquery->root_ptr, pfIR);

		delete rF_BB_repquery; 
		rF_BB_repquery = new CMNDRepqueryRTree(FACILITY_BB_REPQUERY_TREE, cache[9]);
		// construct mndrepquery-rtree for P
		rP_BB_repquery = new CMNDRepqueryRTree(argv[3], LOCATION_BB_REPQUERY_TREE, nBlockSize, cache[8], DIMENSION);				
		rP_BB_repquery->load_root();
		createBBRepqueryRtreeP(rP_BB_repquery->root_ptr, rC_MND_repquery->root_ptr);
		
		delete rP_BB_repquery; 
		rP_BB_repquery = new CMNDRepqueryRTree(LOCATION_BB_REPQUERY_TREE, cache[8]);
#endif

	} else {
		// restore rtrees from existing ".rtree" files
#ifdef ADDING_QUERY
		rC = new RTree(CLIENT_RTREE, cache[0]);
		rF = new RTree(FACILITY_RTREE, cache[2]);
		rP = new RTree(POTENTIAL_LOCATION_RTREE, cache[1]);
		rnnC = new RTree(CLIENT_RNN_TREE, cache[3]);

		rC_MND = new CMNDRTree(CLIENT_MND_TREE, cache[4]);
		rnnC_DNN = new CMNDRTree(CLIENT_RNN_TREE_WITH_DNN, cache[5]);
#endif

#ifdef REPLACEMENT_QUERY
		rC = new RTree(CLIENT_RTREE, cache[0]);
		rF = new RTree(FACILITY_RTREE, cache[2]);
		rP = new RTree(POTENTIAL_LOCATION_RTREE, cache[1]);


		rC_MND_repquery = new CMNDRepqueryRTree(CLIENT_MND_REPQUERY_TREE, cache[7]);
		rP_BB_repquery = new CMNDRepqueryRTree(LOCATION_BB_REPQUERY_TREE, cache[8]);
		rF_BB_repquery = new CMNDRepqueryRTree(FACILITY_BB_REPQUERY_TREE, cache[9]);
#endif
	}

	if (pfd2nn) {
		delete[] pfd2nn;
	}
	
	if (pfIR) {
		delete[] pfIR;
	}
}

/*
	Reads the input from text files, and stores into array of client, potential location, and facility
	Input: argv (data file names), 
		n_c (|C|), n_cPerPage (number of clients per block), 
		n_p (|P|), n_pPerPage (number of potential locations per block), 
		n_f (|F|), n_fPerPage (number of facilities per block), 
		n_cPerPageRepquery (number of clients per block, for the replacement query)
	Output: threee sequential files containing C, F and P
*/
void createSimpleBlockFiles(char** argv, int& n_c, int& n_cPerPage, int& n_p, int& n_pPerPage, int& n_f, int& n_fPerPage, int& n_cPerPageRepquery) {
	
	printf("Reading dataset files...\n");
	// temprorary variables to hold intermediate data
	FILE* fpIn, * fpOut;		// File pointers for input/output
	char s[MAXLEN];				// character array for file input
	int i, j;					
	float temp;
	
	// Block size
	int nBlockSize;
		
	// Arrays to hold a block of each type of objects
	Client* c;
	PotentialLocation* p;
	Facility* f;

	ClientRepQuery* c_repquery;

	// For block padding
	int nPaddingLen;
	char* szPadding = NULL;

	// Read n_c, n_p, and n_f
	fopen_s(&fpIn, argv[1], "r");
	fgets(s, MAXLEN - 1, fpIn);
	sscanf_s(s, "%d %d %d", &n_c, &n_p, &n_f);

	fclose(fpIn);

	// Compute the number of objects held in a block
	nBlockSize = atoi(argv[5]);
	n_cPerPage = (nBlockSize-sizeof(int)) / sizeof(Client);
	n_pPerPage = (nBlockSize-sizeof(int)) / sizeof(PotentialLocation);
	n_fPerPage = (nBlockSize-sizeof(int)) / sizeof(Facility);

	n_cPerPageRepquery = (nBlockSize-sizeof(int)) / sizeof(ClientRepQuery);

	// Create object arrays to hold a block of each type of objects
	c = new Client[n_cPerPage];
	p = new PotentialLocation[n_pPerPage];
	f = new Facility[n_fPerPage];

	c_repquery = new ClientRepQuery[n_cPerPageRepquery];

	// Read every point in C and write it to a sequentially accessible file (as an array)
	nPaddingLen = nBlockSize - sizeof(int) - n_cPerPage * sizeof(Client);
	if(nPaddingLen > 0){
		if (szPadding) {
			delete[] szPadding; 
			szPadding = NULL;
		}
		szPadding = new char[nPaddingLen];
		memset(szPadding, 0, nPaddingLen); 
	}

	fopen_s(&fpIn, argv[2], "r");
	if((fpOut = fopen(CLIENT_ARRAY_FILE, "wb")) == NULL){
		fprintf(stderr, "Cannot construct %s", CLIENT_ARRAY_FILE);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < n_c; i+=j) {
		j = 0;
		while(j < n_cPerPage){
			fgets(s, MAXLEN - 1, fpIn);
			sscanf_s(s, "%d %f %f %f %f", &(c[j].id), &temp, &(c[j].x),
				&temp, &(c[j].y));
			if(feof(fpIn)){
				break;
			}
			j++;
		}
		fwrite((void*)&(j), sizeof(int), 1, fpOut);
		fwrite((void*)c, sizeof(Client), j, fpOut);

		if(feof(fpIn)){
			nPaddingLen = nBlockSize - sizeof(int) - j * sizeof(Client);
			if(szPadding){
				delete[] szPadding;
			}
			szPadding = new char[nPaddingLen];
			memset(szPadding, 0, nPaddingLen); 
		}

		if(nPaddingLen > 0){
			fwrite((void*)szPadding, sizeof(char), nPaddingLen, fpOut);
		}
	}
	fclose(fpIn);
	fclose(fpOut);

	// Read every point in P
	if(nPaddingLen > 0){
		delete[] szPadding;
		szPadding = NULL;
	}
	nPaddingLen = nBlockSize - sizeof(int) - n_pPerPage * sizeof(PotentialLocation);
	if(nPaddingLen > 0){
		szPadding = new char[nPaddingLen];
		memset(szPadding, 0, nPaddingLen); 
	}

	fopen_s(&fpIn, argv[3], "r");
	if((fpOut = fopen(POTENTIAL_LOCATION_ARRAY_FILE, "wb")) == NULL){
		fprintf(stderr, "Cannot construct %s", POTENTIAL_LOCATION_ARRAY_FILE);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < n_p; i+=j) {
		j = 0;
		while(j < n_pPerPage){
			fgets(s, MAXLEN - 1, fpIn);
			sscanf_s(s, "%d %f %f %f %f", &(p[j].id), &temp, &(p[j].x),
				&temp, &(p[j].y));
			if(feof(fpIn)){
				break;
			}
			j++;
		}
		fwrite((void*)&(j), sizeof(int), 1, fpOut);
		fwrite((void*)p, sizeof(PotentialLocation), j, fpOut);

		if(feof(fpIn)){
			nPaddingLen = nBlockSize - sizeof(int) - j * sizeof(PotentialLocation);
			if(szPadding){
				delete[] szPadding;
			}
			szPadding = new char[nPaddingLen];
			memset(szPadding, 0, nPaddingLen); 
		}

		if(nPaddingLen > 0){
			fwrite((void*)szPadding, sizeof(char), nPaddingLen, fpOut);
		}
	}
	fclose(fpIn);
	fclose(fpOut);

	// Read every point in F
	if(nPaddingLen > 0){
		delete[] szPadding;
		szPadding = NULL;
	}
	nPaddingLen = nBlockSize - sizeof(int) - n_fPerPage * sizeof(Facility);
	if(nPaddingLen > 0){
		szPadding = new char[nPaddingLen];
		memset(szPadding, 0, nPaddingLen); 
	}

	fopen_s(&fpIn, argv[4], "r");
	if((fpOut = fopen(FACILITY_ARRAY_FILE, "wb")) == NULL){
		fprintf(stderr, "Cannot construct %s", FACILITY_ARRAY_FILE);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < n_f; i+=j) {
		j = 0;
		while(j < n_fPerPage){
			fgets(s, MAXLEN - 1, fpIn);
			sscanf_s(s, "%d %f %f %f %f", &(f[j].id), &temp, &(f[j].x),
				&temp, &(f[j].y));
			if(feof(fpIn)){
				break;
			}
			j++;
		}
		fwrite((void*)&(j), sizeof(int), 1, fpOut);
		fwrite((void*)f, sizeof(Facility), j, fpOut);

		if(feof(fpIn)){
			nPaddingLen = nBlockSize - sizeof(int) - j * sizeof(Facility);
			if(szPadding){
				delete[] szPadding;
			}
			szPadding = new char[nPaddingLen];
			memset(szPadding, 0, nPaddingLen); 
		}

		if(nPaddingLen > 0){
			fwrite((void*)szPadding, sizeof(char), nPaddingLen, fpOut);
		}
	}
	fclose(fpIn);
	fclose(fpOut);

	// Read every point in C and write it to a sequentially accessible file (as an array), for the replacement query
	if(nPaddingLen > 0){
		delete[] szPadding; 
		szPadding = NULL;
	}
	nPaddingLen = nBlockSize - sizeof(int) - n_cPerPageRepquery * sizeof(ClientRepQuery);
	if(nPaddingLen > 0){
		szPadding = new char[nPaddingLen];
		memset(szPadding, 0, nPaddingLen); 
	}

	fopen_s(&fpIn, argv[2], "r");
	if((fpOut = fopen(CLIENT_ARRAY_FILE_FOR_D2NN_COMPUTATION, "wb")) == NULL){
		fprintf(stderr, "Cannot construct %s", CLIENT_ARRAY_FILE_FOR_D2NN_COMPUTATION);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < n_c; i+=j) {
		j = 0;
		while(j < n_cPerPageRepquery){
			fgets(s, MAXLEN - 1, fpIn);
			sscanf_s(s, "%d %f %f %f %f", &(c_repquery[j].id), &temp, &(c_repquery[j].x),
				&temp, &(c_repquery[j].y));
			if(feof(fpIn)){
				break;
			}
			j++;
		}
		fwrite((void*)&(j), sizeof(int), 1, fpOut);
		fwrite((void*)c_repquery, sizeof(ClientRepQuery), j, fpOut);

		if(feof(fpIn)){
			nPaddingLen = nBlockSize - sizeof(int) - j * sizeof(ClientRepQuery);
			if(szPadding){
				delete[] szPadding;
			}
			szPadding = new char[nPaddingLen];
			memset(szPadding, 0, nPaddingLen); 
		}

		if(nPaddingLen > 0){
			fwrite((void*)szPadding, sizeof(char), nPaddingLen, fpOut);
		}
	}
	fclose(fpIn);
	fclose(fpOut);

	if(nPaddingLen > 0){
		delete[] szPadding;
	}
	delete[] c;
	delete[] p;
	delete[] f;

	delete[] c_repquery;
}

/*
	Make a copy for a file
	Input: srcFileName (source file name), dstFileName (destination file name)
	Output: 
*/
void copyFile(const char* srcFileName, const char* dstFileName) {
	FILE * fpSrc, *fpDst;

	fopen_s(&fpSrc, srcFileName, "rb");
	fopen_s(&fpDst, dstFileName, "wb");

	char c;
	while (!feof(fpSrc)) {
		fread(&c, sizeof(char), 1, fpSrc);
		fwrite(&c, sizeof(char), 1, fpDst);
	}

	fclose(fpSrc);
	fclose(fpDst);
}

////////////////////////////////////////////////////////////////////////////////////
// For NFC
/*
	create an RNN-tree from updating an R-tree
*/
void createRNN(RTNode* n, SBlockFile* sbf_c) {

	int index = 1;
	float* mbr;
	Client c;

	int i, j;
	if (n->level != 0) {
		// For a non-leaf node, update the MBR to bound all its updated child entries
		for (i = 0; i < n->num_entries; i++) {
			createRNN(n->entries[i].get_son(), sbf_c);
			mbr = n->entries[i].get_son()->get_mbr();
			for (j = 0; j < 2 * DIMENSION; j++) {
				n->entries[i].bounces[j] = mbr[j];
			}
			n->entries[i].get_son()->dirty = true;
			n->entries[i].del_son();
			delete[] mbr;
		}
	} else { 
		// For a leaf node, update the MBR of an entry (client) with the corresponding dnn value of the entry (client)
		for (i = 0; i < n->num_entries; i++) {
			if(sbf_c->getItem(&c, n->entries[i].son - 1) == 0){
				printf("Client %d not found\n", n->entries[i].son - 1);
				exit(0);
			}

			n->entries[i].bounces[0] = c.x - c.dnn;
			n->entries[i].bounces[1] = c.x + c.dnn;
			n->entries[i].bounces[2] = c.y - c.dnn;
			n->entries[i].bounces[3] = c.y + c.dnn;

			n->dirty = true;
		}
	}
	n->dirty = true;
}

////////////////////////////////////////////////////////////////////////////////////
// For MND
/*
	create an mnd-rtree
*/
float createMNDRtree(CMNDRTNode *n, float* fParentMbr, float* pfdnn) {
	int i;
	float fmaxDr = 0;
	float fMbr[4];

	if (n->level != 0) {
		// For a non-leaf node, update the mnd based on its child entries' mnd values
		for(i = 0; i < 4; i++){
			fMbr[i] = fParentMbr[i];
		}
		for (i = 0; i < n->num_entries; i++) {
			n->entries[i].fmnd = createMNDRtree(n->entries[i].get_son(), n->entries[i].bounces, pfdnn);

			if(n->entries[i].bounces[0] - n->entries[i].fmnd < fMbr[0]){
				fMbr[0] = n->entries[i].bounces[0] - n->entries[i].fmnd;
			}
			if(n->entries[i].bounces[1] + n->entries[i].fmnd > fMbr[1]){
				fMbr[1] = n->entries[i].bounces[1] + n->entries[i].fmnd;
			}

			if(n->entries[i].bounces[2] - n->entries[i].fmnd < fMbr[2]){
				fMbr[2] = n->entries[i].bounces[2] - n->entries[i].fmnd;
			}
			if(n->entries[i].bounces[3] + n->entries[i].fmnd > fMbr[3]){
				fMbr[3] = n->entries[i].bounces[3] + n->entries[i].fmnd;
			}
		}

		if(fParentMbr[0] - fMbr[0] > fmaxDr){
			fmaxDr = fParentMbr[0] - fMbr[0];
		}
		if(fMbr[1] - fParentMbr[1] > fmaxDr){
			fmaxDr = fMbr[1] - fParentMbr[1];
		}

		if(fParentMbr[2] - fMbr[2] > fmaxDr){
			fmaxDr = fParentMbr[2] - fMbr[2];
		}
		if(fMbr[3] - fParentMbr[3] > fmaxDr){
			fmaxDr = fMbr[3] - fParentMbr[3];
		}
	} else { 
		// For a leaf node, fmnd <- dnn
		for (i = 0; i < n->num_entries; i++) {
			n->entries[i].fmnd = pfdnn[n->entries[i].son-1];

			if(n->entries[i].fmnd - (n->entries[i].bounces[0] - fParentMbr[0]) > fmaxDr){
				fmaxDr = n->entries[i].fmnd - (n->entries[i].bounces[0] - fParentMbr[0]);
			}
			if(n->entries[i].fmnd - (fParentMbr[1] - n->entries[i].bounces[0]) > fmaxDr){
				fmaxDr = n->entries[i].fmnd - (fParentMbr[1] - n->entries[i].bounces[0]);
			}

			if(n->entries[i].fmnd - (n->entries[i].bounces[2] - fParentMbr[2]) > fmaxDr){
				fmaxDr = n->entries[i].fmnd - (n->entries[i].bounces[2] - fParentMbr[2]);
			}
			if(n->entries[i].fmnd - (fParentMbr[3] - n->entries[i].bounces[2]) > fmaxDr){
				fmaxDr = n->entries[i].fmnd - (fParentMbr[3] - n->entries[i].bounces[2]);
			}
		}
	}
	n->dirty = true;

	return fmaxDr;
}

/*
	create an rnn-tree from an mnd-tree with dnn values, for test purpose only
*/
void createRNN_DNN(CMNDRTNode *n, SBlockFile* sbf_c) {

	int index = 1;
	float *mbr;
	Client c;

	int i, j;
	if (n->level != 0) {
		// For a non-leaf node, update the MBR to bound all its updated child entries
		for (i = 0; i < n->num_entries; i++) {
			createRNN_DNN(n->entries[i].get_son(), sbf_c);
			mbr = n->entries[i].get_son()->get_mbr();
			for (j = 0; j < 2 * DIMENSION; j++)
				n->entries[i].bounces[j] = mbr[j];
			n->entries[i].get_son()->dirty = true;
			n->entries[i].del_son();
			delete[] mbr;
		}
	} else { 
		// For a leaf node, update the MBR of an entry (client) with the corresponding dnn value of the entry (client)
		for (i = 0; i < n->num_entries; i++) {
			if(sbf_c->getItem(&c, n->entries[i].son - 1) == 0){
				printf("Client %d not found\n", n->entries[i].son - 1);
				exit(0);
			}

			n->entries[i].bounces[0] = c.x - c.dnn;
			n->entries[i].bounces[1] = c.x + c.dnn;
			n->entries[i].bounces[2] = c.y - c.dnn;
			n->entries[i].bounces[3] = c.y + c.dnn;

			n->dirty = true;
		}
	}
	n->dirty = true;

}

////////////////////////////////////////////////////////////////////////////////////
// For BB
/*
	create a bb-rtree
*/
void createBB(CBBRTNode *n, float& fMaxFDist, int& nDataEntryNum, float* pfdnn) {
	int i;

	if (n->level != 0) {
		// For a non-leaf node, update fMaxFDist and nDataEntryNum based on its child entries
		for (i = 0; i < n->num_entries; i++) {
			n->entries[i].nDataEntryNum = 0;
			n->entries[i].fMaxFDist = 0;
			createBB(n->entries[i].get_son(), n->entries[i].fMaxFDist, n->entries[i].nDataEntryNum, pfdnn);

			if(n->entries[i].fMaxFDist > fMaxFDist){
				fMaxFDist = n->entries[i].fMaxFDist;
			}
			nDataEntryNum += n->entries[i].nDataEntryNum;
		}
	} else { 
		// For a leaf node, fMaxFDist <- dnn
		for (i = 0; i < n->num_entries; i++) {
			n->entries[i].fMaxFDist = pfdnn[n->entries[i].son-1];
			n->entries[i].nDataEntryNum = 0;
			if(n->entries[i].fMaxFDist > fMaxFDist){
				fMaxFDist = n->entries[i].fMaxFDist;
			}
		}
		nDataEntryNum = n->num_entries;
	}
	n->dirty = true;
}

////////////////////////////////////////////////////////////////////////////////////
// For MND, the replacement query
/*
	create an mnd-rtree (C) for the replacement query
*/
void createMNDRepqueryRtreeC(CMNDRepqueryRTNode *n, float* fParentMbr, float* pfdnn, float* pfd2nn, float& fmnd, float &f2mnd) {
	int i;
	float fMbr[4], f2Mbr[4];

	fmnd = f2mnd = 0;
	if (n->level != 0) {
		// For a non-leaf node, update the mnd & 2mnd based on its child entries' mnd & 2mnd values
		for(i = 0; i < 4; i++){
			fMbr[i] = f2Mbr[i] = fParentMbr[i];
		}
		for (i = 0; i < n->num_entries; i++) {
			createMNDRepqueryRtreeC(n->entries[i].get_son(), n->entries[i].bounces, pfdnn, pfd2nn, n->entries[i].fmnd, n->entries[i].f2mnd);
	
			// mnd
			if(n->entries[i].bounces[0] - n->entries[i].fmnd < fMbr[0]){
				fMbr[0] = n->entries[i].bounces[0] - n->entries[i].fmnd;
			}
			if(n->entries[i].bounces[1] + n->entries[i].fmnd > fMbr[1]){
				fMbr[1] = n->entries[i].bounces[1] + n->entries[i].fmnd;
			}

			if(n->entries[i].bounces[2] - n->entries[i].fmnd < fMbr[2]){
				fMbr[2] = n->entries[i].bounces[2] - n->entries[i].fmnd;
			}
			if(n->entries[i].bounces[3] + n->entries[i].fmnd > fMbr[3]){
				fMbr[3] = n->entries[i].bounces[3] + n->entries[i].fmnd;
			}

			// 2mnd
			if(n->entries[i].bounces[0] - n->entries[i].f2mnd < f2Mbr[0]){
				f2Mbr[0] = n->entries[i].bounces[0] - n->entries[i].f2mnd;
			}
			if(n->entries[i].bounces[1] + n->entries[i].f2mnd > f2Mbr[1]){
				f2Mbr[1] = n->entries[i].bounces[1] + n->entries[i].f2mnd;
			}

			if(n->entries[i].bounces[2] - n->entries[i].f2mnd < f2Mbr[2]){
				f2Mbr[2] = n->entries[i].bounces[2] - n->entries[i].f2mnd;
			}
			if(n->entries[i].bounces[3] + n->entries[i].f2mnd > f2Mbr[3]){
				f2Mbr[3] = n->entries[i].bounces[3] + n->entries[i].f2mnd;
			}
		}

		// mnd
		if(fParentMbr[0] - fMbr[0] > fmnd){
			fmnd = fParentMbr[0] - fMbr[0];
		}
		if(fMbr[1] - fParentMbr[1] > fmnd){
			fmnd = fMbr[1] - fParentMbr[1];
		}

		if(fParentMbr[2] - fMbr[2] > fmnd){
			fmnd = fParentMbr[2] - fMbr[2];
		}
		if(fMbr[3] - fParentMbr[3] > fmnd){
			fmnd = fMbr[3] - fParentMbr[3];
		}

		// 2mnd
		if(fParentMbr[0] - f2Mbr[0] > f2mnd){
			f2mnd = fParentMbr[0] - f2Mbr[0];
		}
		if(f2Mbr[1] - fParentMbr[1] > f2mnd){
			f2mnd = f2Mbr[1] - fParentMbr[1];
		}

		if(fParentMbr[2] - f2Mbr[2] > f2mnd){
			f2mnd = fParentMbr[2] - f2Mbr[2];
		}
		if(f2Mbr[3] - fParentMbr[3] > f2mnd){
			f2mnd = f2Mbr[3] - fParentMbr[3];
		}
	} else { 
		// For a leaf node, fmnd <- dnn
		for (i = 0; i < n->num_entries; i++) {
			n->entries[i].fmnd = pfdnn[n->entries[i].son-1];
			n->entries[i].f2mnd = pfd2nn[n->entries[i].son-1];

			// mnd
			if(n->entries[i].fmnd - (n->entries[i].bounces[0] - fParentMbr[0]) > fmnd){
				fmnd = n->entries[i].fmnd - (n->entries[i].bounces[0] - fParentMbr[0]);
			}
			if(n->entries[i].fmnd - (fParentMbr[1] - n->entries[i].bounces[0]) > fmnd){
				fmnd = n->entries[i].fmnd - (fParentMbr[1] - n->entries[i].bounces[0]);
			}

			if(n->entries[i].fmnd - (n->entries[i].bounces[2] - fParentMbr[2]) > fmnd){
				fmnd = n->entries[i].fmnd - (n->entries[i].bounces[2] - fParentMbr[2]);
			}
			if(n->entries[i].fmnd - (fParentMbr[3] - n->entries[i].bounces[2]) > fmnd){
				fmnd = n->entries[i].fmnd - (fParentMbr[3] - n->entries[i].bounces[2]);
			}

			// 2mnd
			if(n->entries[i].f2mnd - (n->entries[i].bounces[0] - fParentMbr[0]) > f2mnd){
				f2mnd = n->entries[i].f2mnd - (n->entries[i].bounces[0] - fParentMbr[0]);
			}
			if(n->entries[i].f2mnd - (fParentMbr[1] - n->entries[i].bounces[0]) > f2mnd){
				f2mnd = n->entries[i].f2mnd - (fParentMbr[1] - n->entries[i].bounces[0]);
			}

			if(n->entries[i].f2mnd - (n->entries[i].bounces[2] - fParentMbr[2]) > f2mnd){
				f2mnd = n->entries[i].f2mnd - (n->entries[i].bounces[2] - fParentMbr[2]);
			}
			if(n->entries[i].f2mnd - (fParentMbr[3] - n->entries[i].bounces[2]) > f2mnd){
				f2mnd = n->entries[i].f2mnd - (fParentMbr[3] - n->entries[i].bounces[2]);
			}
		}
	}
	n->dirty = true;

}

/*
	create an mnd-rtree (P) for the replacement query
*/
void createBBRepqueryRtreeP(CMNDRepqueryRTNode* pNodeP, CMNDRepqueryRTNode* pNodeC) {
	
	float* pMbrC = pNodeC->get_mbr();
	float* pMbrP = pNodeP->get_mbr();
	float fmnd, f2mnd;

	pMbrC[0] = pMbrC[0] - DATA_DOMAIN;
	pMbrC[1] = pMbrC[1] + DATA_DOMAIN;
	pMbrC[2] = pMbrC[2] - DATA_DOMAIN;
	pMbrC[3] = pMbrC[3] + DATA_DOMAIN;

	initialLeafEntriesP(pNodeP, pMbrP, pNodeC, pMbrC);

	initialInnerEntriesP(pNodeP, pMbrP, fmnd, f2mnd);

	delete[] pMbrC;
	delete[] pMbrP;
}

/*
	initialise the data entries of an mnd-rtree (P) for the replacement query
*/
void initialLeafEntriesP(CMNDRepqueryRTNode* pNodeP, float* pMbrP, CMNDRepqueryRTNode* pNodeC, float* pMbrC) {
	
	int i, j;
	float fMbr[4];
	float fDistCP;
	if (pNodeP->level == 0 && pNodeC->level == 0) {
		for (i = 0; i < pNodeP->num_entries; i++) {
			if (ptInMbr(pNodeP->entries[i].bounces[0], pNodeP->entries[i].bounces[2], pMbrC)) {
				for (j = 0; j < pNodeC->num_entries; j++) {
					fDistCP = pointDist(pNodeP->entries[i].bounces[0], pNodeP->entries[i].bounces[2], pNodeC->entries[j].bounces[0], pNodeC->entries[j].bounces[2]);
					if (fDistCP < pNodeC->entries[j].fmnd) {
						pNodeP->entries[i].fmnd += pNodeC->entries[j].fmnd - fDistCP;
						if (fDistCP > pNodeP->entries[i].f2mnd) {
							pNodeP->entries[i].f2mnd = fDistCP;
						}
					}
				}
			}
		}
		pNodeP->dirty = true;
	} else if (pNodeP->level == 0) {
		for (i = 0; i < pNodeC->num_entries; i++) {
			fMbr[0] = pNodeC->entries[i].bounces[0] - pNodeC->entries[i].fmnd;
			fMbr[1] = pNodeC->entries[i].bounces[1] + pNodeC->entries[i].fmnd;
			fMbr[2] = pNodeC->entries[i].bounces[2] - pNodeC->entries[i].fmnd;
			fMbr[3] = pNodeC->entries[i].bounces[3] + pNodeC->entries[i].fmnd;
			if (mbrIntersect(pMbrP, fMbr)) {
				initialLeafEntriesP (pNodeP, pMbrP, pNodeC->entries[i].get_son(), fMbr);
				pNodeC->entries[i].del_son();
			}
		}	
	} else if (pNodeC->level == 0) {
		for (i = 0; i < pNodeP->num_entries; i++) {
			if (mbrIntersect(pNodeP->entries[i].bounces, pMbrC)) {
				initialLeafEntriesP (pNodeP->entries[i].get_son(), pNodeP->entries[i].bounces, pNodeC, pMbrC);
			}
		}
	} else {
		for (i = 0; i < pNodeP->num_entries; i++) {
			if (mbrIntersect(pNodeP->entries[i].bounces, pMbrC)) {
				for (j = 0; j < pNodeC->num_entries; j++) {
					fMbr[0] = pNodeC->entries[j].bounces[0] - pNodeC->entries[j].fmnd;
					fMbr[1] = pNodeC->entries[j].bounces[1] + pNodeC->entries[j].fmnd;
					fMbr[2] = pNodeC->entries[j].bounces[2] - pNodeC->entries[j].fmnd;
					fMbr[3] = pNodeC->entries[j].bounces[3] + pNodeC->entries[j].fmnd;
					if (mbrIntersect(pNodeP->entries[i].bounces, fMbr)) {
						initialLeafEntriesP (pNodeP->entries[i].get_son(), pNodeP->entries[i].bounces, pNodeC->entries[j].get_son(), fMbr);
						pNodeC->entries[j].del_son();
					}
				}
			}
			pNodeP->entries[i].del_son();
		}
	}
}

/*
	initialise the inner node entries of an mnd-rtree (P) for the replacement query
*/
void initialInnerEntriesP(CMNDRepqueryRTNode* pNodeP, float* fParentMbr, float& fmnd, float& f2mnd) {

	int i;
	float fMbr[4];

	if (pNodeP->level != 0) {
		// For a non-leaf node, update the mnd based on its child entries' mnd values
		for(i = 0; i < 4; i++){
			fMbr[i] = fParentMbr[i];
		}
		for (i = 0; i < pNodeP->num_entries; i++) {
			initialInnerEntriesP(pNodeP->entries[i].get_son(), pNodeP->entries[i].bounces, pNodeP->entries[i].fmnd, pNodeP->entries[i].f2mnd);
			pNodeP->entries[i].del_son();

			if (pNodeP->entries[i].fmnd > fmnd) {
				fmnd = pNodeP->entries[i].fmnd;
			}

			if(pNodeP->entries[i].bounces[0] - pNodeP->entries[i].f2mnd < fMbr[0]){
				fMbr[0] = pNodeP->entries[i].bounces[0] - pNodeP->entries[i].f2mnd;
			}
			if(pNodeP->entries[i].bounces[1] + pNodeP->entries[i].f2mnd > fMbr[1]){
				fMbr[1] = pNodeP->entries[i].bounces[1] + pNodeP->entries[i].f2mnd;
			}

			if(pNodeP->entries[i].bounces[2] - pNodeP->entries[i].f2mnd < fMbr[2]){
				fMbr[2] = pNodeP->entries[i].bounces[2] - pNodeP->entries[i].f2mnd;
			}
			if(pNodeP->entries[i].bounces[3] + pNodeP->entries[i].f2mnd > fMbr[3]){
				fMbr[3] = pNodeP->entries[i].bounces[3] + pNodeP->entries[i].f2mnd;
			}
		}

		if(fParentMbr[0] - fMbr[0] > f2mnd){
			f2mnd = fParentMbr[0] - fMbr[0];
		}
		if(fMbr[1] - fParentMbr[1] > f2mnd){
			f2mnd = fMbr[1] - fParentMbr[1];
		}

		if(fParentMbr[2] - fMbr[2] > f2mnd){
			f2mnd = fParentMbr[2] - fMbr[2];
		}
		if(fMbr[3] - fParentMbr[3] > f2mnd){
			f2mnd = fMbr[3] - fParentMbr[3];
		}
		pNodeP->dirty = true;
	} else { 
		// For a leaf node, fmnd <- dnn
		for (i = 0; i < pNodeP->num_entries; i++) {
			if (pNodeP->entries[i].fmnd > fmnd) {
				fmnd = pNodeP->entries[i].fmnd;
			}

			if(pNodeP->entries[i].f2mnd - (pNodeP->entries[i].bounces[0] - fParentMbr[0]) > f2mnd){
				f2mnd = pNodeP->entries[i].f2mnd - (pNodeP->entries[i].bounces[0] - fParentMbr[0]);
			}
			if(pNodeP->entries[i].f2mnd - (fParentMbr[1] - pNodeP->entries[i].bounces[0]) > f2mnd){
				f2mnd = pNodeP->entries[i].f2mnd - (fParentMbr[1] - pNodeP->entries[i].bounces[0]);
			}

			if(pNodeP->entries[i].f2mnd - (pNodeP->entries[i].bounces[2] - fParentMbr[2]) > f2mnd){
				f2mnd = pNodeP->entries[i].f2mnd - (pNodeP->entries[i].bounces[2] - fParentMbr[2]);
			}
			if(pNodeP->entries[i].f2mnd - (fParentMbr[3] - pNodeP->entries[i].bounces[2]) > f2mnd){
				f2mnd = pNodeP->entries[i].f2mnd - (fParentMbr[3] - pNodeP->entries[i].bounces[2]);
			}
		}
	}
}

/*
	create an mnd-rtree (P) for the replacement query
*/
void createBBRepqueryRtreeF(CMNDRepqueryRTNode* pNodeF, influence_region_t* pfIR) {
	
	float* pMbrF = pNodeF->get_mbr();
	float fmnd, f2mnd;

	initialLeafEntriesF(pNodeF, pfIR);

	initialInnerEntriesF(pNodeF, pMbrF, fmnd, f2mnd);

	delete[] pMbrF;
}

/*
	initialise the data entries of an mnd-rtree (F) for the replacement query
*/
void initialLeafEntriesF(CMNDRepqueryRTNode* pNodeF, influence_region_t* pfIR) {
	
	int i;
	if (pNodeF->level == 0) {
		for (i = 0; i < pNodeF->num_entries; i++) {
			pNodeF->entries[i].fmnd = pfIR[pNodeF->entries[i].son-1].fmnd;
			pNodeF->entries[i].f2mnd = pfIR[pNodeF->entries[i].son-1].f2mnd;
		}
		pNodeF->dirty = true;
	} else {
		for (i = 0; i < pNodeF->num_entries; i++) {
			initialLeafEntriesF(pNodeF->entries[i].get_son(), pfIR);
		}
	}
}

/*
	initialise the inner node entries of an mnd-rtree (F) for the replacement query
*/
void initialInnerEntriesF(CMNDRepqueryRTNode* pNodeF, float* fParentMbr, float& fmnd, float& f2mnd) {

	int i;
	float fMbr[4];
	fmnd = FLT_MAX;
	if (pNodeF->level != 0) {
		// For a non-leaf node, update the mnd based on its child entries' mnd values
		for(i = 0; i < 4; i++){
			fMbr[i] = fParentMbr[i];
		}
		for (i = 0; i < pNodeF->num_entries; i++) {
			initialInnerEntriesP(pNodeF->entries[i].get_son(), pNodeF->entries[i].bounces, pNodeF->entries[i].fmnd, pNodeF->entries[i].f2mnd);
			pNodeF->entries[i].del_son();

			if (pNodeF->entries[i].fmnd < fmnd) {
				fmnd = pNodeF->entries[i].fmnd;
			}

			if(pNodeF->entries[i].bounces[0] - pNodeF->entries[i].f2mnd < fMbr[0]){
				fMbr[0] = pNodeF->entries[i].bounces[0] - pNodeF->entries[i].f2mnd;
			}
			if(pNodeF->entries[i].bounces[1] + pNodeF->entries[i].f2mnd > fMbr[1]){
				fMbr[1] = pNodeF->entries[i].bounces[1] + pNodeF->entries[i].f2mnd;
			}

			if(pNodeF->entries[i].bounces[2] - pNodeF->entries[i].f2mnd < fMbr[2]){
				fMbr[2] = pNodeF->entries[i].bounces[2] - pNodeF->entries[i].f2mnd;
			}
			if(pNodeF->entries[i].bounces[3] + pNodeF->entries[i].f2mnd > fMbr[3]){
				fMbr[3] = pNodeF->entries[i].bounces[3] + pNodeF->entries[i].f2mnd;
			}
		}

		if(fParentMbr[0] - fMbr[0] > f2mnd){
			f2mnd = fParentMbr[0] - fMbr[0];
		}
		if(fMbr[1] - fParentMbr[1] > f2mnd){
			f2mnd = fMbr[1] - fParentMbr[1];
		}

		if(fParentMbr[2] - fMbr[2] > f2mnd){
			f2mnd = fParentMbr[2] - fMbr[2];
		}
		if(fMbr[3] - fParentMbr[3] > f2mnd){
			f2mnd = fMbr[3] - fParentMbr[3];
		}
		pNodeF->dirty = true;
	} else { 
		// For a leaf node, fmnd <- dnn
		for (i = 0; i < pNodeF->num_entries; i++) {
			if (pNodeF->entries[i].fmnd < fmnd) {
				fmnd = pNodeF->entries[i].fmnd;
			}

			if(pNodeF->entries[i].f2mnd - (pNodeF->entries[i].bounces[0] - fParentMbr[0]) > f2mnd){
				f2mnd = pNodeF->entries[i].f2mnd - (pNodeF->entries[i].bounces[0] - fParentMbr[0]);
			}
			if(pNodeF->entries[i].f2mnd - (fParentMbr[1] - pNodeF->entries[i].bounces[0]) > f2mnd){
				f2mnd = pNodeF->entries[i].f2mnd - (fParentMbr[1] - pNodeF->entries[i].bounces[0]);
			}

			if(pNodeF->entries[i].f2mnd - (pNodeF->entries[i].bounces[2] - fParentMbr[2]) > f2mnd){
				f2mnd = pNodeF->entries[i].f2mnd - (pNodeF->entries[i].bounces[2] - fParentMbr[2]);
			}
			if(pNodeF->entries[i].f2mnd - (fParentMbr[3] - pNodeF->entries[i].bounces[2]) > f2mnd){
				f2mnd = pNodeF->entries[i].f2mnd - (fParentMbr[3] - pNodeF->entries[i].bounces[2]);
			}
		}
	}
}