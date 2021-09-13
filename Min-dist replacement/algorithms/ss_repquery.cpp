// ss_repquery.cpp
// Sequential Scan Method for the replacement query
//
// Jianzhong Qi
// Date Created : 25/09/2012
// Last Modified: 25/09/2012

#include <string.h>
#include <stdio.h>
#include "ss_repquery.h"

/*
	Sequential scan method
	Input: psbf_c (Client set), psbf_f (facility set), psbf_p (Potential location set)
	Output: Min-dist replacement pair
*/
void ss_repquery(SBlockFile* psbf_c, SBlockFile* psbf_f, SBlockFile* psbf_p) {

	// store the result - optimal replacement pair (facility, potential location)
	opt_facility_location_pair_t optPair;
	optPair.p.maxDr = 0;

	// Temporary variables
	int i, j, k, i1, j1, k1;
	float fDistCP, fDistCF;
	float* fMaxDr = new float[psbf_f->m_nItemPerPage*psbf_p->m_nItemPerPage];

	// find dr(p_i) for each p_i by finding clients who consider p_i as their new NN
	ClientRepQuery* c = new ClientRepQuery[psbf_c->m_nItemPerPage];			// An array to hold a block of clients read from the Client set data file
	Facility* f = new Facility[psbf_f->m_nItemPerPage];						// An array to hold a block of facilities read from the Facility set data file
	PotentialLocation* p = new PotentialLocation[psbf_p->m_nItemPerPage];	// An array to hold a block of potential locations read from the Potential Location set data file
	int nCRecord, nFRecord, nPRecord;										// Number of records read from the data file

	nCRecord = nFRecord = nPRecord = 0;
	psbf_f->rewind();
	for (i = 0; i < psbf_f->m_nItemTotal; i+= nFRecord) {
		nFRecord = psbf_f->readBlock(f);
	
		psbf_p->rewind();
		for (j = 0; j < psbf_p->m_nItemTotal; j += nPRecord) {
			nPRecord = psbf_p->readBlock(p);
			
			memset(fMaxDr, 0, sizeof(float)*psbf_f->m_nItemPerPage*psbf_p->m_nItemPerPage);
			
			psbf_c->rewind();
			for(k = 0; k < psbf_c->m_nItemTotal; k += nCRecord){
				nCRecord = psbf_c->readBlock(c);

				for (i1 = 0; i1 < nFRecord; i1++) {
					for (j1 = 0; j1 < nPRecord; j1++) {
						for (k1 = 0; k1 < nCRecord; k1++) {
							fDistCP = c[k1].dist(p[j1]);
							fDistCF = c[k1].dist(f[i1]);

							if (fDistCP <= c[k1].dnn) {
								fMaxDr[i1*nPRecord+j1] += (c[k1].dnn - fDistCP);
							} else if (fDistCP >= c[k1].d2nn) {
								if (fDistCF > c[k1].dnn) {
									fMaxDr[i1*nPRecord+j1] += 0;
								} else if (fDistCF < c[k1].d2nn) {
									fMaxDr[i1*nPRecord+j1] += (c[k1].dnn - c[k1].d2nn);
								}
							} else {
								if (fDistCF > c[k1].dnn) {
									fMaxDr[i1*nPRecord+j1] += 0;
								} else if (fDistCF < c[k1].d2nn) {
									fMaxDr[i1*nPRecord+j1] += (c[k1].dnn - fDistCP);
								}
							}
						}
						if (fMaxDr[i1*nPRecord+j1] > optPair.p.maxDr) { // update the optimal location
							optPair.f = f[i1];
							optPair.p = p[j1];
							optPair.p.maxDr = fMaxDr[i1*nPRecord+j1];
						}
					}
				}
			}
		}
	}

	// Output
	printf("dr{f[%d]=(%.2f, %.2f), p[%d]=(%.2f, %.2f)} = %.2f\n", 
		optPair.f.id, optPair.f.x, optPair.f.y, optPair.p.id, optPair.p.x, optPair.p.y, optPair.p.maxDr);

	delete[] c;
	delete[] f;
	delete[] p;

	delete[] fMaxDr;
}