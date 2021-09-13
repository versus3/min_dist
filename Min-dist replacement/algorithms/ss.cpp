// ss.cpp
// Sequential Scan Method
//
// Yuan (Andy), Xue && Jianzhong Qi
// Date Created : 01/03/2010
// Last Modified: 16/06/2010

#include <string.h>
#include <stdio.h>
#include "ss.h"

/* Sequential scan method
   Input: psbf_c (Client set), psbf_p (Potential location set)
   Output: Min-dist point
*/
void ss(SBlockFile* psbf_c, SBlockFile* psbf_p) {
	// store the result - optimal potential location
	PotentialLocation r;

	int i, j, m, n;
	float temp;
	r.maxDr = 0;

	// find dr(p_i) for each p_i by finding clients who consider p_i as their new NN
	Client* c = new Client[psbf_c->m_nItemPerPage];							// An array to hold a block of clients read from the Client set data file
	PotentialLocation* p = new PotentialLocation[psbf_p->m_nItemPerPage];	// An array to hold a block of potential locations read from the Potential Location set data file
	int nCRecord, nPRecord;													// Number of records read from the data file

	nCRecord = nPRecord = 0;
	psbf_p->rewind();
	for (i = 0; i < psbf_p->m_nItemTotal; i += nPRecord) {
		nPRecord = psbf_p->readBlock(p);
		
		for(j = 0; j < nPRecord; j++){
			p[j].maxDr = 0;
		}
		
		psbf_c->rewind();
		for(j = 0; j < psbf_c->m_nItemTotal; j += nCRecord){
			nCRecord = psbf_c->readBlock(c);

			for(m = 0; m < nPRecord; m++){
				for(n = 0; n < nCRecord; n++){
					temp = c[n].dist(p[m]);
					if(temp < c[n].dnn){
						p[m].maxDr += c[n].dnn - temp;
					}
				}
			}
		}

		for(j = 0; j < nPRecord; j++){
			if (p[j].maxDr > r.maxDr) { // update the optimal location
				r.id = p[j].id;
				r.x = p[j].x;
				r.y = p[j].y;
				r.maxDr = p[j].maxDr;
			}
		}
	}
	
	// Output
	printf("dr{p[%d]=(%.2f, %.2f)} = %.2f\n", r.id, r.x, r.y, r.maxDr);

	delete[] c;
	delete[] p;

}