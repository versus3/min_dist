// ss_repquery.h
// Sequential Scan Method for the replacement query
//
// Jianzhong Qi
// Date Created : 25/09/2012
// Last Modified: 25/09/2012

#ifndef __SS_REPLACEMENT_QUERY_H
#define __SS_REPLACEMENT_QUERY_H

#include "../point/point.h"
#include "../simpleblockfile/sblockfile.h"

/*
	Sequential scan method
	Input: psbf_c (Client set), psbf_f (facility set), psbf_p (Potential location set)
	Output: Min-dist replacement pair
*/
void ss_repquery(SBlockFile* psbf_c, SBlockFile* psbf_f, SBlockFile* psbf_p);

#endif
