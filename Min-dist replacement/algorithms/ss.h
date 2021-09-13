// ss.h
// Sequential Scan Method
//
// Yuan (Andy), Xue && Jianzhong Qi
// Date Created : 09/04/2010
// Last Modified: 24/09/2010

#ifndef __SS_H
#define __SS_H

#include "../point/point.h"
#include "../simpleblockfile/sblockfile.h"

/*
	Sequential scan method
	Input: psbf_c (Client set), psbf_p (Potential location set)
	Output: Min-dist point
*/
void ss(SBlockFile* psbf_c, SBlockFile* psbf_p);

#endif
