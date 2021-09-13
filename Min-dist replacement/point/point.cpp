// point.cpp
// define the point: clients, facilities, and candidate locations
//
// Yuan (Andy), Xue
// Date Created : 09/04/2010
// Last Modified: 11/04/2010

#include <math.h>
#include <string.h>
#include "point.h"

/* Point */
Point::Point() {
	this->id = 0;
	this->x = 0;
	this->y = 0;

	this->x1 = 0;
	this->y1 = 0;
}

Point::Point(float x, float y) {
	this->id = 0;
	this->x = x;
	this->y = y;

	this->x1 = 0;
	this->y1 = 0;
}

// Euclidean distance bwteen this point and p2
float Point::dist(Point p2) {
	return sqrt(pow(x - p2.x, 2) + pow(y - p2.y, 2));
}

// Euclidean distance bwteen this point and point (x2, y2)
float Point::dist(float x2, float y2) {
	return sqrt(pow(x - x2, 2) + pow(y - y2, 2));
}

/* Client */
Client::Client() {
	this->dnn = 0;
}

Client::Client(float x, float y, float dnnIn) :
	Point(x, y) {
	this->dnn = dnnIn;
}

/* PotentialLocation */
PotentialLocation::PotentialLocation() {
	this->maxDr = 0;
}

PotentialLocation::PotentialLocation(float x, float y, float maxDrIn) :
	Point(x, y) {
	this->maxDr = maxDrIn;
}

// check whether the potential location is inside MndRegion(c)
bool PotentialLocation::isInMndRegion(Client c) {
	return x > c.x - c.dnn && x < c.x + c.dnn && y > c.y - c.dnn && y < c.y
			+ c.dnn;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// For replacement query
ClientRepQuery::ClientRepQuery() {
	this->dnn = 0;
	this->d2nn = 0;
}

ClientRepQuery::ClientRepQuery(float x, float y, float dnnIn, float d2nnIn) :
	Point(x, y) {
	this->dnn = dnnIn;
	this->d2nn = d2nnIn;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Point and MBR distance functions

/*
	Compute point-point distance
	Input: x1, y1 (point 1), x2, y2 (point 2)
	Output: sqrt((x1 - x2)^2 + (y1 - y2)^2)
*/
float pointDist(float x1, float y1, float x2, float y2) {
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

/*
	Check whether a point is in an MBR
	Input: x, y (point), fMbr (MBR)
	Output: 1 if in, 0 if not
*/
int ptInMbr(float x, float y, float* fMbr) {
	return x <= fMbr[1] && x >= fMbr[0] && y <= fMbr[3] && y >= fMbr[2];
}

/*
	Compute point-MBR distance
	Input: P (point), pMbr (MBR)
	Output: the distance between a point and an MBR
*/
float ptMbrDist(Point* p, float* pMbr) {
	float fMinDist;

	if(p->x <= pMbr[1] && p->x >= pMbr[0] && p->y <= pMbr[3] && p->y >= pMbr[2]){
		fMinDist = 0;
	}else{
		if(p->y <= pMbr[3] && p->y >= pMbr[2]){
			if(p->x >= pMbr[1]){
				fMinDist = p->x - pMbr[1];
			}else{
				fMinDist = pMbr[0] - p->x;
			}
		}else if(p->x <= pMbr[1] && p->x >= pMbr[0]){
			if(p->y >= pMbr[3]){
				fMinDist = p->y - pMbr[3];
			}else{
				fMinDist = pMbr[2] - p->y;
			}
		}else{
			if(p->x > pMbr[1] && p->y > pMbr[3]){
				fMinDist = p->dist(Point(pMbr[1], pMbr[3]));
			}
			if(p->x < pMbr[0] && p->y > pMbr[3]){
				fMinDist = p->dist(Point(pMbr[0], pMbr[3]));
			}
			if(p->x < pMbr[0] && p->y < pMbr[2]){
				fMinDist = p->dist(Point(pMbr[0], pMbr[2]));
			}
			if(p->x > pMbr[1] && p->y < pMbr[2]){
				fMinDist = p->dist(Point(pMbr[1], pMbr[2]));
			}
		}
	}

	return fMinDist;
}

/*
	Compute point-MBR distance
	Input: x1, y1 (point), pMbr (MBR)
	Output: the distance between a point and an MBR
*/
float ptMbrDist(float x, float y, float* pMbr)
{
	float fMinDist;

	if(x <= pMbr[1] && x >= pMbr[0] && y <= pMbr[3] && y >= pMbr[2]){
		fMinDist = 0;
	}else{
		if(y <= pMbr[3] && y >= pMbr[2]){
			if(x >= pMbr[1]){
				fMinDist = x - pMbr[1];
			}else{
				fMinDist = pMbr[0] - x;
			}
		}else if(x <= pMbr[1] && x >= pMbr[0]){
			if(y >= pMbr[3]){
				fMinDist = y - pMbr[3];
			}else{
				fMinDist = pMbr[2] - y;
			}
		}else{
			if(x > pMbr[1] && y > pMbr[3]){
				fMinDist = pointDist(x, y, pMbr[1], pMbr[3]);
			}
			if(x < pMbr[0] && y > pMbr[3]){
				fMinDist = pointDist(x, y, pMbr[0], pMbr[3]);
			}
			if(x < pMbr[0] && y < pMbr[2]){
				fMinDist = pointDist(x, y, pMbr[0], pMbr[2]);
			}
			if(x > pMbr[1] && y < pMbr[2]){
				fMinDist = pointDist(x, y, pMbr[1], pMbr[2]);
			}
		}
	}

	return fMinDist;
}

/*
	Check MBR intersection
	Input: fMbr1 (MBR 1), fMbr2 (MBR 2)
	Output: 1 if intersect, 0 if not
*/
int mbrIntersect(float* fMbr1, float* fMbr2) {
	if(fMbr1[1] < fMbr2[0] || fMbr1[0] > fMbr2[1] || fMbr1[3] < fMbr2[2] || fMbr1[2] > fMbr2[3]){
		return 0;
	}
	return 1;
}

/*
	Get MBR for an array of points
	Input: pPts (point array), nLen (size of array)
	Output: fMbr (the MBR for pPts)
*/
void getMbr(Point* pPts, int nLen, float* fMbr) {

	int i;
    int dimension = 2;
    
	fMbr[0] = fMbr[1] = pPts[0].x;
	fMbr[2] = fMbr[3] = pPts[0].y;

	for (i = 1; i < nLen; i++) {
		if (pPts[i].x < fMbr[0]) {
			fMbr[0] = pPts[i].x;
		}
		if (pPts[i].x > fMbr[1]) {
			fMbr[1] = pPts[i].x;
		}
        if (pPts[i].y < fMbr[2]) {
			fMbr[2] = pPts[i].y;
		}
		if (pPts[i].y > fMbr[3]) {
			fMbr[3] = pPts[i].y;
		}
    }
}