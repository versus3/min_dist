// point.h
// define the point: clients, facilities, and candidate locations
//
// Yuan (Andy), Xue
// Date Created : 09/04/2010
// Last Modified: 11/04/2010

#ifndef __POINT_H
#define __POINT_H

#include <stdio.h>
#include "../linlist/linlist.h"

class Client;

// define the points: clients, facilities, and potential locations
class Point {
public:
	int id;			// id
	float x;		// x coordinate
	float y;		// y coordinate
	float x1, y1;	// not used, for equi-space occupation (compared with R-tree entry) only
	float dist(Point p2);				// Euclidean distance bwteen this point and p2
	float dist(float x2, float y2);		// Euclidean distance bwteen this point and p2

	Point();
	Point(float x, float y);
};

class Client: public Point {
public:
	float dnn; // distance to its nearest facility

	Client();
	Client(float x, float y, float dnnIn);
};

class PotentialLocation: public Point {
public:
	float maxDr; // upperbound of distace deduction
	
	PotentialLocation();
	PotentialLocation(float x, float y, float maxDrIn);

	bool isInMndRegion(Client c); // check whether the potential location is inside MndRegion(c)
};

class Facility: public Point {
};

// Used in Replacement query, with a dnn field and a d2nn field
class ClientRepQuery: public Point {
public:
	float dnn;	// distance to its nearest facility
	float d2nn;	// distance to its second nearest facility

	ClientRepQuery();
	ClientRepQuery(float x, float y, float dnnIn, float d2nnIn);
};

typedef struct {
	Facility f;
	PotentialLocation p;
} opt_facility_location_pair_t;


/*
	Compute point-point distance
	Input: x1, y1 (point 1), x2, y2 (point 2)
	Output: sqrt((x1 - x2)^2 + (y1 - y2)^2)
*/
float pointDist(float x1, float y1, float x2, float y2);

/*
	Check whether a point is in an MBR
	Input: x, y (point), fMbr (MBR)
	Output: 1 if in, 0 if not
*/
int ptInMbr(float x, float y, float* fMbr);

/*
	Compute point-MBR distance
	Input: P (point), pMbr (MBR)
	Output: the distance between a point and an MBR
*/
float ptMbrDist(Point* p, float* pMbr);

/*
	Compute point-MBR distance
	Input: x1, y1 (point), pMbr (MBR)
	Output: the distance between a point and an MBR
*/
float ptMbrDist(float x, float y, float* pMbr);

/*
	Check MBR intersection
	Input: fMbr1 (MBR 1), fMbr2 (MBR 2)
	Output: 1 if intersect, 0 if not
*/
int mbrIntersect(float* fMbr1, float* fMbr2);

/*
	Get MBR for an array of points
	Input: pPts (point array), nLen (size of array)
	Output: fMbr (the MBR for pPts)
*/
void getMbr(Point* pPts, int nLen, float* fMbr);

#endif
