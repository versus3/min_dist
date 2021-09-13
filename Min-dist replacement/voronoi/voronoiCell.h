// voronoiCell.h
// Voronoi Cell Class
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 12/06/2010
//		Last Modified: 12/06/2010

#ifndef __VORONOI_CELL
#define __VORONOI_CELL

#include <vector>
#include "../point/point.h"

using namespace std;

class Line{
public:
	Line(float x1, float y1, float x2, float y2);
	Line(Point* pt1, Point* pt2);
	Line(float x1, float y1, float k);

	~Line();

	int Intersect(Line* pLine, Point* ptIntersect);

	Point m_pt1, m_pt2;
	float m_fk, m_fb;								// Y = kX + b
	int m_isUpRight;	
	// float m_fx1 = m_pt1.x;						// X = x1
	float m_fa1, m_fa2, m_fa3;						// a1X + a2Y + a3 = 0
	// float m_fx1 = m_pt1.x, m_fy1 = m_pt1.y;		// Y - y1 = k(X- x1)
};

class VoronoiCell{
	
public:
	// create a VoronoiCell, the first point in pt[] is the center point of the cell
	VoronoiCell(Point* pt[], int n, float* fMbr );

	~VoronoiCell();

	void CreateMBR(float fxMin, float fxMax, float fyMin, float fyMax, float* fMbr);

	void CreateMBR(float* fMbr);
public:
	vector<Point> m_points;
	vector<Point> m_intersectPts;
	vector<Line*> m_lines;
	float m_fMbr[4];

	float m_fDd;
};

#endif