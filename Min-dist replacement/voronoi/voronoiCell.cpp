// voronoiCell.h
// Voronoi Cell Class
//
// Jianzhong Qi
// Modification Histry:
//		Date Created : 12/06/2010
//		Last Modified: 12/06/2010

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "voronoiCell.h"

// class Line
Line::Line(float x1, float y1, float x2, float y2)
{
	m_pt1.x = x1;
	m_pt1.y = y1;

	m_pt2.x = x2;
	m_pt2.y = y2;

	if(m_pt1.x == m_pt2.x){
		m_isUpRight = 1;
		
		// has two types of forms: X = x1, a1X + a2Y + a3 = 0
		m_fa1 = 1;
		m_fa2 = 0;
		m_fa3 = -m_pt1.x;
	}else{
		m_isUpRight = 0;
		
		// has three types of forms: Y = kX + b, a1X + a2Y + a3 = 0, Y - y1 = k(X- x1)
		m_fk = (m_pt2.y - m_pt1.y)/(m_pt2.x - m_pt1.x);
		m_fb = m_pt2.y - m_fk * m_pt2.x;

		m_fa1 = -m_fk;
		m_fa2 = 1;
		m_fa3 = -m_fb;

	}

}

Line::Line(Point* pt1, Point* pt2)
{
	m_pt1 = *pt1;
	
	m_pt2 = *pt2;

	if(m_pt1.x == m_pt2.x){
		m_isUpRight = 1;

		// has two types of forms: X = x1, a1X + a2Y + a3 = 0
		m_fa1 = 1;
		m_fa2 = 0;
		m_fa3 = -m_pt1.x;
	}else{
		m_isUpRight = 0;

		// has three types of forms: Y = kX + b, a1X + a2Y + a3 = 0, Y - y1 = k(X- x1)
		m_fk = (m_pt2.y - m_pt1.y)/(m_pt2.x - m_pt1.x);
		m_fb = m_pt2.y - m_fk * m_pt2.x;

		m_fa1 = -m_fk;
		m_fa2 = 1;
		m_fa3 = -m_fb;

	}
}

Line::Line(float x1, float y1, float k)
{
	m_pt1.x = x1;
	m_pt1.y = y1;

	m_isUpRight = 0;

	m_fk = k;
	m_fb = m_pt1.y - m_fk * m_pt1.x;

	m_fa1 = -m_fk;
	m_fa2 = 1;
	m_fa3 = -m_fb;

	m_pt2.x = m_pt1.x + 1;
	m_pt2.y = m_fk * m_pt2.x + m_fb;
}

Line::~Line()
{
	return;
}

int Line::Intersect(Line* pLine, Point* ptIntersect)
{
	if(this->m_isUpRight){
		if(pLine->m_isUpRight){
			if(this->m_pt1.x == pLine->m_pt1.x){
				ptIntersect->x = this->m_pt1.x;
				ptIntersect->y = this->m_pt1.y;
				return 2;
			}else{
				return 0;
			}
		}else{
			ptIntersect->x = this->m_pt1.x;
			ptIntersect->y = pLine->m_fk * this->m_pt1.x + pLine->m_fb;
			return 1;
		}
	}else{
		if(pLine->m_isUpRight){
			ptIntersect->x = pLine->m_pt1.x;
			ptIntersect->y = this->m_fk * pLine->m_pt1.x + this->m_fb;
			return 1;
		}else{
			if(this->m_fk == pLine->m_fk){
				if(this->m_fb == pLine->m_fb){
					ptIntersect->x = this->m_pt1.x;
					ptIntersect->y = this->m_pt1.y;
					return 2;
				}else{
					return 0;
				}
			}else{
				ptIntersect->y = (pLine->m_fa3*this->m_fa1 - this->m_fa3*pLine->m_fa1) / (pLine->m_fa1*this->m_fa2 - this->m_fa1*pLine->m_fa2);
				if(this->m_fa1 == 0){
					ptIntersect->x = (ptIntersect->y * pLine->m_fa2 + pLine->m_fa3) / (-pLine->m_fa1);
				}else{
					ptIntersect->x = (ptIntersect->y * this->m_fa2 + this->m_fa3) / (-this->m_fa1);
				}
				return 1;
			}
		}
	}
	return 0;
}

/*
	Point m_pt1, m_pt2;
	float m_fk, m_fb;			// Y = kX + b
	int m_isUpRight;	
	// float m_fx1 = m_pt1.x;	// X = x1
	float m_fa1, m_fa2, m_fa3;	// a1X + a2Y + a3 = 0
	// float m_fx1 = m_pt1.x, m_fy1 = m_pt1.y;		// Y - y1 = k(X- x1)
*/

// class VoronoiCell

// create a VoronoiCell, the first point in pt[] is the center point of the cell
VoronoiCell::VoronoiCell(Point* pt[], int n, float* fMbr)
{
	int i;
	Line* pLine;
	float fxMin, fxMax, fyMin, fyMax;
	float fMidX, fMidY, fk;
	Point interPt;
	int isFirst = 1;

	m_fDd = 0;

	i = 0;
	m_points.push_back(*pt[i]);

	i++;
	while(i < n){
		m_points.push_back(*pt[i]);

		// create line
		if(pt[i]->x == pt[0]->x){
			// is Vertical
			pLine = new Line(pt[i]->x, float ((pt[i]->y + pt[0]->y) * 0.5), pt[i]->x+1, float ((pt[i]->y + pt[0]->y) * 0.5));
		}else if(pt[i]->y == pt[0]->y){
			// is Horizontal
			pLine = new Line(float ((pt[i]->x + pt[0]->x) * 0.5), pt[i]->y, float ((pt[i]->x + pt[0]->x) * 0.5), pt[i]->y+1);
		}else{
			fMidX = float ((pt[i]->x + pt[0]->x) * 0.5);
			fMidY = float ((pt[i]->y + pt[0]->y) * 0.5);
			fk = -(pt[0]->x - pt[i]->x) / (pt[0]->y - pt[i]->y);

			pLine = new Line(fMidX, fMidY, fk);
		}
		m_lines.push_back(pLine);
		
		if(i > 1){
			// get intersection point
			if(m_lines[i-1]->Intersect(m_lines[i-2], &interPt)){
				this->m_intersectPts.push_back(interPt);
				if(isFirst){
					fxMin = fxMax = interPt.x;
					fyMin = fyMax = interPt.y;
					isFirst = 0;
				}else{
					if(interPt.x < fxMin){
						fxMin = interPt.x;
					}
					if(interPt.x > fxMax){
						fxMax = interPt.x;
					}
					if(interPt.y < fyMin){
						fyMin = interPt.y;
					}
					if(interPt.y > fyMax){
						fyMax = interPt.y;
					}
				}
			}
		}
		i++;
	}
	
	if(m_lines.size() > 2){
		// get intersection point
		if(m_lines[m_lines.size()-1]->Intersect(m_lines[0], &interPt)){
			this->m_intersectPts.push_back(interPt);
			if(isFirst){
				fxMin = fxMax = interPt.x;
				fyMin = fyMax = interPt.y;
				isFirst = 0;
			}else{
				if(interPt.x < fxMin){
					fxMin = interPt.x;
				}
				if(interPt.x > fxMax){
					fxMax = interPt.x;
				}
				if(interPt.y < fyMin){
					fyMin = interPt.y;
				}
				if(interPt.y > fyMax){
					fyMax = interPt.y;
				}
			}
		}

	}
	
	// get mbr
	if(!isFirst){
		// CreateMBR with no intersection point
		CreateMBR(fxMin, fxMax, fyMin, fyMax, fMbr);
	}else{
		// CreateMBR with at least one intersection point
		CreateMBR(fMbr);
	}
}

VoronoiCell::~VoronoiCell()
{
	unsigned int i;
	m_points.clear();
	m_intersectPts.clear();

	i = 0;
	while(i < m_lines.size()){
		delete m_lines[i];
		i++;
	}
	m_lines.clear();
	//float m_fMbr[4];
}

void VoronoiCell::CreateMBR(float fxMin, float fxMax, float fyMin, float fyMax, float* fMbr)
{
	int i;
	Line* pLine[3];
	float fMbr1[4];
					
	switch(m_lines.size()){
		case 2:
			// two intersecting lines
			pLine[0] = this->m_lines[0];
			pLine[1] = this->m_lines[1];

			m_lines.pop_back();
			this->CreateMBR(fMbr);
			for(i = 0; i < 4; i++){
				fMbr1[i] = this->m_fMbr[i];
			}

			m_lines[0] = pLine[1];
			this->CreateMBR(fMbr);
			m_lines[0] = pLine[0];
			m_lines.push_back(pLine[1]);

			if(fMbr1[0] == fMbr1[1]){
				return;
			}else if(this->m_fMbr[0] == this->m_fMbr[1]){
				for(i = 0; i < 4; i++){
					this->m_fMbr[i] = fMbr1[i];
				}
			}else{
				if(fMbr1[0] < this->m_fMbr[0]){
					this->m_fMbr[0] = fMbr1[0];
				}
				if(fMbr1[1] > this->m_fMbr[1]){
					this->m_fMbr[1] = fMbr1[1];
				}
				if(fMbr1[2] < this->m_fMbr[2]){
					this->m_fMbr[2] = fMbr1[2];
				}
				if(fMbr1[3] > this->m_fMbr[3]){
					this->m_fMbr[3] = fMbr1[3];
				}
			}
			
			break;
		case 3:
			// three intersecting lines
			pLine[0] = this->m_lines[0];
			pLine[1] = this->m_lines[1];
			pLine[2] = this->m_lines[2];

			m_lines.pop_back();
			this->CreateMBR(fxMin, fxMax, fyMin, fyMax, fMbr);
			for(i = 0; i < 4; i++){
				fMbr1[i] = this->m_fMbr[i];
			}

			m_lines[0] = pLine[1];
			m_lines[1] = pLine[2];
			this->CreateMBR(fxMin, fxMax, fyMin, fyMax, fMbr);

			m_lines[0] = pLine[0];
			m_lines[1] = pLine[1];
			m_lines.push_back(pLine[2]);

			if(fMbr1[0] == fMbr1[1]){
				return;
			}else if(this->m_fMbr[0] == this->m_fMbr[1]){
				for(i = 0; i < 4; i++){
					this->m_fMbr[i] = fMbr1[i];
				}
			}else{
				if(fMbr1[0] < this->m_fMbr[0]){
					this->m_fMbr[0] = fMbr1[0];
				}
				if(fMbr1[1] > this->m_fMbr[1]){
					this->m_fMbr[1] = fMbr1[1];
				}
				if(fMbr1[2] < this->m_fMbr[2]){
					this->m_fMbr[2] = fMbr1[2];
				}
				if(fMbr1[3] > this->m_fMbr[3]){
					this->m_fMbr[3] = fMbr1[3];
				}
			}

			break;
		case 4:
			m_fMbr[0] = fxMin;
			m_fMbr[1] = fxMax;
			m_fMbr[2] = fyMin;
			m_fMbr[3] = fyMax;
			return;
	}
	
}

void VoronoiCell::CreateMBR(float* fMbr)
{
	int i;
	Point ptInt[4];
	float fSide;
	float fFlag[4];
	int nIntCount;
	Line* pLine[2];
	float fMbr1[4];

	Line line1(fMbr[0], fMbr[2], fMbr[0], fMbr[3]);
	Line line2(fMbr[0], fMbr[3], fMbr[1], fMbr[3]);
	Line line3(fMbr[1], fMbr[2], fMbr[1], fMbr[3]);
	Line line4(fMbr[0], fMbr[2], fMbr[1], fMbr[2]);

	m_fMbr[0] = m_fMbr[1] = m_points[0].x;
	m_fMbr[2] = m_fMbr[3] = m_points[0].y;
	
	switch(m_lines.size()){
		case 1:
			// only one line
			if(m_lines[0]->m_isUpRight){
				// a vertical line
				if(this->m_points[0].x - m_lines[0]->m_pt1.x > 0){
					if(fMbr[0] >= m_lines[0]->m_pt1.x){
						for(i = 0; i < 4; i++){
							m_fMbr[i] = fMbr[i];
						}
						return; 
					}else if(fMbr[1] > m_lines[0]->m_pt1.x){
						for(i = 0; i < 4; i++){
							m_fMbr[i] = fMbr[i];
						}
						m_fMbr[0] = m_lines[0]->m_pt1.x;
						return; 
					}
				}else{
					if(fMbr[1] <= m_lines[0]->m_pt1.x){
						for(i = 0; i < 4; i++){
							m_fMbr[i] = fMbr[i];
						}
						return; 
					}else if(fMbr[0] < m_lines[0]->m_pt1.x){
						for(i = 1; i < 4; i++){
							m_fMbr[i] = fMbr[i];
						}
						m_fMbr[1] = m_lines[0]->m_pt1.x;
						return; 
					}				
				}
			}else if(m_lines[0]->m_fk == 0){
				// a horizontal line
				if(this->m_points[0].y - m_lines[0]->m_fb > 0){
					if(fMbr[2] >= m_lines[0]->m_fb){
						for(i = 0; i < 4; i++){
							m_fMbr[i] = fMbr[i];
						}
						return; 
					}else if(fMbr[3] > m_lines[0]->m_fb){
						for(i = 0; i < 4; i++){
							m_fMbr[i] = fMbr[i];
						}
						m_fMbr[2] = m_lines[0]->m_fb;
						return; 
					}
				}else{
					if(fMbr[3] <= m_lines[0]->m_fb){
						for(i = 0; i < 4; i++){
							m_fMbr[i] = fMbr[i];
						}
						return; 
					}else if(fMbr[2] < m_lines[0]->m_fb){
						for(i = 1; i < 4; i++){
							m_fMbr[i] = fMbr[i];
						}
						m_fMbr[3] = m_lines[0]->m_fb;
						return; 
					}				
				}
			}else{
				fFlag[0] = m_lines[0]->m_fa1 * fMbr[0] + m_lines[0]->m_fa2 * fMbr[2] + m_lines[0]->m_fa3; 
				fFlag[1] = m_lines[0]->m_fa1 * fMbr[0] + m_lines[0]->m_fa2 * fMbr[3] + m_lines[0]->m_fa3;
				fFlag[2] = m_lines[0]->m_fa1 * fMbr[1] + m_lines[0]->m_fa2 * fMbr[3] + m_lines[0]->m_fa3;
				fFlag[3] = m_lines[0]->m_fa1 * fMbr[1] + m_lines[0]->m_fa2 * fMbr[2] + m_lines[0]->m_fa3;
				fSide = m_lines[0]->m_fa1 * m_points[0].x + m_lines[0]->m_fa2 * m_points[0].y + m_lines[0]->m_fa3;

				nIntCount = 0;
				m_lines[0]->Intersect(&line1, &ptInt[0]);
				m_lines[0]->Intersect(&line2, &ptInt[1]);
				m_lines[0]->Intersect(&line3, &ptInt[2]);
				m_lines[0]->Intersect(&line4, &ptInt[3]);
				if(ptInt[0].y > fMbr[2] && ptInt[0].y <= fMbr[3]){
					nIntCount++;
				}
				if(ptInt[1].x > fMbr[0] && ptInt[1].x <= fMbr[1]){
					nIntCount++;
				}
				if(ptInt[2].y >= fMbr[2] && ptInt[2].y < fMbr[3]){
					nIntCount++;
				}
				if(ptInt[3].x >= fMbr[0] && ptInt[3].x < fMbr[1]){
					nIntCount++;
				}
				if((nIntCount == 0 && fFlag[0] * fSide > 0) ||
					(nIntCount == 1 && (fFlag[0] * fSide > 0 || fFlag[1] * fSide > 0))){
					for(i = 0; i < 4; i++){
						m_fMbr[i] = fMbr[i];
					}
					return;
				}else if(nIntCount == 2){
					for(i = 0; i < 4; i++){
						m_fMbr[i] = fMbr[i];
					}
					
					if(fFlag[0] == 0){
						if(ptInt[1].x > fMbr[0] && ptInt[1].x <= fMbr[1]){
							if(fFlag[1] * fSide > 0){
								m_fMbr[1] = ptInt[1].x;
							}
						}else{
							if(fFlag[3] * fSide > 0){
								m_fMbr[3] = ptInt[2].y;
							}
						}
					}else{
						if(ptInt[1].x >= fMbr[0] && ptInt[1].x <= fMbr[1] && 
							ptInt[3].x >= fMbr[0] && ptInt[3].x <= fMbr[1]){
							if(fFlag[0] * fSide > 0){
								m_fMbr[1] = ptInt[1].x > ptInt[3].x ? ptInt[1].x : ptInt[3].x;
							}else{
								m_fMbr[0] = ptInt[1].x < ptInt[3].x ? ptInt[1].x : ptInt[3].x;
							}
						}else if(ptInt[0].y >= fMbr[2] && ptInt[0].y <= fMbr[3] && 
							ptInt[2].y >= fMbr[2] && ptInt[2].y <= fMbr[3]){
							if(fFlag[0] * fSide > 0){
								m_fMbr[3] = ptInt[0].y > ptInt[2].y ? ptInt[0].y : ptInt[2].y;
							}else{
								m_fMbr[2] = ptInt[0].y < ptInt[2].y ? ptInt[0].y : ptInt[2].y;
							}
						}else if(ptInt[0].y >= fMbr[2] && ptInt[0].y <= fMbr[3]){
							if(ptInt[1].x >= fMbr[0] && ptInt[1].x <= fMbr[1] && fFlag[0] * fSide < 0){
								m_fMbr[1] = ptInt[1].x;
								m_fMbr[2] = ptInt[0].y;
							}else if(ptInt[3].x >= fMbr[0] && ptInt[3].x <= fMbr[1] && fFlag[0] * fSide > 0){
								m_fMbr[1] = ptInt[1].x;
								m_fMbr[3] = ptInt[0].y;
							}
						}else if(ptInt[2].y >= fMbr[2] && ptInt[2].y <= fMbr[3]){
							if(ptInt[1].x >= fMbr[0] && ptInt[1].x <= fMbr[1] && fFlag[0] * fSide < 0){
								m_fMbr[0] = ptInt[1].x;
								m_fMbr[2] = ptInt[2].y;
							}else if(ptInt[3].x >= fMbr[0] && ptInt[3].x <= fMbr[1] && fFlag[0] * fSide < 0){
								m_fMbr[0] = ptInt[1].x;
								m_fMbr[3] = ptInt[2].y;
							}
						}
					}
				}
			}

			break;
		case 2:
			// two parallel lines
			if(m_lines[0]->m_isUpRight){
				// is vertical
				this->m_fMbr[0] = m_lines[0]->m_pt1.x < m_lines[1]->m_pt1.x ? m_lines[0]->m_pt1.x : m_lines[1]->m_pt1.x;
				this->m_fMbr[1] = m_lines[0]->m_pt1.x + m_lines[1]->m_pt1.x - this->m_fMbr[0];
				this->m_fMbr[2] = fMbr[2];
				this->m_fMbr[3] = fMbr[3];
			}else if(m_lines[0]->m_fk == 0){
				// is horizontal
				this->m_fMbr[0] = fMbr[0];
				this->m_fMbr[1] = fMbr[1];
				this->m_fMbr[2] = m_lines[0]->m_fb < m_lines[1]->m_fb ? m_lines[0]->m_fb : m_lines[1]->m_fb;
				this->m_fMbr[3] = m_lines[0]->m_fb + m_lines[1]->m_fb - this->m_fMbr[2];
			}else{
				pLine[0] = this->m_lines[0];
				pLine[1] = this->m_lines[1];
				m_lines.pop_back();
				this->CreateMBR(fMbr);
				for(i = 0; i < 4; i++){
					fMbr1[i] = this->m_fMbr[i];
				}
				
				m_lines[0] = pLine[1];
				this->CreateMBR(fMbr);
				m_lines[0] = pLine[0];
				m_lines.push_back(pLine[1]);
				
				if(fMbr1[0] == fMbr1[1]){
					return;
				}else if(this->m_fMbr[0] == this->m_fMbr[1]){
					for(i = 0; i < 4; i++){
						this->m_fMbr[i] = fMbr1[i];
					}
				}else{
					if(fMbr1[0] > this->m_fMbr[0]){
						this->m_fMbr[0] = fMbr1[0];
					}
					if(fMbr1[1] < this->m_fMbr[1]){
						this->m_fMbr[1] = fMbr1[1];
					}
					if(fMbr1[2] > this->m_fMbr[2]){
						this->m_fMbr[2] = fMbr1[2];
					}
					if(fMbr1[3] < this->m_fMbr[3]){
						this->m_fMbr[3] = fMbr1[3];
					}
				}
			}
			
			break;
	}
}