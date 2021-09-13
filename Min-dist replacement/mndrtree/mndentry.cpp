/* entry.cpp
   implementation of class Entry */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mndentry.h"
#include "mndrtnode.h"
#include "../linlist/linlist.h"
//------------------------------------------------------------
CMNDEntry::CMNDEntry()
  //this ructor does nothing.  remember you should call init_entry
  //to initialize if you use this ructor
{
	son_ptr = NULL;
	bounces = NULL;
	nn_num=-1;

	// qjz
	fmnd = 0;
	// qjz
}
//------------------------------------------------------------
CMNDEntry::CMNDEntry(int _dimension, CMNDRTree *rt)
{
    dimension = _dimension;
    my_tree = rt;
    bounces = new float[2*dimension];
    son_ptr = NULL;
    son = 0;
	level = 0;
	nn_num=-1;

	// qjz
	fmnd = 0;
	// qjz
}
//------------------------------------------------------------
CMNDEntry::~CMNDEntry()
{
    if (bounces)
		delete [] bounces;
    if (son_ptr != NULL)
	    delete son_ptr;
}
//------------------------------------------------------------
void CMNDEntry::del_son()
{
	if (son_ptr != NULL)
	{
		delete son_ptr;
		son_ptr = NULL;
	}
}
//------------------------------------------------------------
Linkable* CMNDEntry::gen_Linkable()
{
	Linkable *new_link = new Linkable(dimension);
	new_link -> son = son;
	//memcpy(new_link -> bounces, bounces, 2 * dimension * sizeof(float));
	for (int i = 0; i < 2 * dimension; i ++)
		new_link -> bounces[i] = bounces[i];
	new_link -> level = level;

	// qjz
	new_link->dd = fmnd;
	// qjz
	return new_link;
}
//------------------------------------------------------------
int CMNDEntry::get_size()
{
    //return 2 * dimension * sizeof(float) + sizeof(int);
	  //for bounces and son

	// qjz
	return 2 * dimension * sizeof(float) + sizeof(int) + sizeof(float);
	//for bounces and son and dd
	// qjz
}
//------------------------------------------------------------
CMNDRTNode* CMNDEntry::get_son()
{
    if (son_ptr == NULL)
	    son_ptr = new CMNDRTNode(my_tree, son);

    return son_ptr;
}
//------------------------------------------------------------
void CMNDEntry::init_entry(int _dimension, CMNDRTree *_rt)
{
	dimension = _dimension;
    my_tree = _rt;
    bounces = new float[2 * dimension];
    son_ptr = NULL;
    son = 0;
	level = 0;

	// qjz
	fmnd = 0;
	// qjz
}
//------------------------------------------------------------
void CMNDEntry::read_from_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(bounces, buffer, i);

	// qjz
	memcpy(&fmnd, &buffer[i], sizeof(float));
	i += sizeof(float);
	// qjz
	memcpy(&son, &buffer[i], sizeof(int));
    i += sizeof(int);
}
//------------------------------------------------------------
SECTION CMNDEntry::section(float *mbr)
{
    bool inside;
    bool overlap;

    overlap = TRUE;
    inside = TRUE;

    for (int i = 0; i < dimension; i++)
    {
		if (mbr[2 * i] > bounces[2 * i + 1] ||  mbr[2 * i + 1] < bounces[2 * i])
			overlap = FALSE;
		if (mbr[2 * i] < bounces[2 * i] ||
			mbr[2 * i + 1] > bounces[2 * i + 1])
			inside = FALSE;
    }
    if (inside)
		return INSIDE;
    else if (overlap)
		return OVERLAP;
    else
		return S_NONE;
}
//------------------------------------------------------------
bool CMNDEntry::section_circle(float *center, float radius)
{
	float r2;

	r2 = radius * radius;

	if ((r2 - MINDIST(center,bounces,dimension)) < FLOATZERO)
		return TRUE;
	else
		return FALSE;
}
//------------------------------------------------------------
void CMNDEntry::set_from_Linkable(Linkable *link)
{
	son = link -> son;
	dimension = link -> dimension;
	memcpy(bounces, link -> bounces, 2 * dimension * sizeof(float));
	level = link -> level;

	// qjz
	fmnd = link->dd;
	// qjz
	my_tree = NULL;
	son_ptr = NULL;
}
//------------------------------------------------------------
void CMNDEntry::write_to_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(buffer, bounces, i);

	// qjz
	memcpy(&buffer[i], &fmnd, sizeof(float));
	i += sizeof(float);
	// qjz
    memcpy(&buffer[i], &son, sizeof(int));
    i += sizeof(int);
}
//------------------------------------------------------------
bool CMNDEntry::operator == (CMNDEntry &_d)
  //this function compares two entries based on (1)son (2)dimension (3)extents
{
	if (son != _d.son) return false;
	if (dimension != _d.dimension) return false;
	// qjz
	if (fmnd != _d.fmnd) return false;
	// qjz

	for (int i = 0; i < 2 * dimension; i++)
		if (fabs(bounces[i] - _d.bounces[i]) > FLOATZERO) return false;
	return true;
}
//------------------------------------------------------------
CMNDEntry& CMNDEntry::operator = (CMNDEntry &_d)
  //this function assigns all fieds of _d with the same values of this entry
{
    dimension = _d.dimension;
    son = _d.son;
    son_ptr = _d.son_ptr;
    memcpy(bounces, _d.bounces, sizeof(float) * 2 * dimension);
    my_tree = _d.my_tree;
	level = _d.level;

	nn_num=_d.nn_num;

	// qjz
	fmnd = _d.fmnd;
	// qjz
    return *this;
}
//------------------------------------------------------------