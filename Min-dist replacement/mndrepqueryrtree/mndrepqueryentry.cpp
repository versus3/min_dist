/* entry.cpp
   implementation of class Entry */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mndrepqueryentry.h"
#include "mndrepqueryrtnode.h"
#include "../linlist/linlist.h"
//------------------------------------------------------------
CMNDRepqueryEntry::CMNDRepqueryEntry()
  //this ructor does nothing.  remember you should call init_entry
  //to initialize if you use this ructor
{
	son_ptr = NULL;
	bounces = NULL;
	nn_num=-1;

	// qjz
	fmnd = 0;
	f2mnd = 0;
	// qjz
}
//------------------------------------------------------------
CMNDRepqueryEntry::CMNDRepqueryEntry(int _dimension, CMNDRepqueryRTree *rt)
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
	f2mnd = 0;
	// qjz
}
//------------------------------------------------------------
CMNDRepqueryEntry::~CMNDRepqueryEntry()
{
    if (bounces)
		delete [] bounces;
    if (son_ptr != NULL)
	    delete son_ptr;
}
//------------------------------------------------------------
void CMNDRepqueryEntry::del_son()
{
	if (son_ptr != NULL)
	{
		delete son_ptr;
		son_ptr = NULL;
	}
}
//------------------------------------------------------------
Linkable* CMNDRepqueryEntry::gen_Linkable()
{
	Linkable *new_link = new Linkable(dimension);
	new_link -> son = son;
	//memcpy(new_link -> bounces, bounces, 2 * dimension * sizeof(float));
	for (int i = 0; i < 2 * dimension; i ++)
		new_link -> bounces[i] = bounces[i];
	new_link -> level = level;

	// qjz
	new_link->dd = fmnd;
	new_link->dd2 = f2mnd;
	// qjz
	return new_link;
}
//------------------------------------------------------------
int CMNDRepqueryEntry::get_size()
{
    //return 2 * dimension * sizeof(float) + sizeof(int);
	  //for bounces and son

	// qjz
	return 2 * dimension * sizeof(float) + sizeof(int) + sizeof(float) + sizeof(float);
	//for bounces and son and mnd and 2mnd
	// qjz
}
//------------------------------------------------------------
CMNDRepqueryRTNode* CMNDRepqueryEntry::get_son()
{
    if (son_ptr == NULL)
	    son_ptr = new CMNDRepqueryRTNode(my_tree, son);

    return son_ptr;
}
//------------------------------------------------------------
void CMNDRepqueryEntry::init_entry(int _dimension, CMNDRepqueryRTree *_rt)
{
	dimension = _dimension;
    my_tree = _rt;
    bounces = new float[2 * dimension];
    son_ptr = NULL;
    son = 0;
	level = 0;

	// qjz
	fmnd = 0;
	f2mnd = 0;
	// qjz
}
//------------------------------------------------------------
void CMNDRepqueryEntry::read_from_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(bounces, buffer, i);

	// qjz
	memcpy(&fmnd, &buffer[i], sizeof(float));
	i += sizeof(float);
	memcpy(&f2mnd, &buffer[i], sizeof(int));
	i += sizeof(float);
	// qjz
	memcpy(&son, &buffer[i], sizeof(int));
    i += sizeof(int);
}
//------------------------------------------------------------
SECTION CMNDRepqueryEntry::section(float *mbr)
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
bool CMNDRepqueryEntry::section_circle(float *center, float radius)
{
	float r2;

	r2 = radius * radius;

	if ((r2 - MINDIST(center,bounces,dimension)) < FLOATZERO)
		return TRUE;
	else
		return FALSE;
}
//------------------------------------------------------------
void CMNDRepqueryEntry::set_from_Linkable(Linkable *link)
{
	son = link -> son;
	dimension = link -> dimension;
	memcpy(bounces, link -> bounces, 2 * dimension * sizeof(float));
	level = link -> level;

	// qjz
	fmnd = link->dd;
	f2mnd = link->dd2;
	// qjz
	my_tree = NULL;
	son_ptr = NULL;
}
//------------------------------------------------------------
void CMNDRepqueryEntry::write_to_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(buffer, bounces, i);

	// qjz
	memcpy(&buffer[i], &fmnd, sizeof(float));
	i += sizeof(float);
	memcpy(&buffer[i], &f2mnd, sizeof(int));
	i += sizeof(float);
	// qjz
    memcpy(&buffer[i], &son, sizeof(int));
    i += sizeof(int);
}
//------------------------------------------------------------
bool CMNDRepqueryEntry::operator == (CMNDRepqueryEntry &_d)
  //this function compares two entries based on (1)son (2)dimension (3)extents
{
	if (son != _d.son) return false;
	if (dimension != _d.dimension) return false;
	// qjz
	if (fmnd != _d.fmnd) return false;
	if (f2mnd != _d.f2mnd) return false;
	// qjz

	for (int i = 0; i < 2 * dimension; i++)
		if (fabs(bounces[i] - _d.bounces[i]) > FLOATZERO) return false;
	return true;
}
//------------------------------------------------------------
CMNDRepqueryEntry& CMNDRepqueryEntry::operator = (CMNDRepqueryEntry &_d)
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
	f2mnd = _d.f2mnd;
	// qjz
    return *this;
}
//------------------------------------------------------------