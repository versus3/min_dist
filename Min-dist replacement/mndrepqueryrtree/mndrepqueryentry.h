/*entry.h
  this file defines class Entry*/
#ifndef __CMNDEntry_H
#define __CMNDEntry_H
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class CMNDRepqueryRTNode;
class CMNDRepqueryRTree;
struct Linkable;
//------------------------------------------------------------
class CMNDRepqueryEntry 
{
public:
//--===on disk===--
	int son;
	float *bounces;
	// qjz
	float fmnd;
	float f2mnd;
	// qjz
//--===others===--
	int dimension;                      
	int level;
    CMNDRepqueryRTree *my_tree;                     
    CMNDRepqueryRTNode *son_ptr;                    
//for TP queries
	float dist; // the threshold of the entry when output
	int nn_num; // the nn that should be replaced 
   
//--===functions===--
	CMNDRepqueryEntry();
	CMNDRepqueryEntry(int dimension, CMNDRepqueryRTree *rt);
    ~CMNDRepqueryEntry();

	void del_son();
	Linkable *gen_Linkable();
	int get_size(); 
	CMNDRepqueryRTNode *get_son();
	void init_entry(int _dimension, CMNDRepqueryRTree *_rt);
	void read_from_buffer(char *buffer);// reads data from buffer
    SECTION section(float *mbr);        // tests, if mbr intersects the box
	bool section_circle(float *center, float radius);
	void set_from_Linkable(Linkable *link);
    void write_to_buffer(char *buffer); // writes data to buffer

    virtual CMNDRepqueryEntry & operator = (CMNDRepqueryEntry &_d);
	bool operator == (CMNDRepqueryEntry &_d);
};

#endif