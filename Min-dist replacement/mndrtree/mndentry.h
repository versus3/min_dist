/*entry.h
  this file defines class Entry*/
#ifndef __CEntry 
#define __CEntry
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class CMNDRTNode;
class CMNDRTree;
struct Linkable;
//------------------------------------------------------------
class CMNDEntry 
{
public:
//--===on disk===--
	int son;
	float *bounces;
	// qjz
	float fmnd;
	// qjz
//--===others===--
	int dimension;                      
	int level;
    CMNDRTree *my_tree;                     
    CMNDRTNode *son_ptr;                    
//for TP queries
	float dist; // the threshold of the entry when output
	int nn_num; // the nn that should be replaced 
   
//--===functions===--
	CMNDEntry();
	CMNDEntry(int dimension, CMNDRTree *rt);
    ~CMNDEntry();

	void del_son();
	Linkable *gen_Linkable();
	int get_size(); 
	CMNDRTNode *get_son();
	void init_entry(int _dimension, CMNDRTree *_rt);
	void read_from_buffer(char *buffer);// reads data from buffer
    SECTION section(float *mbr);        // tests, if mbr intersects the box
	bool section_circle(float *center, float radius);
	void set_from_Linkable(Linkable *link);
    void write_to_buffer(char *buffer); // writes data to buffer

    virtual CMNDEntry & operator = (CMNDEntry &_d);
	bool operator == (CMNDEntry &_d);
};

#endif