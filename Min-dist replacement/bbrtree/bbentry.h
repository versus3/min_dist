/*entry.h
  this file defines class Entry*/
#ifndef __CBBEntry_H
#define __CBBEntry_H
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class CBBRTNode;
class CBBRTree;
struct Linkable;
//------------------------------------------------------------
class CBBEntry 
{
public:
//--===on disk===--
	int son;
	float *bounces;
	// qjz
	float fMaxFDist;
	int nDataEntryNum;
	// qjz
//--===others===--
	int dimension;                      
	int level;
    CBBRTree *my_tree;                     
    CBBRTNode *son_ptr;                    
//for TP queries
	float dist; // the threshold of the entry when output
	int nn_num; // the nn that should be replaced 
   
//--===functions===--
	CBBEntry();
	CBBEntry(int dimension, CBBRTree *rt);
    ~CBBEntry();

	void del_son();
	Linkable *gen_Linkable();
	int get_size(); 
	CBBRTNode *get_son();
	void init_entry(int _dimension, CBBRTree *_rt);
	void read_from_buffer(char *buffer);// reads data from buffer
    SECTION section(float *mbr);        // tests, if mbr intersects the box
	bool section_circle(float *center, float radius);
	void set_from_Linkable(Linkable *link);
    void write_to_buffer(char *buffer); // writes data to buffer

    virtual CBBEntry & operator = (CBBEntry &_d);
	bool operator == (CBBEntry &_d);
};

#endif