#ifndef __CBBRTNode_H
#define __CBBRTNode_H
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class LinList;
class SortedLinList;
class CBBEntry;
class CBBRTree;
//------------------------------------------------------------
class CBBRTNode
{
public:
//--===on disk===--
	char level; 
	int block;
	int num_entries;
	CBBEntry *entries;
//--===others===--
	bool dirty;
	int capacity;
    int dimension;
	CBBRTree *my_tree;  
//--===functions===--
	CBBRTNode(CBBRTree *rt);
    CBBRTNode(CBBRTree *rt, int _block);
    ~CBBRTNode();

    int choose_subtree(float *brm);
	R_DELETE delete_entry(CBBEntry *e); 
	void enter(CBBEntry *de);
	bool FindLeaf(CBBEntry *e);
	float *get_mbr();
	int get_num_of_data();
	R_OVERFLOW insert(CBBEntry *d, CBBRTNode **sn);
	bool is_data_node() { return (level==0) ? TRUE : FALSE ;};
	void print();
	void rangeQuery(float *mbr, SortedLinList *res);
    void read_from_buffer(char *buffer);
	int split(float **mbr, int **distribution);
	void split(CBBRTNode *sn);
    void write_to_buffer(char *buffer); 
	
	//--added for valdity region queries---
	void rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far);
	void rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs);
};

#endif