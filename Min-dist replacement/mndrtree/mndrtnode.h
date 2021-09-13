#ifndef __CRTNode
#define __CRTNode
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class LinList;
class SortedLinList;
class CMNDEntry;
class CMNDRTree;
//------------------------------------------------------------
class CMNDRTNode
{
public:
//--===on disk===--
	char level; 
	int block;
	int num_entries;
	CMNDEntry *entries;
//--===others===--
	bool dirty;
	int capacity;
    int dimension;
	CMNDRTree *my_tree;  
//--===functions===--
	CMNDRTNode(CMNDRTree *rt);
    CMNDRTNode(CMNDRTree *rt, int _block);
    ~CMNDRTNode();

    int choose_subtree(float *brm);
	R_DELETE delete_entry(CMNDEntry *e); 
	void enter(CMNDEntry *de);
	bool FindLeaf(CMNDEntry *e);
	float *get_mbr();
	int get_num_of_data();
	R_OVERFLOW insert(CMNDEntry *d, CMNDRTNode **sn);
	bool is_data_node() { return (level==0) ? TRUE : FALSE ;};
	void print();
	void rangeQuery(float *mbr, SortedLinList *res);
    void read_from_buffer(char *buffer);
	int split(float **mbr, int **distribution);
	void split(CMNDRTNode *sn);
    void write_to_buffer(char *buffer); 
	
	//--added for valdity region queries---
	void rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far);
	void rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs);
};

#endif