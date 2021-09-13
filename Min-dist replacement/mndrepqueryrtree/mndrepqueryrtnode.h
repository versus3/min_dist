#ifndef __CMNDRTNode_REPQUERY_H
#define __CMNDRTNode_REPQUERY_H
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class LinList;
class SortedLinList;
class CMNDRepqueryEntry;
class CMNDRepqueryRTree;
//------------------------------------------------------------
class CMNDRepqueryRTNode
{
public:
//--===on disk===--
	char level; 
	int block;
	int num_entries;
	CMNDRepqueryEntry *entries;
//--===others===--
	bool dirty;
	int capacity;
    int dimension;
	CMNDRepqueryRTree *my_tree;  
//--===functions===--
	CMNDRepqueryRTNode(CMNDRepqueryRTree *rt);
    CMNDRepqueryRTNode(CMNDRepqueryRTree *rt, int _block);
    ~CMNDRepqueryRTNode();

    int choose_subtree(float *brm);
	R_DELETE delete_entry(CMNDRepqueryEntry *e); 
	void enter(CMNDRepqueryEntry *de);
	bool FindLeaf(CMNDRepqueryEntry *e);
	float *get_mbr();
	int get_num_of_data();
	R_OVERFLOW insert(CMNDRepqueryEntry *d, CMNDRepqueryRTNode **sn);
	bool is_data_node() { return (level==0) ? TRUE : FALSE ;};
	void print();
	void rangeQuery(float *mbr, SortedLinList *res);
    void read_from_buffer(char *buffer);
	int split(float **mbr, int **distribution);
	void split(CMNDRepqueryRTNode *sn);
    void write_to_buffer(char *buffer); 
	
	//--added for valdity region queries---
	void rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far);
	void rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs);
};

#endif