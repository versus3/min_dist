/* bbrepquery.h
   this file defines the class CBBRepqeuryRTree*/

#ifndef __CMNDRTree_REPQUERY_H
#define __CMNDRTree_REPQUERY_H
//------------------------------------------------------------
#include "../func/gendef.h"
#include "../heap/heap.h"

//added for TP KNN------------------------------------------------
#define SEQUENCE_SENSITIVE false
  //set it true if you want a new partition point every time the order of the NNs change
#define PAST_PAIR 1000

//------------------------------------------------------------
class LinList;
class SortedLinList;
class Cache;
class CMNDRepqueryRTNode;
class CMNDRepqueryEntry;
//------------------------------------------------------------
class CMNDRepqueryRTree : public Cacheable
{
public:
//--===on disk===--
	int dimension;                       
	int num_of_data;	                 
    int num_of_dnodes;	                 
    int num_of_inodes;	                 
	int root;                            
	bool root_is_data;                   
//--===others===--
	CMNDRepqueryRTNode *root_ptr;
    bool *re_level;  
    LinList *re_data_cands; 
	LinList *deletelist;

//--===added for TP KNN===--
	int last_pair[PAST_PAIR][2]; //records the pairs of points that produce the minimum inflence time
	int lastcnt; //next last pair to be replaced
	Heap *tpheap;

//--===functions===--
	CMNDRepqueryRTree();
	CMNDRepqueryRTree(char *fname, int _b_length, Cache* c, int _dimension);
    CMNDRepqueryRTree(char *fname, Cache* c);
    CMNDRepqueryRTree(char *inpname, char *fname, int _blength, Cache* c, int _dimension);
    ~CMNDRepqueryRTree();
	void del_root();
	bool delete_entry(CMNDRepqueryEntry *d);
	bool FindLeaf(CMNDRepqueryEntry *e);
    int get_num() { return num_of_data; }
	void insert(CMNDRepqueryEntry *d);
	void load_root();  
	void rangeQuery(float *mbr, SortedLinList *res);
	void read_header(char *buffer);      
	void write_header(char *buffer);
	int update_rslt(CMNDRepqueryEntry *_e, float _dist, CMNDRepqueryEntry *_rslt, 
					 float *_key, int _k);

	// This function was added to perform TP-kNN queries by Bobby
	void TPNN_TP(float *_qline, int _k, CMNDRepqueryEntry *_nn, CMNDRepqueryEntry *_rslt, float _max_trvl);
	
	//--added for valdity region queries---
	void rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far);
	void rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs);
	void BFNN(float *_qpt, int _k, CMNDRepqueryEntry *_rslt);
	void BFNNCont(float *qmbr, float *qmbr2, Heap *heap, HeapEntry *e, int k);
};

#endif // __CRTree1
