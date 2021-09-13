#pragma warning(disable:4996)
/*CMNDRTree.cpp
  this file implements the CMNDRTree class*/
#include <math.h>
#include <string.h>
#include "mndrtree.h"
#include "mndentry.h"
#include "mndrtnode.h"
#include "../rtree/distance.h"
#include "../blockfile/cache.h"
#include "../blockfile/blk_file.h"
#include "../linlist/linlist.h"

//------------------------------------------------------------
CMNDRTree::CMNDRTree()
{
}
//------------------------------------------------------------
CMNDRTree::CMNDRTree(char *fname, int _b_length, Cache *c, int _dimension)
  //use this constructor to build a new tree
{
    file = new BlockFile(fname, _b_length);
    cache = c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new CMNDRTNode(this);
	  //note that when a tree is constructed, the root is automatically created
	  //though at this time there is no CMNDEntry in the root yet.
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr -> block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;
}
//------------------------------------------------------------
CMNDRTree::CMNDRTree(char *fname, Cache *c)
  //use this constructor to restore a tree from a file
{
    file = new BlockFile(fname, 0);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    file -> read_header(header);
    read_header(header);
	delete [] header;

    root_ptr = NULL;

	//added for TP KNN----------------------------------------
	tpheap=NULL;
}
//------------------------------------------------------------
CMNDRTree::CMNDRTree(char *inpname, char *fname, int _b_length, Cache *c, int _dimension)
  // construct new R-tree from a specified input textfile with rectangles
{
    CMNDEntry *d;
    FILE *fp;

    file = new BlockFile(fname, _b_length);
    cache =c;

    re_data_cands = new LinList();
	deletelist = new LinList();

    char *header = new char [file->get_blocklength()];
    read_header(header);
	delete [] header;

    dimension = _dimension;
    root = 0;
    root_ptr = NULL;
    root_is_data = TRUE;
    num_of_data = num_of_inodes = num_of_dnodes = 0;

    root_ptr = new CMNDRTNode(this);
    num_of_dnodes++;
    root_ptr -> level = 0;
    root = root_ptr->block;

	//added for TP KNN----------------------------------------
	tpheap=NULL;

	int record_count = 0;

    if((fp = fopen(inpname,"r")) == NULL)
    {
      delete this;
      error("Cannot open R-Tree text file", TRUE);
    }
    else
    {
      while (!feof(fp))
      {
		record_count ++;

		d = new CMNDEntry(dimension, NULL);

    	fscanf(fp, "%d", &(d -> son));
		for (int i = 0; i < 2 * dimension; i ++)
			fscanf(fp, " %f", &(d -> bounces[i]));
		//fscanf(fp, " %f %f %f %f\n", &(d->bounces[0]),
		//	&(d->bounces[2]), &(d->bounces[1]), &(d->bounces[3]));
		fscanf(fp, "\n");

    	insert(d);
		  //d will be deleted in insert()

		if (record_count % 100 == 0)
		{
			for (int i = 0; i < 79; i ++)  //clear a line
				printf("\b");

			printf("inserting object %d", record_count);
		}
      }
    }

	fclose(fp);

	printf("\n");
	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
CMNDRTree::~CMNDRTree()
{
	char *header = new char[file -> get_blocklength()];
    write_header(header);
    file->set_header(header);
    delete [] header;

    if (root_ptr != NULL)
    {
        delete root_ptr;
        root_ptr = NULL;
    }

	if (cache)
      cache -> flush();

    delete file;

    delete re_data_cands;
	delete deletelist;

    //printf("This R-Tree contains %d internal, %d data nodes and %d data\n",
	//   num_of_inodes, num_of_dnodes, num_of_data);
}
//------------------------------------------------------------
void CMNDRTree::del_root()
{
	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
bool CMNDRTree::delete_entry(CMNDEntry *d)
{
	load_root();

	R_DELETE del_ret;
	del_ret=root_ptr->delete_entry(d);

	if (del_ret == NOTFOUND) return false;
	if (del_ret == ERASED) 
		error("CMNDRTree::delete_entry--The root has been deleted\n",true);
 
	if (root_ptr -> level > 0 && root_ptr -> num_entries == 1)
		//there is only one CMNDEntry in the root but the root
		//is not leaf.  in this case, the child of the root is exhalted to root
	{
		root = root_ptr -> entries[0].son;
		delete root_ptr;
		root_ptr = NULL;
		load_root();
		num_of_inodes--;
	}

	//Now will reinsert the entries
	while (deletelist -> get_num() > 0)
	{
		Linkable *e;
		e = deletelist -> get_first();
		CMNDEntry *new_e = new CMNDEntry(dimension, NULL);
		new_e -> set_from_Linkable(e);
		deletelist -> erase();
		insert(new_e);
	}

	delete root_ptr;
	root_ptr = NULL;

	return true;
}
//------------------------------------------------------------
void CMNDRTree::insert(CMNDEntry* d)
{
    int i, j;
    CMNDRTNode *sn;
    CMNDRTNode *nroot_ptr;
    int nroot;
    CMNDEntry *de;
    R_OVERFLOW split_root;
    CMNDEntry *dc;
    float *nmbr;

    // load root into memory
    load_root();

    // no overflow occured until now
    re_level = new bool[root_ptr -> level + 1];
    for (i = 0; i <= root_ptr -> level; i++)
        re_level[i] = FALSE;

    // insert d into re_data_cands as the first CMNDEntry to insert
    // make a copy of d because it should be erased later
    Linkable *new_link;
	new_link = d -> gen_Linkable();
	re_data_cands -> insert(new_link);

	delete d;  //we follow the convention that the CMNDEntry will be deleted when insertion finishes

    j = -1;
    while (re_data_cands -> get_num() > 0)
    {
        // first try to insert data, then directory entries
	    Linkable *d_cand;
		d_cand = re_data_cands -> get_first();
        if (d_cand != NULL)
        {
            // since "erase" deletes the data itself from the
            // list, we should make a copy of the data before
            // erasing it
			dc = new CMNDEntry(dimension, NULL);
            dc -> set_from_Linkable(d_cand);
            re_data_cands -> erase();

            // start recursive insert with root
			split_root = root_ptr -> insert(dc, &sn);
        }
        else
	        error("CMNDRTree::insert: inconsistent list re_data_cands", TRUE);

    	if (split_root == SPLIT)
    	// insert has lead to split --> new root-page with two sons (i.e. root and sn)
    	{
    	    nroot_ptr = new CMNDRTNode(this);
    	    nroot_ptr -> level = root_ptr -> level + 1;
    	    num_of_inodes++;
    	    nroot = nroot_ptr -> block;

    	    de = new CMNDEntry(dimension, this);
    	    nmbr = root_ptr -> get_mbr();
    	    memcpy(de->bounces, nmbr, 2*dimension*sizeof(float));
    	    delete [] nmbr;
    	    de->son = root_ptr->block;
    	    de->son_ptr = root_ptr;
    	    nroot_ptr -> enter(de);

    	    de = new CMNDEntry(dimension, this);
    	    nmbr = sn -> get_mbr();
    	    memcpy(de -> bounces, nmbr, 2*dimension*sizeof(float));
    	    delete [] nmbr;
    	    de -> son = sn -> block;
    	    de -> son_ptr = sn;
    	    nroot_ptr->enter(de);

    	    root = nroot;
            root_ptr = nroot_ptr;

            root_is_data = FALSE;
        }
        j++;
    }

    num_of_data++;

    delete [] re_level;

	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
void CMNDRTree::load_root()
{
    if (root_ptr == NULL)
        root_ptr = new CMNDRTNode(this, root);
}
//------------------------------------------------------------
void CMNDRTree::rangeQuery(float *mbr, SortedLinList *res)
{
    load_root();

    root_ptr -> rangeQuery(mbr,res);

	delete root_ptr;
	root_ptr = NULL;
}
//------------------------------------------------------------
void CMNDRTree::read_header(char *buffer)
{
    int i;

    memcpy(&dimension, buffer, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&num_of_data, &buffer[i], sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&num_of_dnodes, &buffer[i], sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&num_of_inodes, &buffer[i], sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&root_is_data, &buffer[i], sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&root, &buffer[i], sizeof(root));
    i += sizeof(root);
}
//------------------------------------------------------------
int CMNDRTree::update_rslt(CMNDEntry *_e, float _dist, CMNDEntry *_rslt, 
						float *_key, int _k)
{
	for (int i = 0; i < _k; i ++)
	{
		if (_dist < _key[i])
		{
			for (int j = _k - 1; j > i; j --)
			{
				_rslt[j] = _rslt[j - 1];
				_key[j] = _key[j - 1];
			}
			_rslt[i] = *_e;
			_key[i] = _dist;
			return i;
		}
	}
	error("Error in update_rslt\n", true);
	return -1;
}
//------------------------------------------------------------
void CMNDRTree::write_header(char *buffer)
{
    int i;

    memcpy(buffer, &dimension, sizeof(dimension));
    i = sizeof(dimension);

    memcpy(&buffer[i], &num_of_data, sizeof(num_of_data));
    i += sizeof(num_of_data);

    memcpy(&buffer[i], &num_of_dnodes, sizeof(num_of_dnodes));
    i += sizeof(num_of_dnodes);

    memcpy(&buffer[i], &num_of_inodes, sizeof(num_of_inodes));
    i += sizeof(num_of_inodes);

    memcpy(&buffer[i], &root_is_data, sizeof(root_is_data));
    i += sizeof(root_is_data);

    memcpy(&buffer[i], &root, sizeof(root));
    i += sizeof(root);
}
//------------------------------------------------------------
//---added for valdity region queries-------------------------
// perform a window query mbr to get the query result into in_objs, put the outer objects retrieved into out_objs_so_far
void CMNDRTree::rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far)
{
    load_root();

    root_ptr->rect_win_query(mbr, in_objs, out_objs_so_far);

	delete root_ptr;
	root_ptr = NULL;
}

// perform a window query mbr (excluding the window excr) to get the query result into c_inf_objs
void CMNDRTree::rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs)
{
    load_root();

    root_ptr->rect_win_query(mbr, exclmbr, c_inf_objs);

	delete root_ptr;
	root_ptr = NULL;
}

void CMNDRTree::BFNN(float *_qpt, int _k, CMNDEntry *_rslt)
{
	//first init an array for storing the keys of retrieve objects
	float *key = new float [_k];
	for (int i = 0; i < _k; i ++)
		key[i] = MAXREAL; //initially they are infinity
	//init a heap that stores the non-leaf entries to be accessed-
	Heap *heap = new Heap();
	heap->init(dimension);
	//------------------------------------------------------------
	
	int son = root; //this CMNDEntry is to be visited next
	while (son != -1)
	{
		CMNDRTNode *rtn = new CMNDRTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			float edist = MINDIST(_qpt, rtn->entries[i].bounces, dimension);
			if (rtn->level == 0)
			{
				if (edist < key[_k - 1])
					update_rslt(&(rtn->entries[i]), edist, _rslt, key, _k);
			}
			else
			{
				if (edist<key[_k - 1]) 
					//meaning that edist is valid and we insert it to heap
				{
					HeapEntry *he = new HeapEntry();
					he -> key = edist;
					he -> level = rtn -> level;
					he -> son1 = rtn->entries[i].son;
					heap -> insert(he);
					delete he;
				}
			}
		}
	
		delete rtn;

		//get next CMNDEntry from the heap----------------------------
		HeapEntry *he = new HeapEntry();
		bool again = true;
		while (again)
		{
			again = false;
			if (!heap->remove(he))  //heap is empty
				son = -1;
			else
			{
				if (he->key>key[_k - 1]) //the algorithm terminates
					son = -1;
				else if (he->level == 0) 
						//protection. if you see this message, debug
						error("testing... leaf entries found in heap\n", true);
					 else
						son=he->son1;
			}
		}
		delete he;
		//--------------------------------------------------------
	}

	delete [] key;
	delete heap;
}

// The following code was copied from the implementation of TP-kNN queries from Tony
void CMNDRTree::TPNN_TP(float *_qline, int _k, CMNDEntry *_nn, CMNDEntry *_rslt, float _max_trvl)
{
	float key = _max_trvl; 
	  // the minimum distance that the query point must travel 

//we comment these lines to avoid initing the heap everytime--
//this function is called
//Heap *heap = new Heap();
//heap->init(dimension);
//------------------------------------------------------------
	if (tpheap==NULL)
		error("tpheap is not initialized\n", true);
	tpheap->used=0;
	
	int son = root;
	while (son != -1)
	{
		CMNDRTNode *rtn = new CMNDRTNode(this, son);
		for (int i = 0; i < rtn -> num_entries; i ++)
		{
			//first create an array m2 for e to cal dist----------
			float *m2 = new float [4 * dimension];
			memset(m2, 0, 4 * dimension * sizeof(float));
			memcpy(m2, rtn -> entries[i].bounces, 2 * dimension * sizeof(float));
			//----------------------------------------------------
//testing--------------------------
//if (rtn->entries[i].son==573673 && cnt==84-1)
//	printf("testing...\n");
//---------------------------------

			float edist = MAXREAL;
			if (rtn -> level == 0)
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^leaf node case^^^^^^^^^^^^^^^^^^^^^
			{	
				int eNNsub=-1;
				//if e (i.e., m2) is one of the current NN points-
				//eNNsub stores its subsript in _nn; otherwise, 
				//eNNsub=-1
				for (int j=0; j<_k; j++)
					if (rtn->entries[i].son==_nn[j].son)
					{ eNNsub=j; j=_k;}
				//------------------------------------------------
		
				if (eNNsub==-1 || SEQUENCE_SENSITIVE)
					//namely, if sequence insensitive and e is indeed a NN
					//then we do not handle this CMNDEntry
				{

					float *m1 = new float [4 * dimension];
					//find the NN that leads to the minimum-------
					//influence time
					int nn_num=-1; //the subsript of the NN to be found
					for (int j = 0; j < _k; j ++)
						//for each of the NN found in the 1st step
					{
						bool yesdo=true; //whether to compute
						  //the inflence time of nn[j] and e
						if (j==eNNsub) yesdo=false;
						
						//check if this pair has produced --------
						//influence time before
						for (int l=0; l<PAST_PAIR; l++)
							if (min(_nn[j].son, rtn->entries[i].son)==min(last_pair[l][0], last_pair[l][1]) &&
								max(_nn[j].son, rtn->entries[i].son)==max(last_pair[l][0], last_pair[l][1]))
							{	yesdo=false; l=PAST_PAIR; }
						//----------------------------------------

						if (yesdo)
						{
//these codes use NNinf===========================================
/*
							//first create an array m1 (for nn[j])
							//to compute dist
							memset(m1, 0, 4 * dimension * sizeof(float));
							memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
							//get the influence time of m2-------- 
							//(for CMNDEntry e) with respect to(nn[j])
							float this_inf = NNinf(m1, m2, _qline, dimension); 
							//------------------------------------
*/
//================================================================

//these codes use NNinf2==========================================
							//first create an array m1 (for nn[j])
							//to compute dist
							memset(m1, 0, 4 * dimension * sizeof(float));
							memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
							m1[1]=m1[2]; m2[1]=m2[2];
							//create an arry m3 for _qline----------------
							float *m3=new float[2*dimension];
							m3[0]=_qline[0]; m3[1]=_qline[2];
							m3[2]=_qline[4]; m3[3]=_qline[6];
							//--------------------------------------------
							//get the influence time of m2-------- 
							//(for CMNDEntry e) with respect to(nn[j])
							float this_inf = NNinf2(m1, m2, m3); 
							//------------------------------------
							delete []m3;
//================================================================

							if (this_inf>0 && this_inf<edist)
								//this_inf=0 means that there is another point that has the same distance
								//to the current query position as the current NN. In this implementation,
								//we choose to ignore handling such special cases, which, however, may cause
								//problems for datasets with big cardinality
//							if (this_inf>=0 && this_inf<edist)
							{
								edist=this_inf; nn_num=j;
							}
						}  //END if (yesdo)
					}//END checking all neighbors
					//-------------------------------------------------
					//if (edist<key && edist!=0)
					if (edist<key)
					{
						update_rslt(&(rtn->entries[i]), edist, _rslt, &key, 1);
						_rslt->nn_num=nn_num;
					}
					delete []m1;
				}
			}
//^^^^^^^^^^^^^^^^^^^^^^^^^non-leaf node case^^^^^^^^^^^^^^^^^^^^^
			else
				//Next handle non-leaf node case
			{
				float *m1 = new float [4 * dimension];
				for (int j = 0; j < _k; j ++)
				{
					//first create an array m1 to cal dist--------
					memset(m1, 0, 4 * dimension * sizeof(float));
					memcpy(m1, _nn[j].bounces, 2 * dimension * sizeof(float));
					//--------------------------------------------
					float this_mininf = NNmininf(m1, m2, _qline, dimension);
					if (this_mininf < edist)
						edist = this_mininf;
				}
				delete [] m1;

				if (edist < key)
				{
					HeapEntry *he = new HeapEntry();
					he -> key = edist;
					he -> level = rtn -> level;
					he -> son1 = rtn->entries[i].son;
					tpheap -> insert(he);
					delete he;
				}
			}
			delete [] m2;

		}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	
		delete rtn;

		//get next CMNDEntry from the heap
		bool again = true;
		while (again)
		{
			again = false;
			HeapEntry *he = new HeapEntry();
			if (!tpheap -> remove(he))  //heap is empty
				son = -1;
			else
				if (he -> key > key)
					//the algorithm can terminate
					son = -1;
				else
					son = he -> son1;
			delete he;
		}
	}
//delete heap;
}
