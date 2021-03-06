/*CBBRTNode.cpp
  this file implements class CBBRTNode*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "bbrtree.h"
#include "bbentry.h"
#include "bbrtnode.h"
#include "../blockfile/blk_file.h"
#include "../blockfile/cache.h"
#include "../linlist/linlist.h"
//------------------------------------------------------------
CBBRTNode::CBBRTNode(CBBRTree *rt)
  //use this ructor to create a new node on disk.
{
    char *b;
    int header_size;
    CBBEntry * d;
    int i;

    my_tree = rt;
    dimension = rt->dimension;
    num_entries = 0;
	dirty = TRUE;

    d = new CBBEntry();
	d -> init_entry(dimension, NULL);
    header_size = sizeof(char) + sizeof(int);  // level + num_entries
    capacity = (rt -> file -> get_blocklength() - header_size) / d -> get_size();
    delete d;

    entries = new CBBEntry[capacity];
    for (i = 0; i < capacity; i++)
		entries[i].init_entry(dimension, rt);

	//assign a new block on the disk
    b = new char[rt -> file -> get_blocklength()];
    block = rt -> file -> append_block(b);
    delete [] b;
}
//------------------------------------------------------------
CBBRTNode::CBBRTNode(CBBRTree *rt, int _block)
  //use this ructor to restore a node from the disk.
{
    char *b;
    int header_size;
    CBBEntry * d;
    int i;

    my_tree = rt;
    dimension = rt->dimension;
    num_entries = 0;
	dirty = FALSE;

    d = new CBBEntry();
	d -> init_entry(dimension, NULL);
    header_size = sizeof(char) + sizeof(int);
    capacity = (rt -> file -> get_blocklength() - header_size) / d -> get_size();
    delete d;

    entries = new CBBEntry[capacity];
    for (i = 0; i < capacity; i++)
		entries[i].init_entry(dimension, rt);
    
    block = _block;
    b = new char[rt -> file -> get_blocklength()];
    if (rt -> cache == NULL) // no cache
        rt -> file -> read_block(b, block);
    else
        rt -> cache -> read_block(b, block, rt);

    read_from_buffer(b);
    delete [] b;
}
//------------------------------------------------------------
CBBRTNode::~CBBRTNode()
{
    char *b;

    if (dirty)
    {
        b = new char[my_tree->file->get_blocklength()];
        write_to_buffer(b);

        if (my_tree->cache == NULL) // no cache
            my_tree->file->write_block(b, block);
        else
            my_tree->cache->write_block(b, block, my_tree);

        delete [] b;
    }

    delete [] entries;
}
//------------------------------------------------------------
int CBBRTNode::choose_subtree(float *mbr)
{
    int i, j, follow, minindex, *inside, inside_count, *over;
    float *bmbr, old_o, o, omin, a, amin, f, fmin;

    inside_count = 0;
    inside = new int[num_entries];
    over = new int[num_entries];
    for (i = 0; i < num_entries; i++)
    {
    	switch (entries[i].section(mbr))
    	{
        	case INSIDE:
        	    inside[inside_count++] = i;
        	    break;
        }
    }

    if (inside_count == 1)
        // Case 1: There is exactly one dir_mbr that contains mbr
    	follow = inside[0];
    else if (inside_count > 1)
    // Case 2: There are many dir_mbrs that contain mbr
    // choose the one with the minimum area
    {
    	fmin = MAXREAL;
    	for (i = 0; i < inside_count; i++)
    	{
    	    f = area(dimension, entries[inside[i]].bounces);
    	    if (f < fmin)
    	    {
    	    	minindex = i;
          		fmin = f;
       	    }
       	}
    	follow = inside[minindex];
    }
    else
    // Case 3: There are no dir_mbrs that contain mbr
    // choose the one for which insertion causes the minimun overlap if son_is_data
    // else choose the one for which insertion causes the minimun area enlargement
    {
       	if (level == 1) // son_is_data
    	{
            omin = MAXREAL;
    	    fmin = MAXREAL;
    	    amin = MAXREAL;
    	    for (i = 0; i < num_entries; i++)
    	    {
        		enlarge(dimension, &bmbr, mbr, entries[i].bounces);

        		// calculate area and area enlargement
        		a = area(dimension, entries[i].bounces);
        		f = area(dimension, bmbr) - a;

        		// calculate overlap before enlarging CEntry1_i
        		old_o = o = 0.0;
        		for (j = 0; j < num_entries; j++)
        		{
        		    if (j != i)
        		    {
    			        old_o += overlap(dimension,
    					 entries[i].bounces,
    					 entries[j].bounces);
    			        o += overlap(dimension,
    				     bmbr,
    				     entries[j].bounces);
    		        }
    	        }
    	        o -= old_o;

    	        // is this CBBEntry better than the former optimum ?
    	        if ((o < omin) ||
    		    (o == omin && f < fmin) ||
    		    (o == omin && f == fmin && a < amin))
    	        {
    	       	    minindex = i;
        		    omin = o;
        		    fmin = f;
        		    amin = a;
        	    }
    	        delete [] bmbr;
    	    }
        }
        else // son is not a data node
        {
    	    fmin = MAXREAL;
    	    amin = MAXREAL;
    	    for (i = 0; i < num_entries; i++)
    	    {
    	        enlarge(dimension, &bmbr, mbr, entries[i].bounces);

    	        // calculate area and area enlargement
    	        a = area(dimension, entries[i].bounces);
    	        f = area(dimension, bmbr) - a;

    	        // is this CBBEntry better than the former optimum ?
    	        if ((f < fmin) || (f == fmin && a < amin))
    	        {
    	       	    minindex = i;
    		        fmin = f;
    	            amin = a;
    	        }
	            delete [] bmbr;
	        }
        }

    	follow = minindex;

    	dirty = TRUE;
    }

    delete [] inside;
    delete [] over;

    return follow;
}
//------------------------------------------------------------
R_DELETE CBBRTNode::delete_entry(CBBEntry *e)
{
	CBBRTNode *succ;
	float *tmp;
	if (level > 0)
	{
		if (this == my_tree->root_ptr)
			//i.e. this is the root
		{
			for (int i = 0; i < num_entries; i++)
			{
				tmp = overlapRect(dimension, entries[i].bounces, e -> bounces);
				if (tmp != NULL)
				{
					delete [] tmp;
					succ = entries[i].get_son();
					R_DELETE del_ret;
					del_ret = succ -> delete_entry(e);
					if (del_ret != NOTFOUND)
					{
						switch (del_ret)
						{
						case NORMAL:

							float *mbr;
							mbr = succ -> get_mbr();
							memcpy(entries[i].bounces, mbr, sizeof(float) * 2 * dimension);
							dirty = true;
							delete [] mbr;

							delete entries[i].son_ptr;
							entries[i].son_ptr = NULL;

							return NORMAL;
							break;

						case ERASED:
							delete entries[i].son_ptr;
							entries[i].son_ptr = NULL;

							int j;
							for (j = i; j < num_entries - 1; j++)
								entries[j] = entries[j+1];
							for (j = num_entries - 1; j < capacity; j++)
								entries[j].son_ptr = NULL;

							num_entries--;

							dirty = true;
							return NORMAL;
							break;
						}
					}
				}
			}
			//Not found;
			return NOTFOUND;
		}
		else//is not root and not leaf
		{
			for (int i = 0; i < num_entries; i++)
			{
				tmp = overlapRect(dimension, entries[i].bounces, e -> bounces);
				if (tmp != NULL)
				{
					delete [] tmp;
					succ = entries[i].get_son();
					R_DELETE del_ret;
					del_ret = succ->delete_entry(e);
					if (del_ret != NOTFOUND)
					{
						switch (del_ret)
						{
						case NORMAL:

							float *mbr;
							mbr = succ -> get_mbr();
							memcpy(entries[i].bounces, mbr, sizeof(float) * 2 * dimension);
							dirty = true;
							delete [] mbr;

							entries[i].del_son();

							return NORMAL;
							break;

						case ERASED:

							entries[i].del_son();

							int j;
							for (j = i; j < num_entries - 1; j++)
								entries[j] = entries[j+1];
							for (j = num_entries - 1; j < capacity; j++)
								entries[j].son_ptr = NULL;
							
							num_entries--;

							dirty = true;
							delete succ;

							if (num_entries < (int)ceil(0.4 * capacity))
							{
								for (int j = 0; j < num_entries; j++)
								{
									Linkable *e;
									e = entries[j].gen_Linkable();
									my_tree -> deletelist -> insert(e);
								}

								my_tree -> num_of_inodes --;
								return ERASED;
							}
							else
								return NORMAL;
							break;
						}
					}
				}
			}
		}
	}
	else//it's a leaf
	{
		for (int i = 0; i < num_entries; i++)
		{
			if (entries[i] == (*e))
			{
				my_tree -> num_of_data --;

				for (int j = i; j < num_entries-1; j++)
					entries[j] = entries[j+1];
				
				num_entries--;
				dirty = true;

				if (this != my_tree -> root_ptr && num_entries < (int)ceil(0.4 * capacity))
				{
					for (int k = 0; k < num_entries; k++)
					{
						Linkable *en;
					    en = entries[k].gen_Linkable();
						en -> level = 0;
						my_tree -> deletelist -> insert(en);
					}

					my_tree -> num_of_dnodes --;
					return ERASED;
				}
				else
					return NORMAL;
			}
		}
		return NOTFOUND;
	}
}
//------------------------------------------------------------
void CBBRTNode::enter(CBBEntry *de)
  //note that de will be deleted after being entered.
{
    if (num_entries > (capacity-1))
        error("CBBRTNode::enter: called, but node is full", TRUE);

    entries[num_entries] = *de;

    num_entries++;

	dirty = true;

    de->son_ptr = NULL;
    delete de;
}
//------------------------------------------------------------
bool CBBRTNode::FindLeaf(CBBEntry *e)
{
	CBBRTNode *succ;
	if (level > 0)
	{
		for (int i = 0; i < num_entries; i++)
		{
			float *f;
			f = overlapRect(my_tree -> dimension,
				  entries[i].bounces, e -> bounces);
			if (f != NULL)
			{
				delete [] f;
				succ = entries[i].get_son();
				bool find;
				find = succ->FindLeaf(e);
				entries[i].del_son();
				if (find)
					return true;
			}
		}
		return false;
	}
	else
	{
		for (int i = 0; i < num_entries; i++)
		{
			if (entries[i] == (*e))
				return true;
		}
		return false;
	}
	return false;
}
//------------------------------------------------------------
float* CBBRTNode::get_mbr()
{
    int i, j;
    float *mbr;

    mbr = new float[2*dimension];
    for (i = 0; i < 2*dimension; i ++)
        mbr[i] = entries[0].bounces[i];

    for (j = 1; j < num_entries; j++)
    {
    	for (i = 0; i < 2*dimension; i += 2)
    	{
    	    mbr[i]   = min(mbr[i],   entries[j].bounces[i]);
    	    mbr[i+1] = max(mbr[i+1], entries[j].bounces[i+1]);
        }
    }

    return mbr;
}
//------------------------------------------------------------
int CBBRTNode::get_num_of_data()
{
    int i, sum;
    CBBRTNode* succ;

    if (level == 0)
        return num_entries;

    sum = 0;
    for (i = 0; i < num_entries ; i++)
    {
        succ = entries[i].get_son();
        sum += succ->get_num_of_data();
		entries[i].del_son();
    }

    return sum;
}
//------------------------------------------------------------
R_OVERFLOW CBBRTNode::insert(CBBEntry *d, CBBRTNode **sn)
{
    int follow;
    CBBRTNode *succ, *new_succ;
    CBBRTNode *brother;
    CBBEntry *de;
    R_OVERFLOW ret;
    float *mbr,*nmbr;

    int i, last_cand;
    float *center;
    SortMbr *sm;
    CBBEntry *new_entries;

    if (level > 0) // direcrtory node
    {
	  if (level > d -> level)
	  {
        follow = choose_subtree(d -> bounces);

        succ = entries[follow].get_son();

        ret = succ -> insert(d, &new_succ);
    
        mbr = succ -> get_mbr();
        memcpy(entries[follow].bounces, mbr, sizeof(float) * 2 * dimension);
        delete [] mbr;

		entries[follow].del_son();

        if (ret == SPLIT)
        // node has split into itself and *new_succ
        {
            if (num_entries == capacity)
         	    error("CBBRTNode::insert: maximum capacity violation", TRUE);

            de = new CBBEntry(dimension, my_tree);
    	    nmbr = new_succ -> get_mbr();
            memcpy(de -> bounces, nmbr, 2 * dimension * sizeof(float));
    	    delete [] nmbr;
            de -> son = new_succ -> block;
			delete new_succ;
            de -> son_ptr = NULL;
            enter(de);

            if (num_entries == (capacity - 1))
            {
        	    brother = new CBBRTNode(my_tree);
        	    my_tree -> num_of_inodes++;
        	    brother -> level = level;
        	    split(brother);
                *sn = brother;
                ret = SPLIT;
        	}
            else
          	    ret = NONE;
        }
        dirty = TRUE;

        return ret;
	  }
	  else //level==d->level
	  {
		  enter(d);    //note that d will be deleted on return
		    
		  if (num_entries == (capacity - 1))
            // maximun no of entries --> Split
            // this happens already if the node is nearly filled
            // for the algorithms are more easy then
		  {
        	brother = new CBBRTNode(my_tree);
        	my_tree -> num_of_inodes++;
        	brother -> level = level;
        	split(brother);
            *sn = brother;
            ret = SPLIT;
		  }
          else
          	ret = NONE;

		  dirty=true;
		  return ret;
	  }	
    }
    else // data (leaf) node
    {
        if (num_entries == capacity)
        	error("RTDataNode::insert: maximum capacity violation", TRUE);

        enter(d);

        dirty = TRUE;

        if (num_entries == (capacity - 1))
        // maximum # of entries --> Split
        // this happens already if the node is nearly filled
        // for the algorithms are more easy then
        {
            if (my_tree->re_level[0] == FALSE && my_tree -> root_ptr -> level != level)
    	    // there was no reinsert on level 0 during this insertion
			//-------------------------------------------------------
			//Here I changed the condition as if it is already root, no need
			//to reinsert.  Split directly in this case
			//-----------By TAO Yufei
            {
                // calculate center of page
                mbr = get_mbr();
                center = new float[dimension];
                for (i = 0; i < dimension; i++)
                     center[i] = (mbr[2*i] + mbr[2*i+1]) / 2.0;

                new_entries = new CBBEntry[capacity];

				for (i = 0; i < capacity; i ++)
					new_entries[i].init_entry(dimension, my_tree);

        	    sm = new SortMbr[num_entries];
        	    for (i = 0; i < num_entries; i++)
        	    {
            		sm[i].index = i;
            		sm[i].dimension = dimension;
            		sm[i].mbr = entries[i].bounces;
            		sm[i].center = center;
                }

                qsort(sm, num_entries, sizeof(SortMbr), sort_center_mbr);

                last_cand = (int) ((float)num_entries * 0.30);

                // copy the nearest candidates to new array
                for (i = 0; i < num_entries - last_cand; i++)
    	            new_entries[i] = entries[sm[i].index];

                // insert candidates into reinsertion list
                for (; i < num_entries; i++)
                {
					Linkable *nd = entries[sm[i].index].gen_Linkable();
                    my_tree -> re_data_cands -> insert(nd);
                }

                // free and copy data array
                delete [] entries;
        	    entries = new_entries;
				
        	    delete sm;
        	    delete [] mbr;
        	    delete [] center;
        	    my_tree -> re_level[0] = TRUE;

        	    // correct # of entries
        	    num_entries -= last_cand;

        	    // must write page
        	    dirty = TRUE;

                return REINSERT;
        	}
           	else  //there has been reinsertion on this level
           	{
        	    *sn = new CBBRTNode(my_tree);
        	    (*sn) -> level = level;
        	    my_tree -> num_of_dnodes++;
        	    split((CBBRTNode *) *sn);
    	    }
    	    return SPLIT;
        }
        else
            return NONE;
    }
}
//------------------------------------------------------------
void CBBRTNode::print()
{
    int i;

	printf("level %d  Block: %d\n", level, block);
	
    for (i = 0; i < num_entries ; i++)
    {
        printf("(%4.1lf, %4.1lf, %4.1lf, %4.1lf)\n",
	       entries[i].bounces[0],
	       entries[i].bounces[1],
	       entries[i].bounces[2],
	       entries[i].bounces[3]);
    }
}
//------------------------------------------------------------
void CBBRTNode::rangeQuery(float *mbr, SortedLinList *res)
{
    int i, n;
    SECTION s;
    CBBRTNode *succ;

    n = num_entries;
    for (i = 0; i < n; i++)
    {
        s = entries[i].section(mbr);
        if (s == INSIDE || s == OVERLAP)
        {
            if (level == 0)
            {
                Linkable *copy;
				copy = entries[i].gen_Linkable();
        		res -> insert(copy);
            }
            else
            {
                succ = entries[i].get_son();
                succ -> rangeQuery(mbr, res);
				entries[i].del_son();
            }
        }
    }
}
//------------------------------------------------------------
void CBBRTNode::read_from_buffer(char *buffer)
{
    int i, j, s;

    // Level
    memcpy(&level, buffer, sizeof(char));
    j = sizeof(char);

    // num_entries
    memcpy(&num_entries, &(buffer[j]), sizeof(int));
    j += sizeof(int);

    s = entries[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
    	entries[i].read_from_buffer(&buffer[j]);
    	j += s;
    }
}
//------------------------------------------------------------
int CBBRTNode::split(float **mbr, int **distribution)
{
    bool lu;
    int i, j, k, l, s, n, m1, dist, split_axis;
    SortMbr *sml, *smu;
    float minmarg, marg, minover, mindead, dead, over, *rxmbr, *rymbr;

    n = num_entries;

    m1 = (int) ceil((float)n * 0.40);

    sml = new SortMbr[n];
    smu = new SortMbr[n];
    rxmbr = new float[2*dimension];
    rymbr = new float[2*dimension];

    // choose split axis
    minmarg = MAXREAL;
    for (i = 0; i < dimension; i++)
    // for each axis
    {
        for (j = 0; j < n; j++)
        {
            sml[j].index = smu[j].index = j;
            sml[j].dimension = smu[j].dimension = i;
            sml[j].mbr = smu[j].mbr = mbr[j];
        }

        // Sort by lower and upper value perpendicular axis_i
      	qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
        qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

        marg = 0.0;
        // for all possible distributions of sml
        for (k = 0; k < n - 2 * m1 + 1; k++)
        {
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for (l = 0; l < m1 + k; l++)
            {
				for (s = 0; s < 2*dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
				}
			}
			marg += margin(dimension, rxmbr);

			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for (; l < n; l++)
            {
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);
        }

        // for all possible distributions of smu
       	for (k = 0; k < n - 2 * m1 + 1; k++)
        {
            // now calculate margin of R1
			// initialize mbr of R1
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for (l = 0; l < m1+k; l++)
            {
                // calculate mbr of R1
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);

            // now calculate margin of R2
			// initialize mbr of R2
			for (s = 0; s < 2 * dimension; s += 2)
			{
				rxmbr[s] =    MAXREAL;
				rxmbr[s+1] = -MAXREAL;
			}
            for (; l < n; l++)
            {
                // calculate mbr of R1
				for (s = 0; s < 2 * dimension; s += 2)
				{
					rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
					rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
				}
            }
			marg += margin(dimension, rxmbr);
        }

        if (marg < minmarg)
        {
            split_axis = i;
            minmarg = marg;
        }
    }

    // choose best distribution for split axis
    for (j = 0; j < n; j++)
    {
		sml[j].index = smu[j].index = j;
		sml[j].dimension = smu[j].dimension = split_axis;
		sml[j].mbr = smu[j].mbr = mbr[j];
    }

    // Sort by lower and upper value perpendicular split axis
    qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
    qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

    minover = MAXREAL;
    mindead = MAXREAL;
    // for all possible distributions of sml and snu
    for (k = 0; k < n - 2 * m1 + 1; k++)
    {
        dead = 0.0;
		for (s = 0; s < 2 * dimension; s += 2)
		{
			rxmbr[s] =    MAXREAL;
			rxmbr[s+1] = -MAXREAL;
		}
		for (l = 0; l < m1 + k; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
				rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
			}
			dead -= area(dimension, sml[l].mbr);
		}
        dead += area(dimension, rxmbr);
		  //**************note**************
		  //this does not compute the dead space for all the cases.  some overlapping
		  //area may be subtrated twice.
		  //********************************

		for (s = 0; s < 2 * dimension; s += 2)
		{
			rymbr[s] =    MAXREAL;
       		rymbr[s+1] = -MAXREAL;
		}
		for (; l < n; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rymbr[s] =   min(rymbr[s],   sml[l].mbr[s]);
				rymbr[s+1] = max(rymbr[s+1], sml[l].mbr[s+1]);
			}
			dead -= area(dimension, sml[l].mbr);
		}
        dead += area(dimension, rymbr);

		over = overlap(dimension, rxmbr, rymbr);

        if ((over < minover) || (over == minover) && dead < mindead)
        {
            minover = over;
            mindead = dead;
            dist = m1+k;
            lu = TRUE;
        }

		//Now we do the same thing for smu
        dead = 0.0;
		for (s = 0; s < 2*dimension; s += 2)
		{
			rxmbr[s] =    MAXREAL;
			rxmbr[s+1] = -MAXREAL;
		}
		for (l = 0; l < m1+k; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
				rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
			}
			dead -= area(dimension, smu[l].mbr);
		}
        dead += area(dimension, rxmbr);

		for (s = 0; s < 2*dimension; s += 2)
		{
			rymbr[s] =    MAXREAL;
			rymbr[s+1] = -MAXREAL;
		}
		for (; l < n; l++)
		{
			for (s = 0; s < 2*dimension; s += 2)
			{
				rymbr[s] =   min(rymbr[s],   smu[l].mbr[s]);
				rymbr[s+1] = max(rymbr[s+1], smu[l].mbr[s+1]);
			}
			dead -= area(dimension, smu[l].mbr);
		}
		//correcting errors
        dead += area(dimension, rymbr);

		over = overlap(dimension, rxmbr, rymbr);

        if ((over < minover) || (over == minover) && dead < mindead)
        {
            minover = over;
            mindead = dead;
            dist = m1+k;
            lu = FALSE;
        }
    }

    // calculate best distribution
	// the array distribution is deleted in split(CBBRTNode *sn);
    *distribution = new int[n];
    for (i = 0; i < n; i++)
    {
        if (lu)
            (*distribution)[i] = sml[i].index;
        else
            (*distribution)[i] = smu[i].index;
    }

    delete [] sml;
    delete [] smu;
    delete [] rxmbr;
    delete [] rymbr;

    return dist;
}
//------------------------------------------------------------
void CBBRTNode::split(CBBRTNode *sn)
{
    int i, *distribution, dist, n;
    float **mbr_array;
    CBBEntry *new_entries1, *new_entries2;

    n = num_entries;

    mbr_array = new floatptr[n];
    for (i = 0; i < n; i++)
       	mbr_array[i] = entries[i].bounces;

    dist = split(mbr_array, &distribution);

    new_entries1 = new CBBEntry[capacity];
    new_entries2 = new CBBEntry[capacity];

	for (i = 0; i < capacity; i ++)
	{
		new_entries1[i].init_entry(dimension, my_tree);
		new_entries2[i].init_entry(dimension, my_tree);
	}

    for (i = 0; i < dist; i++)
       	new_entries1[i] = entries[distribution[i]];

    for (i = dist; i < n; i++)
       	new_entries2[i-dist] = entries[distribution[i]];

    for (i = 0; i < n; i++)
    {
       	entries[i].son_ptr = NULL;
       	sn->entries[i].son_ptr = NULL;
    }
    delete [] entries;
    delete [] sn->entries;

    entries = new_entries1;
    sn->entries = new_entries2;

    num_entries = dist;
    sn->num_entries = n - dist;

    delete [] mbr_array;
	delete [] distribution;
}
//------------------------------------------------------------
void CBBRTNode::write_to_buffer(char *buffer)
{
    int i, j, s;

    // Level
    memcpy(buffer, &level, sizeof(char));
    j = sizeof(char);

    // num_entries
    memcpy(&buffer[j], &num_entries, sizeof(int));
    j += sizeof(int);

    s = entries[0].get_size();
    for (i = 0; i < num_entries; i++)
    {
    	entries[i].write_to_buffer(&buffer[j]);
       	j += s;
    }
}

//------------------------------------------------------------
void CBBRTNode::rect_win_query(float *mbr, LinList *in_objs, LinList *out_objs_so_far)
{
    int i, n;
    SECTION s;
    CBBRTNode *succ;

    n = num_entries;
    for (i = 0; i < n; i++)
    {
		s = section_new(2, mbr, entries[i].bounces);

		if (level == 0)
        {
			if (s == INSIDE || s == OVERLAP)
			{
				Linkable *copy;
				copy = entries[i].gen_Linkable();
        		in_objs -> insert(copy);
	        }

/*  we remarked this because we assume there is no buffer available to put the outer objects found so far
			else
			{
				Linkable *copy;
				copy = entries[i].gen_Linkable();
        		out_objs_so_far -> insert(copy);
			}
*/		
		}
        else
        {
			if (s == INSIDE || s == OVERLAP)
			{
				succ = entries[i].get_son();
				succ -> rect_win_query(mbr, in_objs, out_objs_so_far);
				entries[i].del_son();
	        }
        }
    }
}

//------------------------------------------------------------
void CBBRTNode::rect_win_query(float *mbr, float *exclmbr, LinList *c_inf_objs)
{
    int i, n;
    SECTION s, s1;
    CBBRTNode *succ;

    n = num_entries;
    for (i = 0; i < n; i++)
    {
		s = section_new(2, mbr, entries[i].bounces);
		s1 = section_new(2, exclmbr, entries[i].bounces);

		if (level == 0)
        {
			if ((s == INSIDE || s == OVERLAP) && s1 != INSIDE && s1 != OVERLAP)
			{
				Linkable *copy;
				copy = entries[i].gen_Linkable();
        		c_inf_objs -> insert(copy);
	        }
		}

/* we remarked this because we assume there is no buffer available to put the outer objects found so far
		else if (level == 1)
        {
			if ((s == INSIDE || s == OVERLAP) && s1 != INSIDE && s1 != OVERLAP)
			{
				succ = entries[i].get_son();
				succ -> rect_win_query(mbr, exclmbr, c_inf_objs);
				entries[i].del_son();
	        }
        }
*/
		else 
        {
			if ((s == INSIDE || s == OVERLAP) && s1 != INSIDE)
			{
				succ = entries[i].get_son();
				succ -> rect_win_query(mbr, exclmbr, c_inf_objs);
				entries[i].del_son();
	        }
        }
    }
}
