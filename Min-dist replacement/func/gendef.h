#ifndef __GENERAL_DEFINITION
#define __GENERAL_DEFINITION

#include <stdio.h>
#include <ctype.h>

//----Valdity region query----------------------------------
#define MAX_REAL 1e+20;


//----BlockFile, CachedBlockFile, Cache---------------------
#define BFHEAD_LENGTH (sizeof(int)*2)    //file header size

#define TRUE 1
#define FALSE 0

#define SEEK_CUR 1
#define SEEK_SET 0
#define SEEK_END 2

typedef char Block[];
//-------------------All-------------------------------------
#define MAXREAL         1e20
#define FLOATZERO       1e-8
#define MAX_DIMENSION   256

#define DIMENSION 2

#define TRUE 1
#define FALSE 0

#define min(a, b) (((a) < (b))? (a) : (b) )
#define max(a, b) (((a) > (b))? (a) : (b) )
//-------------------Class and other data types--------------
class BlockFile;  //for BlockFile definition
class Cache;
class Cacheable   //inherit this class if you wish to design an external
                  //memory structure that can be cached
{
public:
	BlockFile *file;
	Cache *cache;
};
  //==========================================================
class CmdIntrpr  //this is the class of command interpretor.  for a new rtree decescendent
                  //inherit this class to obtain tailored command interpretor
{
public:
	int cnfrm_cmd(char *_msg)
	{ char c = ' ';
	  while (c != 'N' && c != 'Y')
	  { printf("%s (y/n)?", _msg);
	    c = getchar(); 
		char tmp;
		while ((tmp = getchar()) != '\n');
		c = toupper(c); }
	  if (c == 'N') return 0; else return 1; }
  
	void get_cmd(char *_msg, char *_cmd)
	{ printf("%s", _msg);  
	  char *c = _cmd;
	  while (((*c) = getchar()) != '\n')
	    c++;
	  *c = '\0'; } 

	virtual bool build_tree(char *_tree_fname, char *_data_fname, int _b_len, int _dim, int _csize) = 0;
	virtual void free_tree() = 0;
	virtual int qry_sngle(float *_mbr, int *_io_count) = 0;
	/*
	virtual int qry_wrkld(char *_wrkld_fname, int *_io_count, bool _display) = 0;
	*/
	virtual void run() = 0;
	virtual void version() = 0;
};
  //==========================================================
enum SECTION {OVERLAP, INSIDE, S_NONE};
enum R_OVERFLOW {SPLIT, REINSERT, NONE};
enum R_DELETE {NOTFOUND,NORMAL,ERASED};
typedef float *floatptr;
  //==========================================================
struct SortMbr
{
    int dimension;
    float *mbr;
    float *center;
    int index;
};

struct BranchList
{
    int entry_number;
    float mindist;
    float minmaxdist;
    bool section;
};

//-----Global Functions--------------------------------------
void error(char *_errmsg, bool _terminate);

float area(int dimension, float *mbr);
float margin(int dimension, float *mbr);
float overlap(int dimension, float *r1, float *r2);
float* overlapRect(int dimension, float *r1, float *r2);
float objectDIST(float *p1, float *p2, int dim);
float MINMAXDIST(float *Point, float *bounces);
float MINDIST(float *Point, float *bounces, int dim);

bool inside(float &p, float &lb, float &ub);
void enlarge(int dimension, float **mbr, float *r1, float *r2);
bool is_inside(int dimension, float *p, float *mbr);
int pruneBrunchList(float *nearest_distanz, const void *activebrunchList, 
		    int n);
bool section(int dimension, float *mbr1, float *mbr2);
bool section_c(int dimension, float *mbr1, float *center, float radius);

int sort_lower_mbr(const void *d1, const void *d2);
int sort_upper_mbr(const void *d1, const void *d2);
int sort_center_mbr(const void *d1, const void *d2);
int sortmindist(const void *element1, const void *element2);

//---added for validity region---------------------------------
float distinit(float *p, float *e);
float distcont(float *l0, float *l1, float *e);
SECTION section_new(int dimension, float *p, float *q);

#ifdef UNIX
void strupr(char *_msg);
#endif
//-----------------------------------------------------------
#endif