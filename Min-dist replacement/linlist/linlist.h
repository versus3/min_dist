#ifndef __LINLIST
#define __LINLIST

#include <stdio.h>
////////////////////////////////////////////////////////////////////////
// LinList  (SLink)
////////////////////////////////////////////////////////////////////////

class RTree;

struct Linkable
{
public:
	int son;
	int dimension;
	int level;
	float *bounces;
	float distanz;

	// qjz
	float dd;
	int data_entry_num;

	float dd2;
	// qjz

	Linkable(int dim)
	{ dimension = dim;
	  bounces = new float[2 * dim];
    }

	~Linkable()
	{ 
		delete [] bounces;
	}

	void clone(Linkable *old);
};

struct SLink
{
    Linkable *d;          // Pointer to member data
    SLink *next;          // Pointer to next element
    SLink *prev;          // Pointer to previous element

    SLink();
    ~SLink();
};

////////////////////////////////////////////////////////////////////////
// LinList
////////////////////////////////////////////////////////////////////////


class LinList
{
protected:
    SLink *first;         // Rootzeiger des Datenbestands
    SLink *last;          // Pointer to last element
    int anz;              // Number of used items in the list
    SLink *akt;           // points to the current item
    int akt_index;              // Index of the last fetched element
public:
    LinList();
    virtual ~LinList();
    int get_num()               // gibt Anzahl der im Index belegten Elements
        { return anz; }         // zurueck

    void check();               // checked consistency of the list
    void print();
	void printnew();

    void insert(Linkable *f);       // hangs on the front of an item to the list
    bool erase();               // deletes the current item from the list

    Linkable * get(int i);          // provides i-th Element
    Linkable * get_first();         // returns the first element in the index
    Linkable * get_last();          // returns the last element in the index
	Linkable * get_lastnew();		//added by Shirley, get the last linkable but the current poistion is unchanged
	Linkable * get_next();          // returns the next element in the index
    Linkable * get_prev();          // returns the previous element in the index
};


////////////////////////////////////////////////////////////////////////
// SortedLinList
////////////////////////////////////////////////////////////////////////


class SortedLinList : public LinList
{
    bool increasing;
public:
    SortedLinList();

    void set_sorting(bool _increasing); // wenn increasing gleich TRUE, wird
                                // aufsteigend einsortiert
                                // DIESE FUNKTION MUSS VOR DEM ERSTEN EINFUEGEN
                                // GERUFEN WERDEN !!!!!!!!!!!
    void insert(Linkable *f);       // fuegt ein Element durch direktes Einfuegen ein

    void sort(bool _increasing);// sortiert die Liste
};

#endif  // __LINLIST
