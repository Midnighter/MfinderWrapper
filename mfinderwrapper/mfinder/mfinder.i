%module mfinder

%{
    #include "common.h"
    #include "results.h"
    extern Res_tbl* subgraphs_interface(int* edges, int num_edges, int mtf_sz,
            int max);
%}

%inline %{
    void get_members(void* p, int* values,int mtf_sz) {
        unsigned i;
        Member* mems = (Member*)p;
        for (i = 0; i < mtf_sz; ++i) {
            values[i] = mems[i].node;
        }
    }

    Motif* get_motif(void* p) {
        return ((Motif*)p);
    }
%}

#define REAL_NET 1;
#define RAND_NET 2;
typedef long long int64;
extern Res_tbl* subgraphs_interface(int* edges, int num_edges, int mtf_sz,
        int max);

extern void res_tbl_mem_free(Res_tbl* res);

%include <carrays.i>
%array_class(int, IntArray);

typedef struct {
    list64 *real;
    list64 **rand_arr;
} Res_tbl;

typedef struct {
    int64 id;
    double count;
    double prob_count; //sigma (Ri/pi)
    double conc; //concentration in mili motifs (counts/total motif counts)*1000
    int hits; //hits num in prob approach
    list *members; //for uniqueness
    list *all_members; //for list of all members of this motif
    int numberOfSelfEdges; // saves the number of self edges in the motif
    double conv_grade; //convergness grade
} Motif;

typedef struct {
    int node;
    unsigned char role;
} Member;

typedef struct{
	int size;
	list64_item *l;
}list64;

typedef struct _list64_item{
	int64 val;
	void *p;
	struct _list64_item *next;
}list64_item;

typedef struct{
	int size;
	list_item *l;
}list;

typedef struct _list_item{
	int val;
	void *p;
	struct _list_item *next;
}list_item;
