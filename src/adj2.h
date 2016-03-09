#ifndef XPPAUT_ADJ2_H
#define XPPAUT_ADJ2_H

#include <stdio.h>
#include "xpplim.h"

/* --- Types ---*/
typedef struct {
	int here,col0,ncol,colskip;
	int row0,nrow,rowskip;
	float **data;
	char firstcol[11];
} MY_TRANS;

/* --- Data --- */
extern char *coup_string[MAXODE];
extern int AdjRange;
extern int FOUR_HERE;
extern MY_TRANS my_trans;

/* --- Functions --- */
void adj2_init_trans(void);
void adj2_setup_trans(void);
int do_transpose(void);
void alloc_h_stuff(void);
void data_back(void);
void make_adj_com(int com);
void new_h_fun(int silent);
void new_adjoint(void);
void do_liapunov(void);
void alloc_liap(int n);
void do_this_liaprun(int i, double p);

#endif /* XPPAUT_ADJ2_H */
