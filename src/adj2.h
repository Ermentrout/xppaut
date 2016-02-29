#ifndef XPPAUT_ADJ2_H
#define XPPAUT_ADJ2_H

#include <stdio.h>

/* --- Types ---*/
typedef struct {
	int here,col0,ncol,colskip;
	int row0,nrow,rowskip;
	float **data;
	char firstcol[11];
} MY_TRANS;

/* --- Data --- */
extern int AdjRange;
extern int FOUR_HERE;

/* --- Functions --- */
void init_trans(void);
void dump_transpose_info(FILE *fp, int f);
int do_transpose(void);
void alloc_h_stuff(void);
void data_back(void);
void make_adj_com(int com);
void new_h_fun(int silent);
void dump_h_stuff(FILE *fp, int f);
void new_adjoint(void);
void do_liapunov(void);
void alloc_liap(int n);
void do_this_liaprun(int i, double p);

#endif /* XPPAUT_ADJ2_H */
