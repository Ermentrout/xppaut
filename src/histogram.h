#ifndef XPPAUT_HISTOGRAM_H
#define XPPAUT_HISTOGRAM_H

#include <stdio.h>
#include "xpplim.h"

/* --- Types --- */
typedef struct {
	int nbins,nbins2,type,col,col2,fftc;
	double xlo,xhi;
	double ylo,yhi;
	char cond[80];
} HIST_INFO;

/* --- Data --- */
extern double MyData[MAXODE];

extern int post_process;
extern int spec_col;
extern int spec_col2;
extern int spec_wid;
extern int spec_win;

extern HIST_INFO hist_inf;

/* --- Functions --- */
void do_stochast_com(int i);
void post_process_stuff(void);

#endif /* XPPAUT_HISTOGRAM_H */
