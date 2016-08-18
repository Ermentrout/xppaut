#ifndef XPPAUT_FIND_FIXEDPOINT_H
#define XPPAUT_FIND_FIXEDPOINT_H

#include "xpplim.h"


/* --- Data --- */
extern double ShootIC[8][MAXODE];
extern int ShootICFlag;
extern int ShootIndex;
extern int ShootType[8];


/* --- Functions --- */
void do_sing(double *x, double eps, double err, double big, int maxit, int n, int *ierr, float *stabinfo);
void do_sing_info(double *x, double eps, double err, double big, int maxit, int n, double *er, double *em, int *ierr);
void rooter(double *x, double err, double eps, double big, double *work, int *ierr, int maxit, int n);

#endif // XPPAUT_FIND_FIXEDPOINT_H
