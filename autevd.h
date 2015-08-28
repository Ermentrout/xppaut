#ifndef _autevd_h_
#define _autevd_h_

#include "autlim.h"

/* autevd.c */

typedef struct {double r,i;} dcomplex;

typedef struct {
  int pt,br;
  double evr[NAUTO],evi[NAUTO];
} EIGVAL;

double dreal_(dcomplex *z);
void send_eigen(int ibr, int ntot, int n, dcomplex *ev);
void send_mult(int ibr, int ntot, int n, dcomplex *ev);
int get_bif_type(int ibr, int ntot, int lab);
void addbif_(int *ibr, int *ntot, int *itp, int *lab, int *npar, double *a, double *uhigh, double *ulow, double *u0, double *ubar, int *ndim);
double etime_(double *z);
int eigrf_(double *a, int *n, int *m, dcomplex *ecv, double *work, int *ier);
void init_auto(int ndim, int nbc, int ips, int irs, int ilp, int ntst, int isp, int isw, int nmx, int npr, double ds, double dsmin, double dsmax, double rl0, double rl1, double a0, double a1, int ip1, int ip2, int ip3, int ip4, int ip5, int nuzr, double epsl, double epsu, double epss, int ncol);


#endif
