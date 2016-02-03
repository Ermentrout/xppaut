#ifndef _stiff_h_
#define _stiff_h_


/* stiff.c */
void jacobn(double x, double *y, double *dfdx, double *dermat, double eps, double *work, int n);
int adaptive(double *ystart, int nvar, double *xs, double x2, double eps, double *hguess, double hmin, double *work, int *ier, double epjac, int iflag, int *jstart);
int gadaptive(double *ystart, int nvar, double *xs, double x2, double eps, double *hguess, double hmin, double *work, int *ier, double epjac, int iflag, int *jstart);
int stiff(double y[], double dydx[], int n, double *x, double htry, double eps, double yscal[], double *hdid, double *hnext, double *work, double epjac, int *ier);
int rkqs(double *y, double *dydx, int n, double *x, double htry, double eps, double *yscal, double *hdid, double *hnext, double *work, int *ier);
void rkck(double *y, double *dydx, int n, double x, double h, double *yout, double *yerr, double *work);


#endif

