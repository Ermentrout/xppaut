#ifndef XPPAUT_STIFF_H
#define XPPAUT_STIFF_H

/* --- Functions --- */
int adaptive(double *ystart, int nvar, double *xs, double x2, double eps, double *hguess, double hmin, double *work, int *ier, double epjac, int iflag, int *jstart);
int gadaptive(double *ystart, int nvar, double *xs, double x2, double eps, double *hguess, double hmin, double *work, int *ier, double epjac, int iflag, int *jstart);

#endif /* XPPAUT_STIFF_H */
