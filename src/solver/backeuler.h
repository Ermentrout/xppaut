#ifndef XPPAUT_SOLVER_BACKEULER_H
#define XPPAUT_SOLVER_BACKEULER_H

/* --- Functions --- */
int backwards_euler(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);

#endif // XPPAUT_SOLVER_BACKEULER_H
