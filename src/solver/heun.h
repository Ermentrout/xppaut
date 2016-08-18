#ifndef XPPAUT_SOLVER_HEUN_H
#define XPPAUT_SOLVER_HEUN_H

/* --- Functions --- */
int heun(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);

#endif // XPPAUT_SOLVER_HEUN_H
