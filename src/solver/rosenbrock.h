#ifndef XPPAUT_SOLVER_ROSENBROCK_H
#define XPPAUT_SOLVER_ROSENBROCK_H

/* --- Functions --- */
int rb23(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);
int rosen(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);

#endif // XPPAUT_SOLVER_ROSENBROCK_H
