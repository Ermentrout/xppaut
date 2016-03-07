#ifndef XPPAUT_ODESOL2_H
#define XPPAUT_ODESOL2_H

/* --- Functions --- */
int symplect3(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int discrete(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int bak_euler(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int one_bak_step(double *y, double *t, double dt, int neq, double *yg, double *yp, double *yp2, double *ytemp, double *errvec, double *jac, int *istart);
void one_step_discrete(double *y, double dt, double *yp, int neq, double *t);
void one_step_symp(double *y, double h, double *f, int n, double *t);
void one_step_euler(double *y, double dt, double *yp, int neq, double *t);
void one_step_rk4(double *y, double dt, double *yval[3], int neq, double *tim);
void one_step_heun(double *y, double dt, double *yval[2], int neq, double *tim);
int euler(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int mod_euler(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int rung_kut(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work);
int adams(double *y, double *tim, double dt, int nstep, int neq, int *ist, double *work);
int rb23(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);
int rosen(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);

#endif /* XPPAUT_ODESOL2_H */
