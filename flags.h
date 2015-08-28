#ifndef _flags_h_
#define _flags_h_


/* flags.c */
int add_global(char *cond, int sign, char *rest);
void show_flags(void);
int compile_flags(void);
int one_flag_step(double *yold, double *ynew, int *istart, double told, double *tnew, int neq, double *s);
int one_flag_step_symp(double *y, double dt, double *work, int neq, double *tim, int *istart);
int one_flag_step_euler(double *y, double dt, double *work, int neq, double *tim, int *istart);
int one_flag_step_discrete(double *y, double dt, double *work, int neq, double *tim, int *istart);
int one_flag_step_heun(double *y, double dt, double *yval[2], int neq, double *tim, int *istart);
int one_flag_step_rk4(double *y, double dt, double *yval[3], int neq, double *tim, int *istart);
void printflaginfo(void);
int one_flag_step_gear(int neq, double *t, double tout, double *y, double hmin, double hmax, double eps, int mf, double *error, int *kflag, int *jstart, double *work, int *iwork);
int one_flag_step_rosen(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);
int one_flag_step_dp(int *istart, double *y, double *t, int n, double tout, double *tol, double *atol, int flag, int *kflag);
int one_flag_step_adap(double *y, int neq, double *t, double tout, double eps, double *hguess, double hmin, double *work, int *ier, double epjac, int iflag, int *jstart);
int one_flag_step_backeul(double *y, double *t, double dt, int neq, double *yg, double *yp, double *yp2, double *ytemp, double *errvec, double *jac, int *istart);
int one_flag_step_cvode(int *command,double *y,double *t,int n,double tout,int *kflag,double *atol,double *rtol);

#endif

