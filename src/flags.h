#ifndef XPPAUT_FLAGS_H
#define XPPAUT_FLAGS_H

/* --- Macros --- */
#define MAX_EVENTS 20 /*  this is the maximum number of events per flag */

/* --- Types --- */
typedef struct {
  double f0,f1;
  double tstar;
  int lhs[MAX_EVENTS];
  double vrhs[MAX_EVENTS];
  char lhsname[MAX_EVENTS][11];
  char *rhs[MAX_EVENTS];
  int *comrhs[MAX_EVENTS];
  char *cond;
  int *comcond;
  int sign,nevents;
  int hit,type[MAX_EVENTS];
  int anypars;
  int nointerp;
} FLAG;

/* --- Data --- */
extern int NFlags;
extern double STOL;

/* --- Functions --- */
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

#endif /* XPPAUT_FLAGS_H */

