#ifndef _gear_h_
#define _gear_h_

void silent_fixpt(double *x, double eps, double err, double big, int maxit, int n, double *er, double *em, int *ierr);
void do_sing(double *x, double eps, double err, double big, int maxit, int n, int *ierr, float *stabinfo);
void do_sing_info(double *x, double eps, double err, double big, int maxit, int n, double *er, double *em, int *ierr);


void shoot_this_now();
void pr_evec(double *x, double *ev, int n, int pr, double eval,int type);
void get_complex_evec(double *m, double evr, double evm, double *br, double *bm, int n, int maxit, double err, int *ierr);
void get_evec(double *a, double *anew, double *b, double *bp, int n, int maxit, double err, int *ipivot, double eval, int *ierr);
void eigen(int n, double *a, double *ev, double *work, int *ierr);
void hqr(int n, int low, int igh, double *h, double *ev, int *ierr);
void orthes(int n, int low, int igh, double *a, double *ort);
double sign(double x, double y);
int imin(int x, int y);
double amax(double u, double v);
void getjac(double *x, double *y, double *yp, double *xp, double eps, double *dermat, int n);
void rooter(double *x, double err, double eps, double big, double *work, int *ierr, int maxit, int n);
double sqr2(double z);
int gear(int n, double *t, double tout, double *y, double hmin, double hmax, double eps, int mf, double *error, int *kflag, int *jstart, double *work, int *iwork);
int ggear(int n, double *t, double tout, double *y, double hmin, double hmax, double eps, int mf, double *error, int *kflag, int *jstart, double *work, int *iwork);
double sgnum(double x, double y);
double Max(double x, double y);
double Min(double x, double y);
void sgefa(double *a, int lda, int n, int *ipvt, int *info);
void  sgesl(double *a, int lda, int n, int *ipvt, double *b, int job);
void saxpy(int n, double sa, double *sx, int incx, double *sy, int incy);
int isamax(int n, double *sx, int incx);
double sdot(int n, double *sx, int incx, double *sy, int incy);
void sscal(int n, double sa, double *sx, int incx);


#endif
