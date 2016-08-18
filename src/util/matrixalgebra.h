#ifndef XPPAUT_UTIL_MATRIXALGEBRA_H
#define XPPAUT_UTIL_MATRIXALGEBRA_H


/* --- Macros --- */
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define SIGN(a,b) ((b)>=0.0 ? fabs(a):-fabs(a))


/* --- Functions --- */
int bandfac(double *a, int ml, int mr, int n);
void bandsol(double *a, double *b, int ml, int mr, int n);
void eigen(int n, double *a, double *ev, double *work, int *ierr);
void get_evec(double *a, double *anew, double *b, double *bp, int n, int maxit, double err, int *ipivot, double eval, int *ierr);
void getjac(double *x, double *y, double *yp, double *xp, double eps, double *dermat, int n);
void getjactrans(double *x,double *y,double *yp,double *xp, double eps, double *d, int n);
void get_the_jac(double t, double *y, double *yp, double *ypnew, double *dfdy, int neq, double eps, double scal);
void sgefa(double *a, int lda, int n, int *ipvt, int *info);
void sgesl(double *a, int lda, int n, int *ipvt, double *b, int job);
double sqr2(double z);

#endif // XPPAUT_UTIL_MATRIXALGEBRA_H
