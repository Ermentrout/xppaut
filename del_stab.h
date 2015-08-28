#ifndef _del_stab_h_
#define _del_stab_h_

typedef struct{
  double r,i;
}COMPLEX;

/* del_stab.c */
void do_delay_sing(double *x, double eps, double err, double big, int maxit, int n, int *ierr, float *stabinfo);
COMPLEX csum(COMPLEX z, COMPLEX w);
COMPLEX cdif(COMPLEX z, COMPLEX w);
COMPLEX cmlt(COMPLEX z, COMPLEX w);
COMPLEX cdiv(COMPLEX z, COMPLEX w);
COMPLEX cexp2(COMPLEX z);
void switch_rows(COMPLEX *z, int i1, int i2, int n);
COMPLEX rtoc(double x, double y);
void cprintn(COMPLEX z);
void cprint(COMPLEX z);
void cprintarr(COMPLEX *z, int n, int m);
double c_abs(COMPLEX z);
COMPLEX cdeterm(COMPLEX *z, int n);
COMPLEX cxdeterm(COMPLEX *z, int n);
void make_z(COMPLEX *z, double *delay, int n, int m, double *coef, COMPLEX lambda);
int find_positive_root(double *coef, double *delay, int n, int m, double rad, double err, double eps, double big, int maxit, double *rr);
void process_root(double real, double im);
double get_arg(double *delay, double *coef, int m, int n, COMPLEX lambda);
int test_sign(double old, double new);
int plot_args(double *coef, double *delay, int n, int m, int npts, double almax, double wmax);


#endif
