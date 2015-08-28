#ifndef _volterra2_h_
#define _volterra2_h_


/* volterra2.c */
double ker_val(int in);
void alloc_v_memory(void);
void allocate_volterra(int npts, int flag);
void re_evaluate_kernels(void);
void alloc_kernels(int flag);
void init_sums(double t0, int n, double dt, int i0, int iend, int ishift);
double alpha1n(double mu, double dt, double t, double t0);
double alpbetjn(double mu, double dt, int l);
double betnn(double mu, double dt, double t0, double t);
void get_kn(double *y, double t);
int volterra(double *y, double *t, double dt, int nt, int neq, int *istart, double *work);
int volt_step(double *y, double t, double dt, int neq, double *yg, double *yp, double *yp2, double *ytemp, double *errvec, double *jac);


#endif

