#ifndef XPPAUT_VOLTERRA2_H
#define XPPAUT_VOLTERRA2_H

/* --- Data --- */
extern int AutoEvaluate;

/* --- Functions --- */
void alloc_kernels(int flag);
void alloc_v_memory(void);
void allocate_volterra(int npts, int flag);
double ker_val(int in);
void re_evaluate_kernels(void);
int volterra(double *y, double *t, double dt, int nt, int neq, int *istart, double *work);

#endif /* XPPAUT_VOLTERRA2_H */
