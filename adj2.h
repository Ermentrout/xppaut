#ifndef _adj2_h_
#define _adj2_h_

#include <stdio.h>

void init_trans(void);
void dump_transpose_info(FILE *fp, int f);
int do_transpose(void);
int create_transpose(void);
void alloc_h_stuff(void);
void data_back(void);
void adj_back(void);
void h_back(void);
void make_adj_com(int com);
void adjoint_parameters(void);
void new_h_fun(int silent);
void dump_h_stuff(FILE *fp, int f);
int make_h(float **orb, float **adj, float **h, int nt, double dt, int node,int silent);
void new_adjoint(void);
void test_test(void);
void compute_one_orbit(double *ic, double per);
int adjoint(float **orbit, float **adjnt, int nt, double dt, double eps, double minerr, int maxit, int node);
void eval_rhs(double **jac, int k1, int k2, double t, double *y, double *yp, int node);
int rk_interp(double **jac, int k1, int k2, double *y, double *work, int neq, double del, int nstep);
int step_eul(double **jac, int k, int k2, double *yold, double *work, int node, double dt);
void do_liapunov(void);
void alloc_liap(int n);
void do_this_liaprun(int i, double p);
void norm_vec(double *v, double *mu, int n);
int hrw_liapunov(double *liap, int batch, double eps);



#endif

