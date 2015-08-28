#ifndef _integrate_h_
#define _integrate_h_

#include <stdio.h>

void init_ar_ic(void);
void dump_range(FILE *fp, int f);
void init_range(void);
int set_up_eq_range(void);
void cont_integ(void);
int range_item(void);
int range_item2(void);
int set_up_range(void);
int set_up_range2(void);
void init_monte_carlo(void);
void monte_carlo(void);
void do_monte_carlo_search(int append, int stuffbrowse,int ishoot);
void do_eq_range(double *x);
void swap_color(int *col, int rorw);
void set_cycle(int flag, int *icol);
int do_range(double *x, int flag);
void find_equilib_com(int com);
void batch_integrate(void);
void do_batch_dry_run(void);
void batch_integrate_once(void);
int write_this_run(char *file, int i);
void do_init_data(int com);
void run_from_x(double *x);
void run_now(void);
void do_start_flags(double *x, double *t);
void usual_integrate_stuff(double *x);
void do_new_array_ic(char *new, int j1, int j2);
void store_new_array_ic(char *new, int j1, int j2, char *formula);
void evaluate_ar_ic(char *v, char *f, int j1, int j2);
int extract_ic_data(char *big);
void arr_ic_start(void);
int set_array_ic(void);
int form_ic(void);
void get_ic(int it, double *x);
int ode_int(double *y, double *t, int *istart, int ishow);
int integrate(double *t, double *x, double tend, double dt, int count, int nout, int *start);
void send_halt(double *y, double t);
void send_output(double *y, double t);
void do_plot(float *oldxpl, float *oldypl, float *oldzpl, float *xpl, float *ypl, float *zpl);
void export_data(FILE *fp);
void plot_the_graphs(float *xv, float *xvold, int node, int neq, double ddt, int *tc);
void plot_one_graph(float *xv, float *xvold, int node, int neq, double ddt, int *tc);
void restore(int i1, int i2);
void comp_color(float *v1, float *v2, int n, double dt);
void shoot(double *x, double *xg, double *evec, int sgn);
void shoot_easy(double *x);
void stop_integration(void);
int stor_full(void);
int do_auto_range_go();

#endif
