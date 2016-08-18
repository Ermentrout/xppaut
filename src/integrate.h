#ifndef XPPAUT_INTEGRATE_H
#define XPPAUT_INTEGRATE_H

#include <stdio.h>
#include "xpplim.h"

/* --- Macros --- */
#define MAXFP 400


/* --- Types --- */
typedef struct {
	char item[30];
	int steps,shoot,col,movie,mc;
	double plow,phigh;
} EQ_RANGE;


typedef struct {
	char item[30],item2[30];
	int steps,steps2,reset,oldic,index,index2,cycle,type,type2,movie;
	double plow,phigh,plow2,phigh2;
	int rtype;
} RANGE;


typedef struct {
	int nvec,node;
	double *x;
} XPPVEC;


/* --- Data --- */
extern double LastTime;
extern double MyData[MAXODE];
extern double MyTime;

extern int DelayErr;
extern int MakePlotFlag;
extern int MyStart;
extern int RANGE_FLAG;
extern int StableManifoldColor;
extern int UnstableManifoldColor;

extern int SuppressBounds;
extern int SuppressOut;

extern EQ_RANGE eq_range;
extern RANGE range;
extern XPPVEC xpv;

/* --- Functions --- */
void arr_ic_start(void);
void batch_integrate(void);
void comp_color(float *v1, float *v2, int n, double dt);
void cont_integ(void);
void do_init_data(int com);
int do_auto_range_go(void);
int do_range(double *x, int flag);
void export_data(FILE *fp);
int extract_ic_data(char *big);
void find_equilib_com(int com);
void get_ic(int it, double *x);
void init_ar_ic(void);
int integrate(double *t, double *x, double tend, double dt, int count, int nout, int *start);
void integrate_init_range(void);
void integrate_setup_range(void);
int ode_int(double *y, double *t, int *istart);
int one_step_int(double *y, double t0, double t1, int *istart);
void restore(int i1, int i2);
void run_now(void);
void send_halt(double *y, double t);
void send_output(double *y, double t);
void set_cycle(int flag, int *icol);
void shoot(double *x, double *xg, double *evec, int sgn);
void stop_integration(void);
void swap_color(int *col, int rorw);
void usual_integrate_stuff(double *x);

#endif /* XPPAUT_INTEGRATE_H */
