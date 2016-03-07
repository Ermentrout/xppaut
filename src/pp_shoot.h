#ifndef XPPAUT_PP_SHOOT_H
#define XPPAUT_PP_SHOOT_H

#include <stdio.h>

/* --- Types --- */
typedef struct {
	char item[30];
	int steps,side,cycle,movie;
	double plow,phigh;
} SHOOT_RANGE;


/* --- Data --- */
extern SHOOT_RANGE shoot_range;

/* --- Functions -- */
void compile_bvp(void);
void do_bc(double *y__0, double t0, double *y__1, double t1, double *f, int n);
void reset_bvp(void);
void init_shoot_range(char *s);
void find_bvp_com(int com);

#endif /* XPPAUT_PP_SHOOT_H */
