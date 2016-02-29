#ifndef XPPAUT_PP_SHOOT_H
#define XPPAUT_PP_SHOOT_H

#include <stdio.h>

/* --- Types --- */
typedef struct {
	char item[30];
	int steps,side,cycle,movie;
	double plow,phigh;
} SHOOT_RANGE;

/* --- Functions -- */
void do_bc(double *y__0, double t0, double *y__1, double t1, double *f, int n);
void compile_bvp(void);
void reset_bvp(void);
void init_shoot_range(char *s);
void dump_shoot_range(FILE *fp, int f);
void bad_shoot(int iret);
void do_sh_range(double *ystart, double *yend);
int set_up_homoclinic(void);
int set_up_periodic(int *ipar, int *ivar, double *sect, int *ishow);
void find_bvp_com(int com);
void last_shot(int flag);
int set_up_sh_range(void);
void bvshoot(double *y, double *yend, double err, double eps, int maxit, int *iret, int n, int ishow, int iper, int ipar, int ivar, double sect);

#endif /* XPPAUT_PP_SHOOT_H */
