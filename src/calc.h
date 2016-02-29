#ifndef XPPAUT_CALC_H
#define XPPAUT_CALC_H

#include <X11/Xlib.h>

/* --- Types --- */
typedef struct {
	Window base,quit,answer;
	double last_val;
	int use;
} MY_CALC;

/* --- Functions --- */
void draw_calc(Window w);
void make_calc(double z);
void quit_calc(void);
void ini_calc_string(char *name, char *value, int *pos, int *col);
void q_calc(void);
int do_calc(char *temp, double *z);
int has_eq(char *z, char *w, int *where);
double calculate(char *expr, int *ok);

#endif /* XPPAUT_CALC_H */
