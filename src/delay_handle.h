#ifndef XPPAUT_DELAY_HANDLE_H
#define XPPAUT_DELAY_HANDLE_H

#include "xpplim.h"

/* --- Data --- */
extern double AlphaMax;
extern int del_stab_flag;
extern int DelayFlag;
extern int DelayGrid;
extern double delay_list[MAXDELAY];
extern int NDelay;
extern double OmegaMax;
extern double variable_shift[2][MAXODE];
extern int WhichDelay;

/* --- Functions --- */
double delay_stab_eval(double delay, int var);
int alloc_delay(double big);
void free_delay(void);
void stor_delay(double *y);
double get_delay_old(int in, double tau);
void polint(double *xa, double *ya, int n, double x, double *y, double *dy);
double get_delay(int in, double tau);
int do_init_delay(double big);

#endif /* XPPAUT_DELAY_HANDLE_H */
