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
int alloc_delay(double big);
double delay_stab_eval(double delay, int var);
int do_init_delay(double big);
void free_delay(void);
double get_delay(int in, double tau);
void stor_delay(double *y);

#endif /* XPPAUT_DELAY_HANDLE_H */
