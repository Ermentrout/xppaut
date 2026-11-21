
#ifndef _delay_handle_h_
#define _delay_handle_h_


/* delay_handle.c */
double delay_stab_eval(double delay, int var);
int alloc_delay(double big);
void free_delay(void);
void stor_delay(double *y, double t);
int delay_enter_steady(void);
void delay_leave_steady(int old);
double get_delay_old(int in, double tau);
void polint(double *xa, double *ya, int n, double x, double *y, double *dy);
double get_delay(int in, double tau);
int do_init_delay(double big);


#endif
