#ifndef _my_rhs_h_
#define _my_rhs_h_


/* my_rhs.c */
int MAIN__(void);
int main(int argc, char **argv);
void extra(double *y__y, double t, int nod, int neq);
void set_fix_rhs(double t, double *y);
int my_rhs(double t, double *y, double *ydot, int neq);
void update_based_on_current(void);
void fix_only(void);
void rhs_only(double *y, double *ydot);
void vec_rhs(double t, double *y, double *ydot, int neq);


#endif
