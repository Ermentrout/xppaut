#ifndef XPPAUT_MY_RHS_H
#define XPPAUT_MY_RHS_H

/* --- Functions --- */
void extra(double *y__y, double t, int nod, int neq);
void fix_only(void);
int my_rhs(double t, double *y, double *ydot, int neq);
void rhs_only(double *y, double *ydot);
void set_fix_rhs(double t, double *y);
void update_based_on_current(void);

#endif /* XPPAUT_MY_RHS_H */
