#ifndef XPPAUT_GEAR_H
#define XPPAUT_GEAR_H

#include "../xpplim.h"

/* --- Functions --- */
int gear(int n, double *t, double tout, double *y, double hmin, double hmax, double eps, int mf, double *error, int *kflag, int *jstart, double *work);
void gear_err_msg(int kflag);

#endif /* XPPAUT_GEAR_H */
