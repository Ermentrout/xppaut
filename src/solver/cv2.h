#ifndef XPPAUT_SOLVER_CV2_H
#define XPPAUT_SOLVER_CV2_H

void cvode_err_msg(int kflag);
int cvode(int *command, double *y, double *t, int n, double tout, int *kflag, double *atol, double *rtol);
void end_cvode(void);

#endif // XPPAUT_SOLVER_CV2_H

