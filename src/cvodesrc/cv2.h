#ifndef _cv2_h_
#define _cv2_h_


/* cv2.c */
void start_cv(double *y, double t, int n, double tout, double *atol, double *rtol);
void end_cv(void);
void cvode_err_msg(int kflag);
int cvode(int *command, double *y, double *t, int n, double tout, int *kflag, double *atol, double *rtol);
int ccvode(int *command, double *y, double *t, int n, double tout, int *kflag, double *atol, double *rtol);


#endif

