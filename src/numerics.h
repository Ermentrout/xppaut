#ifndef XPPAUT_NUMERICS_H
#define XPPAUT_NUMERICS_H

/* --- Types --- */
typedef enum {
	METHOD_DISCRETE = 0,
	METHOD_EULER,
	METHOD_MODEULER,
	METHOD_RK4,
	METHOD_ADAMS,
	METHOD_GEAR,
	METHOD_VOLTERRA,
	METHOD_BACKEUL,
	METHOD_RKQS,
	METHOD_STIFF,
	METHOD_CVODE,
	METHOD_DP5,
	METHOD_DP83,
	METHOD_RB23,
	METHOD_SYMPLECT,
	NUM_METHODS
} Method;


/* --- Data --- */
extern int cv_bandflag;
extern int cv_bandlower;
extern int cv_bandupper;
extern Method METHOD;

/* --- Functions --- */
void check_delay(void);
void do_meth(void);
void get_num_par(int ch);
void set_delay(void);
void set_total(double total);

#endif /* XPPAUT_NUMERICS_H */
