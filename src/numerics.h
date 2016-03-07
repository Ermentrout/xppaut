#ifndef XPPAUT_NUMERICS_H
#define XPPAUT_NUMERICS_H

/* --- Macros --- */
#define VOLTERRA 6
#define BACKEUL 7
#define RKQS 8
#define STIFF 9
#define CVODE 10
#define GEAR 5
#define DP5 11
#define DP83 12
#define RB23 13
#define SYMPLECT 14


/* --- Data --- */
extern int cv_bandflag;
extern int cv_bandlower;
extern int cv_bandupper;
extern int METHOD;

/* --- Functions --- */
void check_delay(void);
void do_meth(void);
void get_num_par(int ch);
void set_delay(void);
void set_total(double total);

#endif /* XPPAUT_NUMERICS_H */
