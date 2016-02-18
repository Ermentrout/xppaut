#ifndef XPPAUT_NUMERICS_H
#define XPPAUT_NUMERICS_H

/* --- Types --- */
typedef struct {
  double tmod;
  int maxvar,sos,type,sign;
  char section[256];
  int formula[256];
} POINCARE_MAP;

/* --- Data --- */
extern int cv_bandflag;
extern int cv_bandlower;
extern int cv_bandupper;

/* --- Functions --- */
void chk_delay(void);
void chk_volterra(void);
void check_pos(int *j);
void compute_one_period(double period,double *x, char *name);
void do_meth(void);
void get_method(void);
void get_num_par(int ch);
void get_pmap_pars_com(int l);
void init_numerics(void);
void meth_dialog(void);
void quick_num(int com);
void ruelle(void);
void set_col_par_com(int i);
void set_delay(void);
void set_total(double total);
void user_set_color_par(int flag,char *via,double lo,double hi);

#endif /* XPPAUT_NUMERICS_H */
