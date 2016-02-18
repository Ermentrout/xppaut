#ifndef XPPAUT_AUTO_NOX_H
#define XPPAUT_AUTO_NOX_H

#include <stdio.h>
#include "autlim.h"
#include "xpplim.h"

/* --- Macros --- */
#define MAX_AUT_PER 10

/* --- Types --- */
typedef struct {
  int irot;
  int nrot[1000];
  double torper;
} ROTCHK;

typedef struct  {
  int exist;
  int ntst,nmx,npr;
  double ds,dsmax,dsmin,rl0,rl1,a0,a1;
  double xmin,xmax,ymin,ymax;
  double lastx,lasty;
  int wid,hgt,x0,y0,st_wid;
  int nfpar,nbc;
  int ips,irs,ilp,isp,isw,itp;
  int plot,var;
  int icp1,icp2,icp3,icp4,icp5;
  int nper;
  char hinttxt[256];
  double period[MAX_AUT_PER];
  int uzrpar[MAX_AUT_PER];
  double epsl,epsu,epss;
  int ncol;
}BIFUR;

typedef struct {
  int iad;
  int mxbf;
  int iid;
  int itmx;
  int itnw;
  int nwtn;
  int iads;

} ADVAUTO;

typedef struct {
  int package;
  int ibr,ntot,itp,lab;
  double norm,uhi[NAUTO],ulo[NAUTO],u0[NAUTO],ubar[NAUTO];
  double par[20],per,torper;
  int index,nfpar,icp1,icp2,icp3,icp4,icp5;
  int flag;
} GRABPT;


typedef struct diagram {
  int package;
  int ibr,ntot,itp,lab,calc;
  double norm,*uhi,*ulo,*u0,*ubar,*evr,*evi;
  double par[20],per,torper;
  int index,nfpar;
  int icp1,icp2,icp3,icp4,icp5,flag2;
  struct diagram *prev;
  struct diagram *next;
} DIAGRAM;

typedef struct {
  int plot,var,icp1,icp2,icp3,icp4,icp5;
  double xmin,ymin,xmax,ymax;
}  AUTOAX;

/* --- Data --- */
extern char fort3[200];
extern char fort7[200];
extern char fort8[200];
extern char fort9[200];

extern double auto_a0;
extern double auto_a1;
extern double auto_ds;
extern double auto_dsmax;
extern double auto_dsmin;
extern double auto_epsl;
extern double auto_epsu;
extern double auto_epss;
extern double auto_rl0;
extern double auto_rl1;
extern double auto_xmax;
extern double auto_xmin;
extern double auto_ymax;
extern double auto_ymin;
extern double homo_l[100];
extern double homo_r[100];
extern double outperiod[20];

extern int auto_ncol;
extern int auto_nmx;
extern int auto_npr;
extern int auto_ntst;
extern int auto_var;
extern int AutoPar[8];
extern int Auto_index_to_array[8];
extern int AutoTwoParam;
extern int HomoFlag;
extern int load_all_labeled_orbits;
extern int NAutoPar;
extern int NAutoUzr;
extern int NewPeriodFlag;
extern int SEc;
extern int SPc;
extern int UEc;
extern int UPc;
extern int sparity;
extern int TypeOfCalc;

extern unsigned int DONT_XORCross;

extern ADVAUTO aauto;
extern BIFUR Auto;
extern GRABPT grabpt;
extern ROTCHK blrtn;

/* --- Functions --- */
void colset(int type);
void pscolset2(int flag2);
void colset2(int flag2);
void get_auto_str(char *xlabel, char *ylabel);
void draw_ps_axes(void);
void draw_svg_axes(void);
void draw_bif_axes(void);
int IXVal(double x);
int IYVal(double y);
int chk_auto_bnds(int ix, int iy);
void renamef(char *old, char *new);
void copyf(char *old, char *new);
void appendf(char *old, char *new);
void deletef(char *old);
void close_auto(int flag);
void open_auto(int flag);
void do_auto(int iold, int isave, int itp);
void set_auto(void);
int auto_name_to_index(char *s);
int auto_par_to_name(int index, char *s);
void auto_per_par(void);
void auto_params(void);
void auto_num_par(void);
void auto_plot_par(void);
void auto_fit(void);
void auto_default(void);
void auto_zoom_in(int i1, int j1, int i2, int j2);
void auto_zoom_out(int i1, int j1, int i2, int j2);
void auto_xy_plot(double *x, double *y1, double *y2, double par1, double par2, double per, double *uhigh, double *ulow, double *ubar, double a);
int plot_point(int flag2, int icp1, int icp2);
void add_ps_point(double *par, double per, double *uhigh, double *ulow, double *ubar, double a, int type, int flag, int lab, int npar, int icp1, int icp2, int flag2, double *evr, double *evi);
void auto_line(double x1i, double y1i, double x2i, double y2i);
void add_point(double *par, double per, double *uhigh, double *ulow, double *ubar, double a, int type, int flag, int lab, int npar, int icp1, int icp2, int flag2,int icp3, int icp4, double *evr, double *evi);
void get_bif_sym(char *at, int itp);
void info_header(int flag2, int icp1, int icp2);
void new_info(int ibr, int pt, char *ty, int lab, double *par, double norm, double u0, double per, int flag2, int icp1, int icp2);
void traverse_out(DIAGRAM *d, int *ix, int *iy, int dodraw);
void do_auto_win(void);
void load_last_plot(int flag);
void keep_last_plot(int flag);
void init_auto_win(void);
void plot_stab(double *evr, double *evi, int n);
int yes_reset_auto(void);
int reset_auto(void);
void auto_grab(void);

void auto_start_diff_ss(void);
void auto_start_at_bvp(void);
void auto_start_at_per(void);
void find_best_homo_shift(int n);
void get_start_period(double *p);
void get_start_orbit(double *u, double t, double p, int n);
void get_shifted_orbit(double *u, double t, double p, int n);
void auto_new_ss(void);
void auto_new_discrete(void);
void auto_extend_ss(void);
void auto_start_choice(void);
void torus_choice(void);
void per_doub_choice(void);
void periodic_choice(void);
void hopf_choice(void);
void auto_new_per(void);
void auto_start_at_homoclinic(void);
int get_homo_info(int flg,int *nun,int *nst,double *ul, double *ur);
void auto_extend_homoclinic(void);
void auto_extend_bvp(void);
void auto_switch_per(void);
void auto_switch_bvp(void);
void auto_switch_ss(void);
void auto_2p_limit(int ips);
void auto_twopar_double(void);
void auto_torus(void);
void auto_2p_branch(int ips);
void auto_branch_choice(int ibr, int ips);
void auto_homo_choice(int itp);
void auto_2p_fixper(void);
void auto_2p_hopf(void);
void auto_period_double(void);
void auto_err(char *s);
void auto_run(void);
void load_auto_orbit(void);
void load_auto_orbitx(int ibr,int flag, int lab, double per);
void save_auto(void);
void save_auto_numerics(FILE *fp);
void load_auto_numerics(FILE *fp);
void save_auto_graph(FILE *fp);
void load_auto_graph(FILE *fp);
void save_q_file(FILE *fp);
void make_q_file(FILE *fp);
int noinfo(char *s);
void load_auto(void);
int move_to_label(int mylab, int *nrow, int *ndim, FILE *fp);
void get_a_row(double *u, double *t, int n, FILE *fp);
void auto_file(void);
int check_plot_type(int flag2,int icp1, int icp2);

#endif /* XPPAUT_AUTO_NOX_H */
