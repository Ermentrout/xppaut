#ifndef XPPAUT_AUTO_NOX_H
#define XPPAUT_AUTO_NOX_H

#include <stdio.h>
#include "autlim.h"
#include "xpplim.h"

/* --- Macros --- */
#define MAX_AUT_PER 10

enum BifurcationType {
	STABLE_EQUILIBRIUM = 1,
	UNSTABLE_EQUILIBRIUM,
	STABLE_PERIODIC,
	UNSTABLE_PERIODIC,
	SPECIAL
};

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
	char hinttxt[DEFAULT_STRING_LENGTH];
	double period[MAX_AUT_PER];
	int uzrpar[MAX_AUT_PER];
	double epsl,epsu,epss;
	int ncol;
} BIFUR;

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
int IXVal(double x);
int IYVal(double y);

void add_point(double *par, double per, double *uhigh, double *ulow, double *ubar, double a, int type, int flag, int lab, int npar, int icp1, int icp2,int icp3, int icp4, int flag2, double *evr, double *evi);
void add_ps_point(double *par, double per, double *uhigh, double *ulow, double *ubar, double a, int type, int flag, int lab, int npar, int icp1, int icp2, int flag2, double *evr, double *evi);
void auto_file(void);
void auto_grab(void);
void auto_num_par(void);
void auto_params(void);
void auto_per_par(void);
void auto_plot_par(void);
void auto_run(void);
void auto_xy_plot(double *x, double *y1, double *y2, double par1, double par2, double per, double *uhigh, double *ulow, double *ubar, double a);
int check_plot_type(int flag2,int icp1, int icp2);
void do_auto_win(void);
void draw_bif_axes(void);
void draw_ps_axes(void);
void draw_svg_axes(void);
void get_bif_sym(char *at, int itp);
void get_shifted_orbit(double *u, double t, double p, int n);
void get_start_orbit(double *u, double t, double p, int n);
void get_start_period(double *p);
void init_auto_win(void);
void load_auto_orbitx(int ibr,int flag, int lab, double per);
void new_info(int ibr, int pt, char *ty, int lab, double *par, double norm, double u0, double per, int flag2, int icp1, int icp2);
void plot_stab(double *evr, double *evi, int n);
void set_auto(void);
void storeautopoint(double x,double y);
int yes_reset_auto(void);

#endif /* XPPAUT_AUTO_NOX_H */
