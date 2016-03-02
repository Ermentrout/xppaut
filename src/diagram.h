#ifndef XPPAUT_DIAGRAM_H
#define XPPAUT_DIAGRAM_H

#include <stdio.h>

/* --- Types --- */
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

/* --- Data --- */
extern DIAGRAM *bifd;
extern int NBifs;

/* --- Functions --- */
void add_diagram(int ibr, int ntot, int itp, int lab, int nfpar, double a, double *uhi, double *ulo, double *u0, double *ubar, double *par, double per, int n, int icp1, int icp2, int icp3,int icp4,int flag2, double *evr, double *evi);
void bound_diagram(double *xlo, double *xhi, double *ylo, double *yhi);
void edit_start(int ibr, int ntot, int itp, int lab, int nfpar, double a, double *uhi, double *ulo, double *u0, double *ubar, double *par, double per, int n, int icp1, int icp2, int icp3, int icp4,double *evr, double *evi);
void kill_diagrams(void);
void load_browser_with_branch(int ibr,int pts,int pte);
int load_diagram(FILE *fp, int node);
void post_auto(void);
void redraw_diagram(void);
int save_diagram(FILE *fp, int n);
void start_diagram(int n);
void svg_auto(void);
void write_init_data_file(void);
void write_info_out(void);
void write_pts(void);

#endif /* XPPAUT_DIAGRAM_H */
