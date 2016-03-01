#ifndef XPPAUT_NULLCLINE_H
#define XPPAUT_NULLCLINE_H

#include <stdio.h>

/* --- Types --- */
typedef struct {
	float x,y,z;
} Pt;

typedef struct nclines {
	float *xn,*yn;
	int nmx,nmy;
	int n_ix,n_iy;
	struct nclines *n,*p;
}  NCLINES;

/* --- Data --- */
extern double ColorViaHi;
extern double ColorViaLo;
extern int ColorizeFlag;
extern char ColorVia[15];
extern int DF_FLAG;
extern int DF_GRID;
extern int DFBatch;
extern int DOING_DFIELD;
extern int NCBatch;
extern int XNullColor;
extern int YNullColor;

/* --- Functions --- */
void add_froz_cline(float *xn, int nmx, int n_ix, float *yn, int nmy, int n_iy);
void clear_froz_cline(void);
void create_new_cline();
void direct_field_com(int c);
void do_batch_dfield();
void do_batch_nclines();
void do_cline(int ngrid, double x1, double y1, double x2, double y2);
void do_range_clines(void);
void dump_clines(FILE *fp, float *x, int nx, float *y, int ny);
void dump_clines_old(FILE *fp, float *x, int nx, float *y, int ny);
float fnull(double x, double y);
void froz_cline_stuff_com(int i);
void get_max_dfield(double *y, double *ydot, double u0, double v0, double du, double dv, int n, int inx, int iny, double *mdf);
int get_nullcline_floats(float **v, int *n, int who, int type);
int interpolate(Pt p1, Pt p2, double z, float *x, float *y);
void new_clines_com(int c);
void new_nullcline(int course, double xlo, double ylo, double xhi, double yhi, float *stor, int *npts);
void quad_contour(Pt p1, Pt p2, Pt p3, Pt p4);
void redraw_dfield(void);
void redraw_froz_cline(int flag);
void restor_null(float *v, int n, int d);
void restore_nullclines(void);
void save_frozen_clines(char *fn);
void save_the_nullclines(void);
void silent_dfields(void);
void silent_nullclines(void);
void set_colorization_stuff(void);
void start_ncline(void);
void stor_null(double x1, double y1, double x2, double y2);
void triangle_contour(Pt p1, Pt p2, Pt p3);

#endif /* XPPAUT_NULLCLINE_H */
