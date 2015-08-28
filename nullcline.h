#ifndef _nullcline_h_
#define _nullcline_h_


#include <stdio.h>

typedef struct {
		float x,y,z;
		} Pt;

typedef struct nclines {
                float *xn,*yn;
                int nmx,nmy;
                int n_ix,n_iy;
                struct nclines *n,*p;
}  NCLINES;


void create_new_cline();
void froz_cline_stuff_com(int i);
void do_range_clines(void);
void start_ncline(void);
void clear_froz_cline(void);
int get_nullcline_floats(float **v, int *n, int who, int type);
void save_frozen_clines(char *fn);
void redraw_froz_cline(int flag);
void add_froz_cline(float *xn, int nmx, int n_ix, float *yn, int nmy, int n_iy);
void get_max_dfield(double *y, double *ydot, double u0, double v0, double du, double dv, int n, int inx, int iny, double *mdf);
void redraw_dfield(void);
void direct_field_com(int c);
void save_the_nullclines(void);
void restore_nullclines(void);
void dump_clines(FILE *fp, float *x, int nx, float *y, int ny);
void dump_clines_old(FILE *fp, float *x, int nx, float *y, int ny);
void restor_null(float *v, int n, int d);
void new_clines_com(int c);
void new_nullcline(int course, double xlo, double ylo, double xhi, double yhi, float *stor, int *npts);
void stor_null(double x1, double y1, double x2, double y2);
float fnull(double x, double y);
int interpolate(Pt p1, Pt p2, double z, float *x, float *y);
void quad_contour(Pt p1, Pt p2, Pt p3, Pt p4);
void triangle_contour(Pt p1, Pt p2, Pt p3);
void do_cline(int ngrid, double x1, double y1, double x2, double y2);
void do_batch_nclines();
void do_batch_dfield();

#endif
