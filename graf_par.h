#ifndef _graf_par_h_
#define _graf_par_h_


#define RUBBOX 0
#define RUBLINE 1

#define SCRNFMT 0
#define PSFMT 1
#define SVGFMT 2


#define REAL_SMALL 1.e-6
#define MAX_LEN_SBOX 25
#define MAXBIFCRV 100
#define lmax(a,b) ((a>b) ? a : b)

#include <X11/Xlib.h>
#include <stdio.h>

typedef struct {
  char angle[20];
  char yes[3];
  double start;
  double incr;
  int nclip;
} MOV3D;

typedef struct {
  float *x[MAXBIFCRV],*y[MAXBIFCRV];
  int color[MAXBIFCRV],npts[MAXBIFCRV],nbifcrv;
  Window w;
} BD;



void change_view_com(int com);
void ind_to_sym(int ind, char *str);
void check_flags(void);
void get_2d_view(int ind);
void axes_opts(void);
void get_3d_view(int ind);
void check_val(double *x1, double *x2, double *xb, double *xd);
void get_max(int index, double *vmin, double *vmax);
void pretty(double *x1, double *x2);
void corner_cube(double *xlo, double *xhi, double *ylo, double *yhi);
void fit_window(void);
void check_windows(void);
void user_window(void);
void xi_vs_t(void);
void redraw_the_graph(void);
void movie_rot(double start, double increment, int nclip, int angle);
void test_rot(void);
void get_3d_par_com(void);
void get_3d_par_noper(void);
void window_zoom_com(int c);
void zoom_in(int i1, int j1, int i2, int j2);
void zoom_out(int i1, int j1, int i2, int j2);
void graph_all(int *list, int n, int type);
int find_color(int in);
int alter_curve(char *title, int in_it, int n);
void edit_curve(void);
void new_curve(void);
void create_ps(void);
void change_cmap_com(int i);
void freeze_com(int c);
void set_key(int x, int y);
void draw_freeze_key(void);
void key_frz_com(int c);
void edit_frz(void);
void delete_frz_crv(int i);
void delete_frz(void);
void kill_frz(void);
int freeze_crv(int ind);
void auto_freeze_it(void);
int create_crv(int ind);
void edit_frz_crv(int i);
void draw_frozen_cline(int index, Window w);
void draw_freeze(Window w);
void init_bd(void);
void draw_bd(Window w);
void free_bd(void);
void add_bd_crv(float *x, float *y, int len, int type, int ncrv);
void frz_bd(void);
void read_bd(FILE *fp);
int get_frz_index(Window w);
void export_graf_data(void);
void add_a_curve_com(int c);
void default_window();
void dump_ps( int i);

#endif
