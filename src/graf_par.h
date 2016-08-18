#ifndef XPPAUT_GRAF_PAR_H
#define XPPAUT_GRAF_PAR_H

#include <stdio.h>
#include <X11/Xlib.h>

/* --- Macros --- */
#define SCRNFMT 0
#define PSFMT 1
#define SVGFMT 2

/* type values */
#define REAL_SMALL 1.e-6
#define MAXBIFCRV 100
#define lmax(a,b) ((a>b) ? a : b)


/* --- Data --- */
extern int AutoFreezeFlag;

/* --- Functions --- */
void add_a_curve_com(int c);
void auto_freeze_it(void);
void change_cmap_com(int i);
void change_view_com(int com);
void check_windows(void);
void default_window(void);
void draw_freeze(Window w);
void freeze_com(int c);
void get_3d_par_com(void);
void get_max(int index, double *vmin, double *vmax);
void graph_all(int *list, int n, int type);
void ind_to_sym(int ind, char *str);
void init_bd(void);
void key_frz_com(int c);
void redraw_the_graph(void);
void window_zoom_com(int c);
void xi_vs_t(void);

#endif /* XPPAUT_GRAF_PAR_H */
