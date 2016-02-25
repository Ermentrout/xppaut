#ifndef XPPAUT_COLOR_H
#define XPPAUT_COLOR_H

#include <X11/Xlib.h>

/* --- Macros --- */
#define FIRSTCOLOR 30

/* --- Data --- */
extern int COLOR;
extern int color_total;
extern int custom_color;
extern int periodic;

/* --- Functions --- */
void tst_color(Window w);
void set_scolor(int col);
void set_color(int col);
void make_cmaps(int *r, int *g, int *b, int n, int type);
int rfun(double y, int per);
int gfun(double y, int per);
int bfun(double y, int per);
void NewColormap(int type);
void get_ps_color(int i, float *r, float *g, float *b);
void get_svg_color(int i,int *r,int *g,int *b);
void MakeColormap(void);
int ColorMap(int i);

#endif /* XPPAUT_COLOR_H */
