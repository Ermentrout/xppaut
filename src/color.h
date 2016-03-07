#ifndef XPPAUT_COLOR_H
#define XPPAUT_COLOR_H

#include <X11/Xlib.h>

/* --- Macros --- */
#define FIRSTCOLOR 30

/* --- Data --- */
extern int COLOR;
extern int ColorizeFlag;
extern int color_total;
extern double ColorViaHi;
extern double ColorViaLo;
extern char ColorVia[15];
extern int custom_color;
extern int periodic;

extern char *color_names[12];
extern int colorline[12];

/* --- Functions --- */
int ColorMap(int i);
void get_ps_color(int i, float *r, float *g, float *b);
void get_svg_color(int i,int *r,int *g,int *b);
void MakeColormap(void);
void NewColormap(int type);
void set_color(int col);
void set_colorization_stuff(void);
void set_scolor(int col);
void user_set_color_par(int flag,char *via,double lo,double hi);

#endif /* XPPAUT_COLOR_H */
