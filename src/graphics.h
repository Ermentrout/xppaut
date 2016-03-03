#ifndef XPPAUT_GRAPHICS_H
#define XPPAUT_GRAPHICS_H

#include <X11/Xlib.h>

/* --- Macros --- */
/* TextJustify values */
#define TJ_LEFT 0
#define TJ_CENTER 1
#define TJ_RIGHT 2

/* --- Data --- */
extern double XMin, YMin, XMax, YMax;
extern double PHI0;
extern double THETA0;

extern int avromfonts[5];
extern int avsymfonts[5];
extern int DBottom, DLeft, DRight, DTop;
extern int HChar, VChar;
extern int HTic, VTic;
extern int PointRadius;
extern int PointType;
extern int TextAngle;
extern int TextJustify;

extern XFontStruct *romfonts[5];
extern XFontStruct *symfonts[5];

/* --- Functions --- */
void bead_abs(double x1, double y1);
void change_current_linestyle(int new_int, int *old);
int clip(double x1, double x2, double y1, double y2, float *x1_out, float *y1_out, float *x2_out, float *y2_out);
void eq_symb(double *x, int type);
void fancy_text_abs(double x, double y, char *old, int size, int font);
void fillintext(char *old, char *new_char);
void frect_abs(double x1, double y1, double w, double h);
void get_draw_area(void);
void get_draw_area_flag(int flag);
void get_graph(void);
void get_scale(double *x1, double *y1, double *x2, double *y2);
void init_all_graph(void);
void init_ps(void);
void init_svg(void);
void init_x11(void);
void line(int x1, int y1, int x2, int y2);
void copy_graph(int i, int l);
void line_3d(double x, double y, double z, double xp, double yp, double zp);
void line_abs(double x1, double y1, double x2, double y2);
void line3d(double x01, double y01, double z01, double x02, double y02, double z02);
void line3dn(double xs1, double ys1, double zs1, double xsp1, double ysp1, double zsp1);
void make_rot(double theta, double phi);
void point_3d(double x, double y, double z);
void point_abs(double x1, double y1);
void put_text(int x, int y, char *str);
void put_text_x11(int x, int y, char *str);
void reset_all_line_type();
void reset_graph(void);
void scale3d(double x, double y, double z, float *xp, float *yp, float *zp);
void scale_dxdy(double x, double y, double *i, double *j);
void scale_to_screen(double x, double y, int *i, int *j);
void scale_to_real(int i, int j, float *x, float *y);
void set_extra_graphs(void);
void set_linestyle(int ls);
void set_normal_scale(void);
void set_scale(double x1, double y1, double x2, double y2);
void special_put_text_x11(int x, int y, char *str, int size);
void text3d(double x, double y, double z, char *s);
int threed_proj(double x, double y, double z, float *xp, float *yp);
int threedproj(double x2p, double y2p, double z2p, float *xp, float *yp);

#endif /* XPPAUT_GRAPHICS_H */
