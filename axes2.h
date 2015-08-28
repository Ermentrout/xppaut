
#ifndef _axes2_h_
#define _axes2_h_

#include <X11/Xlib.h>

void re_title(void);
void get_title_str(char *s1, char *s2, char *s3);
void make_title(char *str);
double dbl_raise(double x, int y);
double make_tics(double tmin, double tmax);
void find_max_min_tic(double *tmin, double *tmax, double tic);
void redraw_cube_pt(double theta, double phi);
void do_axes(void);
void redraw_cube(double theta, double phi);
void draw_unit_cube(void);
void Frame_3d(void);
void Box_axis(double x_min, double x_max, double y_min, double y_max, char *sx, char *sy, int flag);
void draw_ytics(char *s1, double start, double incr, double end);
void draw_xtics(char *s2, double start, double incr, double end);

#endif
