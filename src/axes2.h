#ifndef XPPAUT_AXES2_H
#define XPPAUT_AXES2_H

#include <X11/Xlib.h>

/* --- Data --- */
extern int DOING_AXES;
extern int DOING_BOX_AXES;

/* --- Functions --- */
void Box_axis(double x_min, double x_max, double y_min, double y_max, char *sx, char *sy, int flag);
void do_axes(void);
void redraw_cube(double theta, double phi);
void redraw_cube_pt(double theta, double phi);

#endif /* XPPAUT_AXES2_H */
