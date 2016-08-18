#ifndef XPPAUT_TXTREAD_H
#define XPPAUT_TXTREAD_H

#include <X11/Xlib.h>

void init_txtview(void);
void make_txtview(void);
void redraw_txtview(Window w);
void txt_view_events(XEvent ev);

#endif /* XPPAUT_TXTREAD_H */
