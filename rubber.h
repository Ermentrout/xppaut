
#ifndef _rubber_h
#define _rubber_h

#include <X11/Xlib.h>


int rubber(int *x1, int *y1, int *x2, int *y2, Window w, int f);
void rbox(int i1, int j1, int i2, int j2, Window w, int f);

#endif
