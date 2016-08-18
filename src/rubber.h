#ifndef XPPAUT_RUBBER_H
#define XPPAUT_RUBBER_H

#include <X11/Xlib.h>

/* --- Macros --- */
/* f values */
#define RUBBOX 0
#define RUBLINE 1

/* --- Functions --- */
int rubber(int *x1, int *y1, int *x2, int *y2, Window w, int f);

#endif /* XPPAUT_RUBBER_H */
