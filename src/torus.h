#ifndef XPPAUT_TORUS_H
#define XPPAUT_TORUS_H

#include <X11/Xlib.h>

#include "xpplim.h"

/* --- Types --- */
typedef struct {
    Window base,done,cancel;
    Window w[MAXODE];
} Torbox;

void do_torus_com(int c);
void draw_tor_var(int i);
void draw_torus_box(Window win);
void choose_torus(void);
void make_tor_box(char *title);
void do_torus_events(void);

#endif /* XPPAUT_TORUS_H */
