#ifndef XPPAUT_EIG_LIST_H
#define XPPAUT_EIG_LIST_H

#include <X11/Xlib.h>

/* --- Functions --- */
void create_eq_box(int cp, int cm, int rp, int rm, int im, double *y, double *ev, int n);
void create_eq_list(void);
void draw_eq_box(Window w);
void draw_eq_list(Window w);
void enter_eq_stuff(Window w, int b);
void eq_list_button(XEvent ev);
void eq_list_keypress(XEvent ev, int *used);
void get_new_size(Window win, unsigned int *wid, unsigned int *hgt);
void resize_eq_list(Window win);

#endif /* XPPAUT_EIG_LIST_H */
