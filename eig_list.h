#ifndef _eig_list_h_
#define _eig_list_h_

#include <X11/Xlib.h>

void draw_eq_list(Window w);
void create_eq_list(void);
void eq_list_keypress(XEvent ev, int *used);
void enter_eq_stuff(Window w, int b);
void eq_list_button(XEvent ev);
void eq_list_up(void);
void eq_list_down(void);
void eq_box_import(void);
void get_new_size(Window win, unsigned int *wid, unsigned int *hgt);
void resize_eq_list(Window win);
void eq_box_button(Window w);
void create_eq_box(int cp, int cm, int rp, int rm, int im, double *y, double *ev, int n);
void draw_eq_box(Window w);


#endif
