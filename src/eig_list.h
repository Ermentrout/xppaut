#ifndef XPPAUT_EIG_LIST_H
#define XPPAUT_EIG_LIST_H

#include <X11/Xlib.h>
#include "xpplim.h"

/* --- Types --- */
typedef struct {
  Window base,stab,rest,top,close,import;
  double y[MAXODE],ev[MAXODE+MAXODE];
  int n,flag;
  int info[5];
  char type[15];
} EQ_BOX;

typedef struct{
  Window base,up,down,list,main,close;
  int istart,nlines,flag;
} EQ_LIST;

/* --- Functions --- */
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

#endif /* XPPAUT_EIG_LIST_H */
