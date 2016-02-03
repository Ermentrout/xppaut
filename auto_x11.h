#ifndef _auto_x11_h_
#define _auto_x11_h_

#include <X11/Xlib.h>

typedef struct  {
  Window canvas, axes,numerics,grab,next,run,clear,redraw,base,per;
  Window info,param,file,abort,stab,hint,kill;
} AUTOWIN;


void ALINE(int a, int b, int c, int d);
void DLINE(double a, double b, double c, double d);
void ATEXT(int a, int b, char *c);
void clr_stab(void);
void auto_stab_line(int x, int y, int xp, int yp);
void clear_auto_plot(void);
void redraw_auto_menus(void);
void traverse_diagram(void);
void clear_auto_info(void);
void draw_auto_info(char *bob, int x, int y);
void refreshdisplay(void);
int byeauto_(int *iflag);
void Circle(int x, int y, int r);
void autocol(int col);
void autobw(void);
int auto_rubber(int *i1, int *j1, int *i2, int *j2, int flag);
int auto_pop_up_list(char *title, char **list, char *key, int n, int max, int def, int x, int y, char **hints, char *httxt);
void MarkAuto(int x,int y);
void XORCross(int x, int y);
void FillCircle(int x, int y, int r);
void LineWidth(int wid);
void auto_motion(XEvent ev);
void display_auto(Window w);
Window lil_button(Window root, int x, int y, char *name);
void make_auto(char *wname, char *iname);
void resize_auto_window(XEvent ev);
void a_msg(int i, int v);
void auto_enter(Window w, int v);
void auto_button(XEvent ev);
void auto_kill(void);
void auto_keypress(XEvent ev, int *used);
int query_special(char* title,char *nsymb);
void clear_msg();
void find_point(int ibr,int pt);
void auto_get_info( int *n, char *pname);
void auto_set_mark(int i);
void do_auto_range();
void RedrawMark();

#endif
