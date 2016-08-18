#ifndef XPPAUT_AUTO_X11_H
#define XPPAUT_AUTO_X11_H

#include <X11/Xlib.h>

/* --- Data --- */
extern int AutoRedrawFlag;
extern int mark_flag;
extern int mark_ibrs;
extern int mark_ipte;
extern int mark_ipts;

/* --- Functions --- */
void ALINE(int a, int b, int c, int d);
void ATEXT(int a, int b, char *c);
void auto_button(XEvent ev);
void auto_enter(Window w, int v);
void auto_get_info( int *n, char *pname);
void auto_keypress(XEvent ev, int *used);
void auto_motion(XEvent ev);
int auto_pop_up_list(char *title, char **list, char *key, int n, int max, int def, int x, int y, char **hints, char *httxt);
int auto_rubber(int *i1, int *j1, int *i2, int *j2, int flag);
void auto_scroll_window(void);
void auto_set_mark(int i);
void auto_stab_line(int x, int y, int xp, int yp);
void autobw(void);
void autocol(int col);
int byeauto_(int *iflag);
void Circle(int x, int y, int r);
void clear_auto_info(void);
void clear_auto_plot(void);
void clr_stab(void);
void display_auto(Window w);
void DLINE(double a, double b, double c, double d);
void do_auto_range(void);
void draw_auto_info(char *bob, int x, int y);
void FillCircle(int x, int y, int r);
void LineWidth(int wid);
void make_auto(char *wname, char *iname);
void redraw_auto_menus(void);
void refreshdisplay(void);
void resize_auto_window(XEvent ev);
void traverse_diagram(void);

#endif /* XPPAUT_AUTO_X11_H */
