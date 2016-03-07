#ifndef XPPAUT_MANY_POPS_H
#define XPPAUT_MANY_POPS_H

#include <X11/Xlib.h>
#include "struct.h"

/* --- Data --- */
extern int ActiveWinList[MAXPOP];
extern int current_pop;
extern int num_pops;
extern int SimulPlotFlag;

extern CURVE frz[MAXFRZ];
extern GRAPH graph[MAXPOP];
extern GRAPH *MyGraph;

/* --- Functions --- */
void add_label(char *s, int x, int y, int size, int font);
void BaseCol(void);
void canvas_xy(char *buf);
void change_plot_vars(int k);
int check_active_plot(int k);
void check_draw_button(XEvent ev);
void create_a_pop(void);
void do_expose(XEvent ev);
void do_gr_objs_com(int com);
void do_motion_events(XEvent ev);
void do_windows_com(int c);
void draw_label(Window w);
void edit_object_com(int com);
void get_intern_set(void);
void gtitle_text(char *string, Window win);
void GrCol(void);
void hi_lite(Window wi);
void init_grafs(int x, int y, int w, int h);
void make_active(int i,int flag);
void make_icon(char *icon, int wid, int hgt, Window w);
void resize_all_pops(int wid, int hgt);
void restore_off(void);
void restore_on(void);
int rotate3dcheck(XEvent ev);
int select_table(void);
void set_active_windows(void);
void SmallBase(void);
void SmallGr(void);
void svg_restore(void);
void title_text(char *string);

#endif /* XPPAUT_MANY_POPS_H */
