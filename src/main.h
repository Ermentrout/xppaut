#ifndef XPPAUT_MAIN_H
#define XPPAUT_MAIN_H

#include <X11/Xlib.h>

/* --- Data --- */
extern Display *display;
extern GC gc;
extern GC font_gc;
extern GC gc_graph;
extern GC small_gc;
extern Window command_pop;
extern Window info_pop;

/* --- Functions --- */
void do_main(int argc, char **argv);
void check_for_quiet(int argc, char **argv);
void do_vis_env(void);
void init_X(void);
void set_big_font(void);
void set_small_font(void);
void xpp_events(XEvent report, int min_wid, int min_hgt);
void do_events(unsigned int min_wid, unsigned int min_hgt);
void bye_bye(void);
void clr_scrn(void);
void redraw_all(void);
void commander(int ch);
Window init_win(unsigned int bw, char *icon_name, char *win_name, int x, int y, unsigned int min_wid, unsigned int min_hgt, int argc, char **argv);
void top_button_draw(Window w);
void top_button_cross(Window w, int b);
void top_button_press(Window w);
void top_button_events(XEvent report);
void make_top_buttons(void);
void getGC(GC *gc);
void load_fonts(void);
void make_pops(void);
void FixWindowSize(Window w, int width, int height, int flag);
int getxcolors(XWindowAttributes *win_info, XColor **colors);
void test_color_info(void);

#endif /* XPPAUT_MAIN_H */
