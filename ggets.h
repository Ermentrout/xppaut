
#ifndef _ggets_h
#define _ggets_h

#include <X11/Xlib.h>

#define MaxIncludeFiles 10
#define ClickTime 200

void ping(void);
void reset_graphics(void);
void blank_screen(Window w);
void set_fore(void);
void set_back(void);
void showchar(int ch, int col, int row, Window or);
void chk_xor(void);
void set_gcurs(int y, int x);
void clr_command(void);
void draw_info_pop(Window win);
void bottom_msg(int line, char *msg);
void gputs(char *string, Window win);
void err_msg(char *string);
int plintf(char *fmt, ...);
int show_position(XEvent ev, int *com);
void gpos_prn(char *string, int row, int col);
void put_command(char *string);
int get_key_press(XEvent *ev);
void cput_text(void);
int get_mouse_xy(int *x, int *y, Window w);
void Ftext(int x, int y, char *string, Window o);
void bar(int x, int y, int x2, int y2, Window w);
void rectangle(int x, int y, int x2, int y2, Window w);
void setfillstyle(int type, int color);
void circle(int x, int y, int radius, Window w);
void xline(int x0, int y0, int x1, int y1, Window w);
int new_float(char *name, double *value);
int new_int(char *name, int *value);
void display_command(char *name, char *value, int pos, int col);
void clr_line_at(Window w, int col0, int pos, int n);
void put_cursor_at(Window w, int col0, int pos);
void put_string_at(Window w, int col, char *s, int off);
void movmem(char *s1, char *s2, int len);
void memmov(char *s1, char *s2, int len);
void edit_window(Window w, int *pos, char *value, int *col, int *done, int ch);
void do_backspace(int *pos, char *value, int *col, Window w);
void edit_command_string(XEvent ev, char *name, char *value, int *done, int *pos, int *col);
int new_string(char *name, char *value);


#endif
 
