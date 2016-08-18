#ifndef XPPAUT_GGETS_H
#define XPPAUT_GGETS_H

#include <X11/Xlib.h>

/* --- Macros --- */
#define ClickTime 200

/* --- Data --- */
extern char *info_message;
extern int MSStyle;
extern int xor_flag;

/* --- Functions --- */
void bar(int x, int y, int x2, int y2, Window w);
void blank_screen(Window w);
void bottom_msg(int line, char *msg);
void chk_xor(void);
void clr_command(void);
void cput_text(void);
void display_command(char *name, char *value, int pos, int col);
void edit_command_string(XEvent ev, char *name, char *value, int *done, int *pos, int *col);
void edit_window(Window w, int *pos, char *value, int *col, int *done, int ch);
void err_msg(const char *string);
void Ftext(int x, int y, const char *string, Window o);
int get_key_press(XEvent *ev);
int get_mouse_xy(int *x, int *y, Window w);
int new_float(char *name, double *value);
int new_int(char *name, int *value);
int new_string(char *name, char *value);
void ping(void);
int plintf(char *fmt, ...);
void put_command(char *string);
void put_cursor_at(Window w, int col0, int pos);
void rectangle(int x, int y, int x2, int y2, Window w);
void reset_graphics(void);
void set_back(void);
void set_fore(void);
int show_position(XEvent ev, int *com);
void showchar(int ch, int col, int row, Window win);
void xline(int x0, int y0, int x1, int y1, Window w);

#endif /* XPPAUT_GGETS_H */
