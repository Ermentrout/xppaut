#ifndef XPPAUT_USERBUT_H
#define XPPAUT_USERBUT_H

#include <X11/Xlib.h>

/* --- Types --- */
typedef struct {
  Window w;
  char bname[10];
  int com;
} USERBUT;

/* --- Functions --- */
void user_button_events(XEvent report);
void user_button_press(Window w);
void user_button_draw(Window w);
void user_button_cross(Window w, int b);
int get_button_info(char *s, char *bname, char *sc);
int find_kbs(char *sc);
void add_user_button(char *s);
void create_user_buttons(int x0, int y0, Window base);

#endif /* XPPAUT_USERBUT_H */
