#ifndef XPPAUT_USERBUT_H
#define XPPAUT_USERBUT_H

#include <X11/Xlib.h>


/* --- Functions --- */
void add_user_button(char *s);
void create_user_buttons(int x0, int y0, Window base);
void user_button_events(XEvent report);

#endif /* XPPAUT_USERBUT_H */
