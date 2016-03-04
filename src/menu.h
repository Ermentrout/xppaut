#ifndef XPPAUT_MENU_H
#define XPPAUT_MENU_H

#include <X11/Xlib.h>

/* --- Data --- */
extern int help_menu;

/* --- Functions --- */
void create_the_menus(Window base);
void draw_help(void);
void help(void);
void help_file(void);
void help_num(void);
void menu_button(Window win);
void menu_crossing(Window win, int yn);
void menu_expose(Window win);

#endif /* XPPAUT_MENU_H */
