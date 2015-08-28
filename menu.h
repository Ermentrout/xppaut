#ifndef _xppmenu_h_
#define _xppmenu_h_

#include <X11/Xlib.h>


void flash(int num);
void add_menu(Window base, int j, int n, char **names, char *key, char **hint);
void create_the_menus(Window base);
void show_menu(int j);
void unshow_menu(int j);
void help(void);
void help_num(void);
void help_file(void);
void menu_crossing(Window win, int yn);
void menu_expose(Window win);
void menu_button(Window win);
void draw_help(void);

#endif
