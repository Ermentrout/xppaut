#ifndef _choice_box_h_
#define _choice_box_h_

#include <X11/Xlib.h>
#include "struct.h"

void destroy_choice(CHOICE_BOX p);
void display_choice(Window w, CHOICE_BOX p);
void do_checks(CHOICE_BOX p);
void base_choice(char *wname, int n, int mcc, char **names, int *check, int type);
int do_choice_box(Window root, char *wname, int n, int mcc, char **names, int *check, int type);
int choice_box_event_loop(CHOICE_BOX p);

#endif
