
#ifndef _dialog_box_h
#define _dialog_box_h

#include <X11/Xlib.h>
#include "struct.h"

int get_dialog(char *wname, char *name, char *value, char *ok, char *cancel, int max);
int dialog_event_loop(DIALOG *d, int max, int *pos, int *col);
void display_dialog(Window w, DIALOG d, int pos, int col);


#endif
