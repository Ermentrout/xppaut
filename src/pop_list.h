#ifndef XPPAUT_POP_LIST_H
#define XPPAUT_POP_LIST_H

#include <X11/Xlib.h>
#include "xpplim.h"


/* --- Macros --- */
#define ALL_DONE 2
#define FORGET_ALL 0
#define DONE_THIS 1
#define DONE_ALL 2
#define FORGET_THIS 3

#define SB_PLOTTABLE 0
#define SB_VARIABLE 1
#define SB_PARAMETER 2
#define SB_PARVAR 3
#define SB_COLOR 4
#define SB_MARKER 5
#define SB_METHOD 6

#define EV_MASK    (ButtonPressMask 	|\
	KeyPressMask		|\
	ExposureMask		|\
	StructureNotifyMask)

#define BUT_MASK   (ButtonPressMask 	|\
	KeyPressMask		|\
	ExposureMask		|\
	StructureNotifyMask	|\
	EnterWindowMask		|\
	LeaveWindowMask)



/* --- Functions --- */
void do_hilite_text(char *name, char *value, int flag, Window w, int pos, int col);
int do_string_box(int n, int row, int col, char *title, char **names, char values[][25], int maxchar);
Window make_fancy_window(Window root, int x, int y, int width, int height, int bw, int fc, int bc);
Window make_icon_window(Window root,int x,int y,int width,int height,int bw,int icx,int icy,unsigned char* icdata);
Window make_plain_unmapped_window(Window root, int x, int y, int width, int height, int bw);
Window make_plain_window(Window root, int x, int y, int width, int height, int bw);
Window make_unmapped_window(Window root, int x, int y, int width, int height, int bw);
Window make_window(Window root, int x, int y, int width, int height, int bw);
void make_scrbox_lists(void);
void message_box(Window *w, int x, int y, char *message);
int pop_up_list(Window *root, char *title, char **list, char *key, int n, int max, int def, int x, int y, char **hints, Window hwin, char *httxt);
void respond_box(const char *button, const char *message);
void set_window_title(Window win, char *string);
int two_choice(char *choice1, char *choice2, char *string, char *key, int x, int y, Window w,char *title);
int yes_no_box(void);

#endif /* XPPAUT_POP_LIST_H */
