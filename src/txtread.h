#ifndef XPPAUT_TXTREAD_H
#define XPPAUT_TXTREAD_H

#include <X11/Xlib.h>

/* --- Types ---*/
typedef struct {
	Window up,down,pgup,pgdn,kill,home,end,base,text,src,action;
	int here,first,hgt,wid,nlines,which;
	int dh,dw;
} TXTVIEW;

void txt_view_events(XEvent ev);
void txtview_keypress(XEvent ev);
void enter_txtview(Window w, int val);
void do_txt_action(char *s);
void resize_txtview(int w, int h);
void txtview_press(Window w, int x, int y);
void redraw_txtview(Window w);
void redraw_txtview_text(void);
void init_txtview(void);
void make_txtview(void);

#endif /* XPPAUT_TXTREAD_H */
