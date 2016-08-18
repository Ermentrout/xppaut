#include "abort.h"

#include <stdio.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "aniparse.h"
#include "ggets.h"
#include "main.h"
#include "many_pops.h"
#include "mykeydef.h"

int get_command_width(void) {
	int x,y;
	unsigned int w,h,bw,de;

	Window root;
	XGetGeometry(display,command_pop,&root,&x,&y,&w,&h,&bw,&de);
	XClearWindow(display,command_pop);
	return(w);
}

void plot_command(int nit, int icount, int cwidth) {
	int i;
	float dx;
	if(nit==0) {
		return;
	}
	dx=(float)icount*(float)cwidth/(float)nit;
	i=(int)dx;

	XDrawPoint(display,command_pop,gc,i,5);
}


int my_abort(void) {
	int ch;
	XEvent event;
	while(XPending(display)>0) {
		XNextEvent(display,&event);
		if(check_ani_pause(event)==ESC) {
			return(ESC);
		}
		switch(event.type) {
		case Expose:
			do_expose(event);
			break;
		case ButtonPress:
			break;
		case KeyPress:
			ch=get_key_press(&event);
			return(ch);
		}
		return(0);
	}
	return(64);
}
