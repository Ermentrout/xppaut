#include "choice_box.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/cursorfont.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>

#include "ggets.h"
#include "main.h"
#include "struct.h"

/* --- Functions --- */
void display_choice(w,p)
Window w;
CHOICE_BOX p;
{
	int i;
	int n=p.n;
	XSetFillStyle(display,gc,FillSolid);
	XSetForeground(display,gc,MyForeColor);

	if(w==p.ok)XDrawString(display,w,gc,0,CURY_OFF,"Ok",2);
	if(w==p.cancel)
	XDrawString(display,w,gc,0,CURY_OFF,"Cancel",6);
	for(i=0;i<n;i++)
	{
		if(w!=p.cw[i])continue;
		XDrawString(display,w,gc,0,CURY_OFF,p.name[i],strlen(p.name[i]));
		if(p.flag[i]==1)set_fore();
		else set_back();
		XDrawString(display,w,gc,(p.mc+1)*DCURX,CURY_OFF,"X",1);
	}
	set_fore();
}
