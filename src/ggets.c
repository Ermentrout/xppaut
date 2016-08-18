#include "ggets.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <X11/XKBlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "axes2.h"
#include "calc.h"
#include "color.h"
#include "graphics.h"
#include "load_eqn.h"
#include "main.h"
#include "many_pops.h"
#include "menudrive.h"
#include "mykeydef.h"
#include "newhome.h"
#include "pop_list.h"
#include "util/strutil.h"
#include "util/timeutil.h"

/* --- Forward declarations --- */
static void clr_line_at(Window w, int col0, int pos, int n);
static void put_string_at(Window w, int col, char *s, int off);


/* --- Data --- */
static int CURS_X,CURS_Y;

char *info_message;
int MSStyle=0;
int xor_flag;


/* --- Functions --- */
void ping(void) {
	if(tfBell && !XPPBatch) {
		/*
		XkbBell allows window managers to react
		to bell events using possibly user-specified
		accessiblity options (e.g. visual bell)
		*/
		XkbBell(display,command_pop,100,(Atom)NULL);
	}
}


void reset_graphics(void) {
	blank_screen(draw_win);
	do_axes();
	hi_lite(draw_win);
}


void blank_screen(Window w) {
	CURS_X=0;
	CURS_Y=0;
	xor_flag=0;
	XClearWindow(display,w);
}

void set_fore(void) {
	XSetForeground(display,gc,MyForeColor);
}


void set_back(void) {
	XSetForeground(display,gc,MyBackColor);
}

void showchar(int ch, int col, int row, Window win) {
	char bob[2];
	bob[0]=ch;
	chk_xor();
	XDrawString(display,win,gc,col,row+CURY_OFF,bob,1);
}

void chk_xor(void) {
	if(xor_flag==1) {
		XSetFunction(display,gc,GXxor);
	} else {
		XSetFunction(display,gc,GXcopy);
	}
}


void clr_command(void) {
	blank_screen(command_pop);
}


void bottom_msg(int line, char *msg) {
	XClearWindow(display,info_pop);
	BaseCol();
	strcpy(info_message,msg);
	XDrawString(display,info_pop,gc,5,CURY_OFF,msg,strlen(msg));
}

void err_msg(const char *string) {
	if(Xup) {
		respond_box("OK",string);
	} else {
		plintf("%s\n",string);
	}
}


int plintf(char *fmt,...) {
	int nchar=0;
	va_list arglist;

	if(!XPPVERBOSE) {
		return(nchar);/*Don't print at all!*/
	}
	if(logfile == NULL) {
		printf("The log file is NULL!\n");
		logfile = stdout;
	}

	va_start(arglist,fmt);
	nchar=vfprintf(logfile,fmt,arglist);
	va_end(arglist);
	/*Makes sense to flush to the output file to
	prevent loss of log info if program crashes.
	Then maybe user can figure out what happened and when.*/
	fflush(logfile);

	return(nchar);
}

int show_position(XEvent ev, int *com) {
	check_draw_button(ev);
	return(0);
}


void put_command(char *string) {
	clr_command();
	Ftext(0,0,string,command_pop);
	CURS_X=strlen(string);
}


int get_key_press(XEvent *ev) {
	int maxlen=64;
	char buf[65];
	XComposeStatus comp;
	KeySym ks;

	XLookupString((XKeyEvent *)ev,buf,maxlen,&ks,&comp);
	/*       printf(" ks=%d buf[0]=%d char=%c \n",ks,(int)buf[0],buf[0]); */

	/*    plintf("h=%d e=%d ks=%d \n",XK_Home,XK_End,ks); */
	if(ks==XK_Escape) {
		return(ESC);
	}
	if((ks==XK_Return)  ||  (ks==XK_KP_Enter)  ||  (ks==XK_Linefeed)) {
		return(FINE);
	} else if(((ks>=XK_KP_Space) && (ks<=XK_KP_9))  ||
			  ((ks>=XK_space) && (ks<=XK_asciitilde))) {
		return((int)buf[0]);
	} else if(ks==XK_BackSpace) {
		return(BACKSPACE);
	} else if(ks==XK_Delete) {
		return(DELETE);
	} else if(ks==XK_Tab) {
		return(TAB);
	} else if(ks==XK_Home) {
		return(HOME);
	} else if(ks==XK_End) {
		return(END);
	} else if(ks==XK_Left) {
		return(LEFT);
	} else if(ks==XK_Right) {
		return(RIGHT);
	} else if(ks==XK_Up) {
		return(UP);
	} else if(ks==XK_Down) {
		return(DOWN);
	} else if(ks==XK_PgUp) {
		return(PGUP);
	} else if(ks==XK_PgDn) {
		return(PGDOWN);
	} else {
		return(BADKEY);
	}
}


void cput_text(void) {
	char string[DEFAULT_STRING_LENGTH],new_char[256];
	int x,y,size=2,font=0;
	Window temp;
	temp=main_win;
	strcpy(string,"");
	if(new_string("Text: ",string)==0) {
		return;
	}
	if(string[0]=='%') {
		fillintext(&string[1],new_char);
		strcpy(string,new_char);
	} /* this makes it permanent */

	new_int("Size 0-4 :",&size);
	/* new_int("Font  0-times/1-symbol :",&font); */
	if(size>4) {
		size=4;
	}
	if(size<0) {
		size=0;
	}
	message_box(&temp,0,SCALEY-5*DCURY,"Place text with mouse");
	if(GetMouseXY(&x,&y)) {
		GrCol();
		fillintext(string,new_char);
		special_put_text_x11(x,y,new_char,size);
		add_label(string,x,y,size,font);
		BaseCol();
	}
	waitasec(ClickTime);
	XDestroyWindow(display,temp);
}


int get_mouse_xy(int *x, int *y, Window w) {
	int no_but=1;
	char ch;
	XEvent ev;
	*x=0;
	*y=0;
	while(no_but) {
		XNextEvent(display,&ev);
		switch(ev.type) {
		case Expose:
			do_expose(ev);
			break;
		case KeyPress:
			ch=get_key_press(&ev);
			if(ch==ESC) {
				return 0;
			}
			if(ch==FINE) {
				return -2;
			}
			if(ch==TAB) {
				return -3;
			}
			break;
		case ButtonPress:
			if(ev.xbutton.window!=w) {
				return(0);
			}
			no_but=0;
			*x=ev.xbutton.x;
			*y=ev.xbutton.y;
			return(1);
		}
	}
	return(0);
}


void Ftext(int x, int y, const char *string, Window o) {
	chk_xor();
	XDrawString(display,o,gc,x,y+CURY_OFF,string,strlen(string));
}


void bar(int x, int y, int x2, int y2, Window w) {
	XFillRectangle(display,w,gc,x,y,x2-x,y2-y);
}


void rectangle(int x, int y, int x2, int y2, Window w) {
	XDrawRectangle(display,w,gc,x,y,x2-x,y2-y);
}


void xline(int x0, int y0, int x1, int y1, Window w) {
	XDrawLine(display,w,gc_graph, x0,y0,x1,y1);
}

int new_float(char *name, double *value) {
	int done;
	int flag;
	double newz;
	char tvalue[200];
	sprintf(tvalue,"%.16g",*value);
	done=new_string(name,tvalue);
	if(done==0 || strlen(tvalue)==0) {
		return -1;
	}

	if(tvalue[0]=='%') {
		flag=do_calc(&tvalue[1],&newz);
		if(flag!=-1) {
			*value=newz;
		}
		return(0);
	}
	*value=atof(tvalue);
	return(0);
}


int new_int(char *name, int *value) {
	char svalue[200];
	sprintf(svalue,"%d",*value);
	if(new_string(name,svalue)==0 || strlen(svalue)==0) {
		return(-1);
	}

	*value=atoi(svalue);
	return(0);
}


void display_command(char *name, char *value, int pos, int col) {
	int l=strlen(name);
	int m=strlen(value);

	set_fore();
	bar(0,0,l*DCURX,DCURY+4,command_pop);
	set_back();
	XDrawString(display,command_pop,gc,0,CURY_OFF,name,l);
	set_fore();
	if(m>0) {
		XDrawString(display,command_pop,gc,l*DCURX,CURY_OFF,value,m);
		/* showchar('_',DCURX*(l+m),0,command_pop); */
		put_cursor_at(command_pop,DCURX*l,pos);
	}
}


void put_cursor_at(Window w, int col0, int pos) {
	int x1=col0+pos*DCURX;
	int x2=x1+1;
	int y1=DCURY-2,y2=2;
	/* XDrawString(display,w,gc,col0+pos*DCURX-1,DCURY,"^",1);*/
	XDrawLine(display,w,gc,x1,y1,x1,y2);
	XDrawLine(display,w,gc,x2,y1,x2,y2);
}


void edit_window(Window w, int *pos, char *value, int *col, int *done, int ch) {
	int col0=*col-*pos*DCURX;
	*done=0;
	/* plintf(" po=%d cl=%d ch=%d  || %s ||  c0=%d\n",*pos,*col,ch,value,col0); */
	switch(ch) {
	case LEFT:
		if(*pos>0) {
			*pos=*pos-1;
			*col-=DCURX;
		} else {
			ping();
		}
		break;
	case RIGHT:
		if(*pos<strlen(value)) {
			*pos=*pos+1;
			*col+=DCURX;
		} else {
			ping();
		}
		break;
	case HOME:
		*pos=0;
		*col=col0;
		break;
	case END:
		*pos=strlen(value);
		*col=*pos*DCURX+col0;
		break;
	case BADKEY:
	case DOWN:
	case UP:
	case PGUP:
	case PGDOWN:
		return;    /* junk key  */
	case ESC:
		*done=-1;  /* quit without saving */
		return;
	case FINE:
		if(MSStyle==0) {
			*done=1;
		} else {
			*done=2;
		}
		return;   /* save this guy */
	case BACKSPACE:
	case DELETE:
		if(*pos>0) {
			memmov(&value[*pos-1],&value[*pos],strlen(value)-*pos+1);
			*pos=*pos-1;
			*col-=DCURX;
		} else {
			ping();
		}
		break;
	case TAB:
		if(MSStyle==0) {
			*done=2;
		} else {
			*done=1;
		}
		return;
	default:
		if( (ch>=' ') && (ch <= '~')) {
			movmem(&value[*pos+1],&value[*pos],strlen(value)-*pos+1);
			value[*pos]=ch;
			*pos=*pos+1;
			*col+=DCURX;
		}
		break;
	} /* end key cases */
	/* Now redraw everything !!  */
	clr_line_at(w,col0,0,strlen(value));
	put_string_at(w,col0,value,0);
	put_cursor_at(w,col0,*pos);
	XFlush(display);
}


void edit_command_string(XEvent ev, char *name, char *value, int *done, int *pos, int *col) {
	char ch;
	switch(ev.type) {
	case ConfigureNotify:
	case Expose:
	case MapNotify:
		do_expose(ev);
		if(ev.xexpose.window==command_pop) {
			display_command(name,value,*pos,0);
		}
		break;
	case ButtonPress:
		if(ev.xbutton.window==command_pop) {
			XSetInputFocus(display,command_pop,
						   RevertToParent,CurrentTime);
		}
		break;
	case KeyPress:
		ch=get_key_press(&ev);
		edit_window(command_pop,pos,value,col,done,ch);
	} /* end event cases */
}


int new_string(char *name, char *value) {
	char old_value[80];
	int done=0;
	int pos=strlen(value);
	int col=(pos+strlen(name))*DCURX;

	XEvent ev;
	strcpy(old_value,value);
	clr_command();
	display_command(name,value,pos,0);
	while(done==0) {
		XNextEvent(display,&ev);
		edit_command_string(ev,name,value,&done,&pos,&col);
	}
	clr_command();
	if(done==1 || done==2) {
		return(done);
	}
	strcpy(value,old_value);
	return(0);
}

/* --- Static functions --- */
static void clr_line_at(Window w, int col0, int pos, int n) {
	XClearArea(display,w,col0+pos*DCURX,0,(n+2)*DCURX,2*DCURY,False);
}


static void put_string_at(Window w, int col, char *s, int off) {
	int l=strlen(s)-off;
	XDrawString(display,w,gc,col,CURY_OFF,s+off,l);
}
