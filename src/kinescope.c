/*    Kinescope for X  windows       */
#include "kinescope.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>

#include "aniparse.h"
#include "browse.h"
#include "ggets.h"
#include "load_eqn.h"
#include "main.h"
#include "mykeydef.h"
#include "pop_list.h"
#include "scrngif.h"
#include "util/timeutil.h"

/* --- Macros --- */
#define MAXFILM 250


/* --- Types --- */
typedef struct {
	unsigned int h,w;
	Pixmap xi;
} MOVIE;


/* --- Forward declarations --- */
static void auto_play(void);
static void make_anigif(void);
static void play_back(void);
static void save_kine(void);
static void save_movie(char *basename, int fmat);
static int show_frame(int i, int h, int w);
static void too_small(void);


/* --- Data --- */
static int ks_ncycle=1;
static int ks_speed=50;
static MOVIE movie[MAXFILM];

int mov_ind;


/* --- Functions --- */
void do_movie_com(int c) {
	switch(c) {
	case 0:
		if(film_clip()==0) {
			respond_box("Okay","Out of film!");
		}
		break;
	case 1:
		reset_film();
		break;
	case 2:
		play_back();
		break;
	case 3:
		auto_play();
		break;
	case 4:
		save_kine();
		break;
	case 5:
		make_anigif();
		break;
	case 6: /* test_keys(); */
		break;
	}
}


int film_clip(void) {
	int x,y;
	unsigned int h,w,bw,d;
	Window root;
	if(mov_ind>=MAXFILM) {
		return(0);
	}
	XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
	movie[mov_ind].h=h;
	movie[mov_ind].w=w;
	movie[mov_ind].xi=XCreatePixmap(display,RootWindow(display,screen),w,h,
									DefaultDepth(display,screen));
	XCopyArea(display,draw_win,movie[mov_ind].xi,gc_graph,0,0,w,h,0,0);
	mov_ind++;
	return 1;
}


void reset_film(void) {
	int i;
	if(mov_ind==0) {
		return;
	}
	for(i=0;i<mov_ind;i++) {
		XFreePixmap(display,movie[i].xi);
	}
	mov_ind=0;
}


/* --- Static functions --- */
static void auto_play(void) {
	int x,y;
	unsigned int h,w,bw,d,key;
	Window root;

	int dt=20;
	int smax=500;
	XEvent ev;
	int i=0,cycle=0;

	new_int("Number of cycles",&ks_ncycle);
	new_int("Msec between frames",&ks_speed);
	if(ks_speed<0) {
		ks_speed=0;
	}
	if(ks_ncycle<=0) {
		return;
	}
	XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
	if(mov_ind==0) {
		return;
	}
	if(h<movie[i].h || w<movie[i].w) {
		too_small();
		return;
	}
	XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
	XFlush(display);
	while(1) {
		/* check for events    */
		if(XPending(display)>0) {
			XNextEvent(display,&ev);
			switch(ev.type) {
			case ButtonPress:
				return;
				break;
			case  KeyPress:
				key=get_key_press(&ev);
				if(key==ESC) {
					return;
				}
				if(key==',') {
					ks_speed-=dt;
					if(ks_speed<dt) {
						ks_speed=dt;
					}
				}
				if(key=='.') {
					ks_speed+=dt;
					if(ks_speed>smax) {
						ks_speed=smax;
					}
				}
				break;
			}
		}  /* done checking  now increment pix   */

		waitasec(ks_speed);
		i++;
		if(i>=mov_ind) {
			cycle++;
			i=0;
		}
		if(h<movie[i].h || w<movie[i].w) {
			too_small();
			return;
		}
		XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
		XFlush(display);
		if(cycle>=ks_ncycle) {
			return;
		}
	}  /*  Big loop   */
}


static void make_anigif(void) {
	int i=0;
	int x,y;
	FILE *fp;
	Window root;
	unsigned int h,w,bw,d;
	XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
	if(mov_ind==0) {
		return;
	}
	if(h<movie[i].h || w<movie[i].w) {
		too_small();
		return;
	}
	h=movie[0].h;
	w=movie[0].w;
	for(i=0;i<mov_ind;i++) {
		if((movie[i].h!=h) || (movie[i].w!=w)) {
			err_msg("All clips must be same size");
			return;
		}
	}
	fp=fopen("anim.gif","wb");
	set_global_map(1);
	for(i=0;i<mov_ind;i++) {
		XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
		XFlush(display);
		add_ani_gif(movie[i].xi,fp,i);
	}

	end_ani_gif(fp);
	fclose(fp);
	set_global_map(0);
}


static void play_back(void) {
	int x,y;
	unsigned int h,w,bw,d;

	Window root;
	XEvent ev;
	int i=0;
	XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
	if(mov_ind==0) {
		return;
	}
	if(h<movie[i].h || w<movie[i].w) {
		too_small();
		return;
	}

	XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
	XFlush(display);
	while(1) {
		XNextEvent(display,&ev);
		switch(ev.type) {
		case ButtonPress:
			i++;
			if(i>=mov_ind) {
				i=0;
			}
			if(show_frame(i,h,w)) {
				return;
			}
			break;
		case KeyPress:
			switch(get_key_press(&ev)) {
			case ESC:
				return;
			case RIGHT:
				i++;
				if(i>=mov_ind) {
					i=0;
				}
				if(show_frame(i,h,w)) {
					return;
				}
				break;
			case  LEFT:
				i--;
				if(i<0) {
					i=mov_ind-1;
				}
				if(show_frame(i,h,w)) {
					return;
				}
				break;
			case  HOME:
				i=0;
				if(show_frame(i,h,w)) {
					return;
				}
				break;
			case END:
				i=mov_ind-1;
				if(show_frame(i,h,w)) {
					return;
				}
				break;
			}
		}
	}
}


static void save_kine(void) {
	char base[128];
	int fmat=2;
	sprintf(base,"frame");
	new_string("Base file name",base);
	if(strlen(base)>0) {
		save_movie(base,fmat);
	}
}


static void save_movie(char *basename, int fmat) {
	char file[XPP_MAX_NAME];
	int i=0;
	int x,y;
	FILE *fp;
	Window root;
	Pixmap xi;
	unsigned int h,w,bw,d;
	XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
	if(mov_ind==0) {
		return;
	}
	if(h<movie[i].h || w<movie[i].w) {
		too_small();
		return;
	}
	for(i=0;i<mov_ind;i++) {
		if(fmat==1) {
			sprintf(file,"%s_%d.ppm",basename,i);
		} else {
			sprintf(file,"%s_%d.gif",basename,i);
		}
		XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
		XFlush(display);
		if(fmat==1) {
			writeframe(file,draw_win,w,h);
		}
#ifndef NOGIF
		else{
			XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
			xi=XCreatePixmap(display,RootWindow(display,screen),w,h,
							 DefaultDepth(display,screen));
			XCopyArea(display,draw_win,xi,gc_graph,0,0,w,h,0,0);

			fp=fopen(file,"wb");
			screen_to_gif(xi,fp);
			fclose(fp);
		}
#endif
	}
}


static int show_frame(int i, int h, int w) {
	if(h<movie[i].h || w<movie[i].w) {
		too_small();
		return 1;
	}
	XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
	XFlush(display);
	return 0;
}


static void too_small(void) {
	respond_box("Okay","Window too small for film!");
}
