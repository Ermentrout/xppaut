/***************   NOTES ON MPEG STUFF   ********************
To prepare for mpeg encoding in order to make your movies
permanent, I have to do some image manipulation - the main
routine is writeframe()

The current version works for most 8 bit color servers.  I have
a version also working for TrueColor 16 bit and I think it works on
24 bit color as well but havent tried it.  I really dont know
how all colors are organized.  For my machine the 15 lowest order bits
code color as
	 xrrrrrgggggbbbbb
in binary so lobits are blue etc. If the colors seem screwy, then you might
want to alter the ordering below

************************************************************/
#include "aniparse.h"

#ifndef HAVE_WCTYPE_H
# include <ctype.h>
#else
# include <wctype.h>
#endif
#include <fcntl.h>
#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <X11/Xlib.h>
#include <X11/Xproto.h>
#include <X11/Xutil.h>

#include "browse.h"
#include "color.h"
#include "dialog_box.h"
#include "form_ode.h"
#include "ggets.h"
#include "graf_par.h"
#include "graphics.h"
#include "init_conds.h"
#include "integrate.h"
#include "load_eqn.h"
#include "main.h"
#include "many_pops.h"
#include "markov.h"
#include "menudrive.h"
#include "mykeydef.h"
#include "my_rhs.h"
#include "nullcline.h"
#include "parserslow.h"
#include "pop_list.h"
#include "scrngif.h"
#include "util/strutil.h"
#include "toons.h"
#include "util/timeutil.h"
#include "bitmap/aniwin.bitmap"

/* --- Macros --- */
/* maximum variables you can change per grabbable */
#define MAX_GEVENTS 20
/* max grabbable objects  */
#define MAX_ANI_GRAB 50
#define INIT_C_SHIFT 0

/* who knows how the colors are ordered */
#ifdef BGR
#define MY_BLUE hibits
#define MY_GREEN midbits
#define MY_RED lobits
#else
#define MY_BLUE lobits
#define MY_GREEN midbits
#define MY_RED hibits
#endif

#define LINE 0
#define RLINE 1
#define CIRC 2
#define FCIRC 3
#define RECT 4
#define FRECT 5
#define TEXT 6
#define VTEXT 7
#define ELLIP 9
#define FELLIP 10
#define COMET 11
#define PCURVE 12
#define AXNULL 13
#define AYNULL 14
#define GRAB 25
/*  not for drawing */

#define SETTEXT 8

/*  Not in command list   */

#define TRANSIENT 20
#define PERMANENT 21
#define DIMENSION 22
#define COMNT 30
#define SPEED 23


/* --- Types --- */
typedef struct {
	int n;
	int *x,*y,*col;
	int i;
} Comet;

typedef struct {  /* tasks have the form {name1=formula1;name2=formula2;...} */
	double vrhs[MAX_GEVENTS];
	char lhsname[MAX_GEVENTS][11];
	int lhsivar[MAX_GEVENTS];
	int *comrhs[MAX_GEVENTS];
	int runnow;
	int n; /* number of tasks <= MAX_GEVENTS */
}GRAB_TASK;

typedef struct {
	double x0,y0;
	double x,y;
	double ox,oy;
	double t1,t2,tstart;
	double vx,vy;
	double vax,vay;
} ANI_MOTION_INFO;

typedef struct {
	int flag;
	int skip;
	char root[100];
	char filter[DEFAULT_STRING_LENGTH];
	int aviflag,filflag;
} MPEG_SAVE;

typedef struct {
	Window base, wfile,wgo,wpause,wreset,wfast,wslow,wmpeg;
	Window wfly,kill,slider;
	Window wup,wdn,wskip;
	Window view,wgrab;
	int hgt,wid,iexist,ok;
	int pos,inc;
	int slipos,sliwid;
	char file[XPP_MAX_NAME];
} VCR;

typedef struct {
	Comet c;
	int type, flag;
	int *col,*x1,*y1,*x2,*y2,*who;
	double zcol,zx1,zy1,zx2,zy2,zrad,zval;
	int zthick,tfont,tsize,tcolor;
} ANI_COM;

typedef struct {
	int ok;
	double zx,zy,tol;
	int *x,*y;
	GRAB_TASK start,end;
} ANI_GRAB;


/* --- Forward Declarations --- */
static int add_ani_circle(ANI_COM *a, char *x1, char *y1, char *x2, char *col,
						  char *thick);
static int add_ani_com(int type, char *x1, char *y1, char *x2, char *y2,
					   char *col, char *thick);
static int add_ani_comet(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						 char *col, char *thick);
static int add_ani_ellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						 char *col, char *thick);
static int add_ani_expr(char *x, int *c);
static int add_ani_fellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						  char *col, char *thick);
static int add_ani_frect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						 char *col, char *thick);
static int add_ani_line(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						char *col, char *thick);
static int add_ani_null(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						char *col, char *who);
static int add_ani_rect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						char *col, char *thick);
static int add_ani_rline(ANI_COM *a, char *x1, char *y1, char *col,
						 char *thick);
static int add_ani_settext(ANI_COM *a, char *x1, char *y1, char *col);
static int add_ani_text(ANI_COM *a, char *x1, char *y1, char *y2);
static int add_ani_vtext(ANI_COM *a, char *x1, char *y1, char *x2, char *y2);
static int add_grab_command(char *xs, char *ys, char *ts, FILE *fp);
static int add_grab_task(char *lhs, char *rhs, int igrab, int which);
static void ani_border(Window w, int i);
static void ani_button(Window w);
static void ani_buttonx(XEvent ev, int flag);
static void ani_create_mpeg(void);
static void ani_disk_warn(void);
static void ani_flip(void);
static void ani_flip1(int n);
static void ani_frame(int task);
static int ani_grab_tasks(char *line, int igrab, int which);
static void ani_ij_to_xy(int ix, int iy, double *x, double *y);
static void ani_motion_stuff(Window w, int x, int y);
static int ani_new_file(char *filename);
static void ani_newskip(void);
static void ani_rad2scale(double rx, double ry, int *ix, int *iy);
static void ani_radscale(double rad, int *ix, int *iy);
static void ani_resize(int x, int y);
static void ani_xyscale(double x, double y, int *ix, int *iy);
static void check_on_the_fly(void);
static int chk_ani_color(char *s, int *index);
static void create_vcr(char *name);
static void destroy_vcr(void);
static void do_ani_slider_motion(Window w, int x);
static void do_grab_tasks(int which);
static void draw_ani_circ(int j);
static void draw_ani_comet(int j);
static void draw_ani_ellip(int j);
static void draw_ani_fcirc(int j);
static void draw_ani_fellip(int j);
static void draw_ani_frect(int j);
static void draw_ani_line(int j);
static void draw_ani_null(int j, int id);
static void draw_ani_rect(int j);
static void draw_ani_rline(int j);
static void draw_ani_slider(Window w,int x);
static void draw_ani_text(int j);
static void draw_ani_vtext(int j);
static void draw_grab_points(void);
static void eval_ani_color(int j);
static void eval_ani_com(int j);
static void free_ani(void);
static void free_grabber(void);
static double get_current_time(void);
static void init_ani_stuff(void);
static int load_ani_file(FILE *fp);
static int parse_ani_string(char *s, FILE *fp);
static void read_ani_line(FILE *fp, char *s);
static void redraw_ani_slider(void);
static void render_ani(void);
static void reset_comets(void);
static void roll_comet(ANI_COM *a, int xn, int yn, int col);
static int run_now_grab(void);
static int search_for_grab(double x, double y);
static void set_ani_col(int j);
static void set_ani_dimension(char *x1, char *y1, char *x2, char *y2);
static void set_ani_font_stuff(int size, int font, int color);
static void set_ani_perm(void);
static void set_ani_thick(int t);
static void set_from_init_data(void);
static void set_to_init_data(void);
static void tst_pix_draw(void);
static void update_ani_motion_stuff(int x, int y);
static void xset_ani_col(int icol);

/* --- Data --- */
static int n_ani_grab = 0;
static int show_grab_points = 0;
static int ani_grab_flag = 0;
static int who_was_grabbed;

static int on_the_fly_speed = 10;

static int aniflag;
static int LastAniColor;

static int n_anicom;
static int ani_line;
static int ani_speed = 10;
static int ani_speed_inc = 2;
static int ani_text_size;
static int ani_text_color;
static int ani_text_font;

static double ani_xlo = 0, ani_xhi = 1, ani_ylo = 0, ani_yhi = 1;
static double ani_lastx, ani_lasty;
static Pixmap ani_pixmap;

static ANI_COM my_ani[MAX_ANI_LINES];
static ANI_GRAB ani_grab[MAX_ANI_GRAB];
static ANI_MOTION_INFO ami;
static GC ani_gc;
static MPEG_SAVE mpeg;
static VCR vcr;

int animation_on_the_fly = 0;

/* --- Functions --- */
void ani_expose(Window w) {
	if(vcr.iexist==0) {
		return;
	}
	if(w==vcr.wgrab) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Grab",4);
	}
	if(w==vcr.view) {
		XCopyArea(display,ani_pixmap,vcr.view,ani_gc,0,0,vcr.wid,vcr.hgt,0,0);
	}
	if(w==vcr.wgo) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Go  ",4);
	}
	if(w==vcr.wup) {
		XDrawString(display,w,small_gc,5,CURY_OFFs," >>>>",5);
	}
	if(w==vcr.wskip) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Skip",4);
	}
	if(w==vcr.wdn) {
		XDrawString(display,w,small_gc,5,CURY_OFFs," <<<<",5);
	}
	if(w==vcr.wfast) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Fast",4);
	}
	if(w==vcr.wslow) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Slow",4);
	}
	if(w==vcr.slider) {
		draw_ani_slider(w,vcr.slipos);
	}
	if(w==vcr.wpause) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Pause",5);
	}
	if(w==vcr.wreset) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Reset",5);
	}
	if(w==vcr.kill) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Close",5);
	}
	if(w==vcr.wfile) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"File",4);
	}
	if(w==vcr.wmpeg) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"MPEG",4);
	}
	if(w==vcr.wfly) {
		check_on_the_fly();
	}
}


void ani_zero(void) {
	vcr.iexist=0;
	vcr.ok=0;
	vcr.inc=1;
	vcr.pos=0;
	n_anicom=0;
	ani_speed=10;
	aniflag=TRANSIENT;
	ani_grab_flag=0;
	if(use_ani_file) {
		strcpy(vcr.file,anifile);
	} else {
		strcpy(vcr.file,this_file);
		sprintf(vcr.file,"%s/",dirname(vcr.file));
	}
}


int check_ani_pause(XEvent ev) {
	if((vcr.iexist==0) || (!animation_on_the_fly)) {
		return 0;
	}
	if(ev.type==ButtonPress  &&  ev.xbutton.window==vcr.wpause) {
		return(ESC);
	}
	return 0;
}


void do_ani_events(XEvent ev) {
	int x,y;

	if(vcr.iexist==0) {
		return;
	}
	switch(ev.type) {
	case ConfigureNotify:
		if(ev.xconfigure.window!=vcr.base) {
			return;
		}
		x=ev.xconfigure.width;
		y=ev.xconfigure.height;
		x=(x)/8;
		x=8*x;
		y=(y)/8;
		y=y*8;
		ani_resize(x,y);
		break;
	case EnterNotify:
		ani_border(ev.xexpose.window,2);
		break;
	case LeaveNotify:
		ani_border(ev.xexpose.window,1);
		break;
	case MotionNotify:
		do_ani_slider_motion(ev.xmotion.window,ev.xmotion.x);
		if(ani_grab_flag == 0) {
			break;
		}
		ani_motion_stuff(ev.xmotion.window,ev.xmotion.x,ev.xmotion.y);
		break;
	case ButtonRelease:
		if(ani_grab_flag==0) {
			break;
		}
		ani_buttonx(ev,0);
		break;
	case ButtonPress:
		ani_buttonx(ev,1);
		break;
	}
}


void get_ani_file(char *fname) {
	int status;
	int err;

	if (fname == NULL) {
		status=file_selector("Load animation",vcr.file,"*.ani");
		if(status==0) {
			return;
		}
	} else {
		strcpy(vcr.file,fname);
	}
	err=ani_new_file(vcr.file);
	if(err>=0) {
		vcr.ok=1; /* loaded and compiled */
		plintf("Loaded %d lines successfully!\n",n_anicom);
		ani_grab_flag=0;
	}
}


int getppmbits(Window window,int *wid,int *hgt, unsigned char *out) {
	XImage *ximage;
	Colormap cmap;
	unsigned long value;
	int i;
	int CMSK=0,CSHIFT=0,CMULT=0;
	int bbp=0,bbc=0;
	int lobits,midbits,hibits;
	unsigned x,y;
	XColor palette[DEFAULT_STRING_LENGTH];
	XColor pix;
	unsigned char *dst,*pixel;

	cmap = DefaultColormap(display,screen);
	ximage=XGetImage(display,window,0,0,*wid,*hgt,AllPlanes,ZPixmap);
	if(!ximage) {
		return -1;
	}
	/* this is only good for 256 color displays */
	for(i = 0; i < 256; i++) {
		palette[i].pixel = i;
	}
	XQueryColors(display, cmap, palette, 256);
	if(TrueColorFlag==1) {
		bbp=ximage->bits_per_pixel; /* is it 16 or 24 bit */
		if(bbp>24){
			bbp=24;
		}
		bbc=bbp/3;  /*  divide up the 3 colors equally to bbc bits  */
		CMSK=(1<<bbc)-1;  /*  make a mask  2^bbc  -1  */
		CSHIFT=bbc;       /*  how far to shift to get the next color */
		CMULT=8-bbc;       /* multiply 5 bit color to get to 8 bit */
	}
	*wid=ximage->width;
	*hgt=ximage->height;
	pixel=(unsigned char*)ximage->data;
	dst=out;
	for(y=0;y < (unsigned)(ximage->height); y++) {
		for (x = 0; x < (unsigned)(ximage->width); x++) {
			if(TrueColorFlag==1) {
				/*  use the slow way to get the pixel but then you dont need
				 * to screw around with byte order etc
				 */
				value=XGetPixel(ximage,x,y)>>INIT_C_SHIFT;
				/*vv=value; Not used?*/
				/*  get the 3 colors   hopefully  */
				lobits=value&CMSK;
				value=value>>CSHIFT;
				if(bbc==5) {
					value=value>>1;
				}
				midbits=value&CMSK;
				value=value>>CSHIFT;
				hibits=value&CMSK;

				/* store them for ppm dumping  */
				*dst++=(MY_RED<<CMULT);
				*dst++=(MY_GREEN<<CMULT);
				*dst++=(MY_BLUE<<CMULT);
			} else {
				/* 256 color is easier sort of  */
				pix = palette[*pixel++];
				*dst++ = pix.red;
				*dst++ = pix.green;
				*dst++ = pix.blue;
			}
		}
	}
	/* XDestroyImage(ximage); */
	return(1);
}


void new_vcr(void) {
	int tt,i;
	if(vcr.iexist==1) {
		return;
	}
	tt=gettimenow();
	i=(10+(tt%10))%10;
	if(i>=0 && i<10) {
		create_vcr(toons[i]);
	} else {
		create_vcr("Wanna be a member");
	}
}


void on_the_fly(int task) {
	if(vcr.iexist==0 || n_anicom==0) {
		return;
	}
	ani_frame(task);
	waitasec(on_the_fly_speed);
}


int writeframe(char *filename, Window window, int wid, int hgt) {
	int fd;
	XImage *ximage;
	Colormap cmap;
	unsigned long value;
	int i;
	int CMSK=0,CSHIFT=0,CMULT=0;
	int bbp=0,bbc=0;
	int lobits,midbits,hibits;

	unsigned x,y;
	char head[100];
	XColor palette[DEFAULT_STRING_LENGTH];
	XColor pix;
	unsigned char *pixel;
	unsigned area;
	unsigned char *out,*dst;

	cmap = DefaultColormap(display,screen);
	ximage=XGetImage(display,window,0,0,wid,hgt,AllPlanes,ZPixmap);
	if(!ximage) {
		return -1;
	}
	/* this is only good for 256 color displays */
	for(i = 0; i < 256; i++) {
		palette[i].pixel = i;
	}
	XQueryColors(display, cmap, palette, 256);
	fd=creat(filename,0666);
	if(fd==-1) {
		return -1;
	}
	/*    this worked for me - but you may want to change it for your machine */
	if(TrueColorFlag==1) {
		bbp=ximage->bits_per_pixel; /* is it 16 or 24 bit */
		if(bbp>24) {
			bbp=24;
		}
		bbc=bbp/3;  /*  divide up the 3 colors equally to bbc bits  */
		CMSK=(1<<bbc)-1;  /*  make a mask  2^bbc  -1  */
		CSHIFT=bbc;       /*  how far to shift to get the next color */
		CMULT=8-bbc;       /* multiply 5 bit color to get to 8 bit */
		/* plintf(" bbp=%d CMSK=%d CSHIFT=%d CMULT=%d \n",
	  bbp,CMSK,CSHIFT,CMULT); */
	}
	sprintf(head,"P6\n%d %d\n255\n",ximage->width,ximage->height);
	if(write(fd,head,strlen(head))<0) {
		plintf("Error writing file %s", fd);
	}
	area=ximage->width*ximage->height;
	pixel=(unsigned char*)ximage->data;
	out=(unsigned char *)malloc(3*area);
	dst=out;
	for(y=0;y < (unsigned)(ximage->height); y++) {
		for (x = 0; x < (unsigned)(ximage->width); x++) {
			if(TrueColorFlag==1) {
				/* use the slow way to get the pixel
		* but then you dont need to screw around
		* with byte order etc
	*/
				value=XGetPixel(ximage,x,y)>>INIT_C_SHIFT;
				/* vv=value*/
				/* get the 3 colors   hopefully  */
				lobits=value&CMSK;
				value=value>>CSHIFT;
				if(bbc==5) {
					value=value>>1;
				}
				midbits=value&CMSK;
				value=value>>CSHIFT;
				hibits=value&CMSK;
				/* store them for ppm dumping  */
				*dst++=(MY_RED<<CMULT);
				*dst++=(MY_GREEN<<CMULT);
				*dst++=(MY_BLUE<<CMULT);
			} else {
				/* 256 color is easier sort of  */
				pix = palette[*pixel++];
				*dst++ = pix.red;
				*dst++ = pix.green;
				*dst++ = pix.blue;
			}
		}
	}
	if(write(fd,out,area*3)<0) {
		plintf("Error writing file %s", fd);
	}
	close(fd);
	free(out);
	free(ximage);
	return 1;
}

/* --- Static functions --- */
static void create_vcr(char *name) {
	unsigned int valuemask=0;
	XGCValues values;
	Window base;
	int wid=280,hgt=350;
	XSizeHints size_hints;
	XTextProperty winname,iconname;

	base=make_plain_window(RootWindow(display,screen),0,0,5*12*DCURXs+8*DCURXs+4,20*(DCURYs+6),1);
	vcr.base=base;
	size_hints.flags=PPosition|PSize|PMinSize;
	size_hints.min_width=51*DCURXs;
	size_hints.min_height=300;
	XStringListToTextProperty(&name,1,&winname);
	XStringListToTextProperty(&name,1,&iconname);
	XSetWMProperties(display,base,&winname,&iconname,NULL,0,&size_hints,NULL,NULL);
	make_icon((char *)aniwin_bits,aniwin_width,aniwin_height,base);
	vcr.wfile   = br_button(base,0,0,"File",0);
	vcr.wgo = br_button(base,0,1,"Go",0);
	vcr.wreset = br_button(base,0,2,"Reset",0);
	vcr.wskip=br_button(base,0,3,"Skip",0);
	vcr.wfast   = br_button(base,1,0,"Fast",0);
	vcr.wslow = br_button(base,1,1,"Slow",0);
	vcr.wup = br_button(base,1,2,">>>>",0);
	vcr.wdn = br_button(base,1,3,"<<<<",0);
	vcr.wgrab=br_button(base,2,3,"Grab",0);
	vcr.slider=make_window(base,DCURXs,7+4*DCURYs,48*DCURXs,DCURYs+4,1);
	vcr.slipos=0;
	vcr.sliwid=48*DCURXs;
	vcr.wpause = br_button(base,2,0,"Pause",0);
	vcr.wmpeg = br_button(base,2,1,"MPeg",0);
	vcr.kill=br_button(base,2,2,"Close",0);

	vcr.wfly=make_window(base,4*12*DCURXs,4,5+DCURXs+5,(DCURYs+6)-4,1);
	vcr.view=make_plain_window(base,10,100,wid,hgt,2);
	ani_gc=XCreateGC(display,vcr.view,valuemask,&values);
	vcr.hgt=hgt;
	vcr.wid=wid;
	ani_pixmap=  XCreatePixmap(display,RootWindow(display,screen),vcr.wid,vcr.hgt,
							   DefaultDepth(display,screen));
	if(ani_pixmap==0) {
		err_msg("Failed to get the required pixmap");
		XFlush(display);
		waitasec(ClickTime);
		XDestroySubwindows(display,base);
		XDestroyWindow(display,base);
		vcr.iexist=0;
		return;
	}
	vcr.iexist=1;

	XSetFunction(display,ani_gc,GXcopy);
	XSetForeground(display,ani_gc,WhitePixel(display,screen));
	XFillRectangle(display,ani_pixmap,ani_gc,0,0,vcr.wid,vcr.hgt);
	XSetForeground(display,ani_gc,BlackPixel(display,screen));
	XSetFont(display,ani_gc,romfonts[0]->fid);
	tst_pix_draw();
	get_global_colormap(ani_pixmap);
	mpeg.flag=0;
	mpeg.filflag=0;
	strcpy(mpeg.root,"frame");
	mpeg.filter[0]=0;
	mpeg.skip=1;
	vcr.pos=0;
	if(use_ani_file) {
		get_ani_file(vcr.file);
	}
}


static void ani_border(Window w, int i) {
	if( w==vcr.wgrab  ||
		w==vcr.wgo	  ||
		w==vcr.wreset ||
		w==vcr.wpause ||
		w==vcr.wfast  ||
		w==vcr.wfile  ||
		w==vcr.wslow  ||
		w==vcr.wmpeg  ||
		w==vcr.wup    ||
		w==vcr.wdn	  ||
		w==vcr.wskip  ||
		w==vcr.kill) {
		XSetWindowBorderWidth(display,w,i);
	}
}


static void destroy_vcr(void) {
	vcr.iexist=0;
	XDestroySubwindows(display,vcr.base);
	XDestroyWindow(display,vcr.base);
}

/*************************  NEW ANIMaTION STUFF ***********************/
static void ani_motion_stuff(Window w,int x,int y) {
	if(w==vcr.view) {
		update_ani_motion_stuff(x,y);
	}
}


static double get_current_time(void) {
	double t1;
	struct timeval tim;
	gettimeofday(&tim,NULL);
	t1=tim.tv_sec+(tim.tv_usec/1000000.0);
	return t1;
}


static void update_ani_motion_stuff(int x,int y) {
	double dt;
	ami.t2=ami.t1;
	ami.t1=get_current_time();
	ami.ox=ami.x;
	ami.oy=ami.y;
	ani_ij_to_xy(x,y,&ami.x,&ami.y);
	dt=ami.t1-ami.t2;
	if(dt==0.0)dt=10000000000;
	ami.vx=(ami.x-ami.ox)/dt;
	ami.vy=(ami.y-ami.oy)/dt;

	dt=ami.tstart-ami.t2;
	if(dt==0.0)dt=100000000000;
	ami.vax=(ami.x0-ami.x)/dt;
	ami.vay=(ami.y0-ami.y)/dt;
	set_val("mouse_x",ami.x);
	set_val("mouse_y",ami.y);
	set_val("mouse_vx",ami.vx);
	set_val("mouse_vy",ami.vy);
	do_grab_tasks(1);
	fix_only();
	ani_frame(0);
}


/*************************** End motion & speed stuff   ****************/
static void ani_buttonx(XEvent ev,int flag) {
	Window w=ev.xbutton.window;
	/*   ADDED FOR THE GRAB FEATURE IN ANIMATOR  This is BUTTON PRESS */
	if((w==vcr.view) && (ani_grab_flag==1)) {
		if(flag==1) {
			ami.t1=get_current_time();
			ami.tstart=ami.t1;
			ani_ij_to_xy(ev.xbutton.x,ev.xbutton.y,&ami.x,&ami.y);
			ami.x0=ami.x;
			ami.y0=ami.y;
			who_was_grabbed=search_for_grab(ami.x,ami.y);
			if(who_was_grabbed<0) {
				printf("Nothing grabbed\n");
			}
		}
		if(flag==0) { /* This is BUTTON RELEASE  */
			/*  update_ani_motion_stuff(ev.xbutton.x,ev.xbutton.y); */
			if(who_was_grabbed<0) {
				/*  ani_grab_flag=0; */
				return;
			}
			do_grab_tasks(2);
			set_to_init_data();
			ani_grab_flag=0;
			redraw_params();
			if(run_now_grab()) {
				run_now();
				ani_grab_flag=0;
			}
		}
		return;
	}
	if(flag==0) {
		return;
	}
	/*   END OF ADDED STUFF  ************************/
	ani_button(w);
}


static void ani_button(Window w) {
	if((ani_grab_flag==1)) {
		return;
	}
	/* Grab button resets and shows first frame */
	if(w==vcr.wgrab) {
		if(n_ani_grab==0) {
			return;
		}
		if(vcr.ok) {
			vcr.pos=0;
			show_grab_points=1;
			/* ani_flip1(0); */
			ani_frame(1);
			ani_frame(0);
			ani_grab_flag=1;
		}
	}
	if(w==vcr.wmpeg) {
		ani_create_mpeg();
	}
	if(w==vcr.wgo){
		ani_flip();
	}
	if(w==vcr.wskip) {
		ani_newskip();
	}
	if(w==vcr.wup) {
		ani_flip1(1);
	}
	if(w==vcr.wdn) {
		ani_flip1(-1);
	}
	if(w==vcr.wfile) {
		get_ani_file(NULL);
	}
	if(w==vcr.wfly) {
		animation_on_the_fly=1-animation_on_the_fly;
		check_on_the_fly();
	}
	if(w==vcr.wreset) {
		vcr.pos=0;
		reset_comets();
		redraw_ani_slider();
		ani_flip1(0);
	}
	if(w==vcr.kill) {
		destroy_vcr();
	}
}


static void ani_create_mpeg(void) {
	static char *n[]={"PPM 0/1","Basename","AniGif(0/1)" };
	char values[3][MAX_LEN_SBOX];
	int status;

	mpeg.flag=0;
	sprintf(values[0],"%d",mpeg.flag);
	sprintf(values[1],"%s",mpeg.root);
	sprintf(values[2],"%d",mpeg.aviflag);
	status=do_string_box(3,3,1,"Frame saving",n,values,28);
	if(status!=FORGET_ALL) {
		mpeg.flag=atoi(values[0]);
		if(mpeg.flag>0) {
			mpeg.flag=1;
		}
		mpeg.aviflag=atoi(values[2]);
		sprintf(mpeg.root,"%s",values[1]);
		if(mpeg.aviflag==1) {
			mpeg.flag=0;
		}
	} else {
		mpeg.flag=0;
	}
	if(mpeg.flag==1) {
		ani_disk_warn();
	}
}


static void do_ani_slider_motion(Window w, int x) {
	int l=48*DCURXs,x0=x;
	int mr=my_browser.maxrow;
	int k;

	if(w!=vcr.slider || mr<2) {
		return;
	}
	if(x0>l-2) {
		x0=l-2;
	}
	vcr.slipos=x0;
	draw_ani_slider(w,x0);
	k=x0*mr/l;
	vcr.pos=0;
	ani_flip1(0);
	ani_flip1(k);
}


static void redraw_ani_slider(void) {
	int k=vcr.pos;
	int l=48*DCURXs;
	int xx;

	int mr=my_browser.maxrow;
	if(mr<2) {
		return;
	}
	xx=(k*l)/mr;
	draw_ani_slider(vcr.slider,xx);
}


static void draw_ani_slider(Window w,int x) {
	int hgt=DCURYs+4,l=48*DCURXs;
	int x0=x-2,i;

	if(x0<0) {
		x0=0;
	}
	if(x0>(l-4)) {
		x0=l-4;
	}
	XClearWindow(display,w);
	for(i=0;i<4;i++) {
		XDrawLine(display,w,small_gc,x0+i,0,x0+i,hgt);
	}
}


static void ani_resize(int x, int y) {
	int ww=x-(2*4);
	int hh=y-((2.5*(DCURYs+6))+5);
	if(ww==vcr.wid && hh==vcr.hgt) {
		return;
	}
	XFreePixmap(display,ani_pixmap);

	vcr.hgt=5*((y-((4.5*(DCURYs+6))+5))/5);
	vcr.wid=4*((x-(2*4))/4);

	/*This little safety check prevents a <X Error of failed request:  BadValue>
 from occuring if the user shrinks the window size smaller than the vcr.hgt | vcr.wid
 */
	if (vcr.hgt < 1) {
		vcr.hgt = 1;
	}
	if (vcr.wid < 1) {
		vcr.wid = 1;
	}

	XMoveResizeWindow(display,vcr.view,4,4.5*(DCURYs+6),vcr.wid,vcr.hgt);
	ani_pixmap=  XCreatePixmap(display,RootWindow(display,screen),vcr.wid,vcr.hgt,
							   DefaultDepth(display,screen));
	if(ani_pixmap==0) {
		err_msg("Failed to get the required pixmap");
		XFlush(display);
		XDestroySubwindows(display,vcr.base);
		XDestroyWindow(display,vcr.base);
		vcr.iexist=0;
		return;
	}

	XSetFunction(display,ani_gc,GXcopy);
	XSetForeground(display,ani_gc,WhitePixel(display,screen));
	XFillRectangle(display,ani_pixmap,ani_gc,0,0,vcr.wid,vcr.hgt);
	XSetForeground(display,ani_gc,BlackPixel(display,screen));
	tst_pix_draw();
}


static void ani_newskip(void) {
	char bob[20];
	Window w;
	int rev,status;

	XGetInputFocus(display,&w,&rev);
	sprintf(bob,"%d",vcr.inc);
	status=get_dialog("Frame skip","Increment:",bob,"Ok","Cancel",20);
	if(status!=FORGET_ALL) {
		vcr.inc=atoi(bob);
		if(vcr.inc<=0) {
			vcr.inc=1;
		}
	}
	XSetInputFocus(display,w,rev,CurrentTime);
}


static void check_on_the_fly(void) {
	XClearWindow(display,vcr.wfly);
	if(animation_on_the_fly) {
		XDrawString(display,vcr.wfly,small_gc,5,1.5*CURY_OFFs,"*",1);
	}
}


static void ani_frame(int task) {
	XSetForeground(display,ani_gc,WhitePixel(display,screen));
	XFillRectangle(display,ani_pixmap,ani_gc,0,0,vcr.wid,vcr.hgt);
	XSetForeground(display,ani_gc,BlackPixel(display,screen));
	if(task==1) {
		set_ani_perm();
		reset_comets();
		return;
	}
	/* now draw the stuff  */
	render_ani();

	/*  done drawing   */
	XCopyArea(display,ani_pixmap,vcr.view,ani_gc,0,0,vcr.wid,vcr.hgt,0,0);
	XFlush(display);
}


static void set_to_init_data(void) {
	int i;
	for(i=0;i<NODE;i++) {
		last_ic[i]=get_ivar(i+1);
	}
	for(i=NODE+FIX_VAR;i<NODE+FIX_VAR+NMarkov;i++) {
		last_ic[i-FIX_VAR]=get_ivar(i+1);
	}
	redraw_ics();
}


static void set_from_init_data(void) {
	double y[MAXODE];
	int i;

	for(i=0;i<NODE+NMarkov;i++) {
		y[i]=last_ic[i];
	}
	set_fix_rhs(T0,y);
}


static void ani_flip1(int n) {
	int row;
	float **ss;
	double y[MAXODE];
	double t;
	int i;
	if(n_anicom==0  ||  my_browser.maxrow<2) {
		return;
	}
	ss=my_browser.data;
	XSetForeground(display,ani_gc,WhitePixel(display,screen));
	XFillRectangle(display,ani_pixmap,ani_gc,0,0,vcr.wid,vcr.hgt);
	XSetForeground(display,ani_gc,BlackPixel(display,screen));
	if(vcr.pos==0) {
		set_ani_perm();
	}
	vcr.pos=vcr.pos+n;
	if(vcr.pos>=my_browser.maxrow) {
		vcr.pos=my_browser.maxrow-1;
	}
	if(vcr.pos<0) {
		vcr.pos=0;
	}
	row=vcr.pos;

	t=(double)ss[0][row];
	for(i=0;i<NODE+NMarkov;i++) {
		y[i]=(double)ss[i+1][row];
	}
	set_fix_rhs(t,y);

	/* now draw the stuff  */
	render_ani();

	/*  done drawing   */
	XCopyArea(display,ani_pixmap,vcr.view,ani_gc,0,0,vcr.wid,vcr.hgt,0,0);
	XFlush(display);
}


static void ani_flip(void) {
	double y[MAXODE];
	double t;
	char fname[DEFAULT_STRING_LENGTH];
	FILE *angiffile=NULL;
	float **ss;
	int i,row,done;
	int mpeg_frame=0,mpeg_write=0,count=0;
	XEvent ev;
	Window w;

	done=0;
	if(n_anicom==0  ||  my_browser.maxrow<2) {
		return;
	}
	ss=my_browser.data;
	set_ani_perm(); /* evaluate all permanent structures  */
	/* check avi_flags for initialization */
	if(mpeg.aviflag==1) {
		angiffile=fopen("anim.gif","wb");
		set_global_map(1);
	}
	count=0;
	while(!done) { /* Ignore all events except the button presses */
		if(XPending(display)>0) {
			XNextEvent(display,&ev);
			switch(ev.type) {
			case ButtonPress:
				w=ev.xbutton.window;
				if(w==vcr.wpause) {
					done=1;
					break;
				}
				if(w==vcr.wfast) {
					ani_speed=ani_speed-ani_speed_inc;
					if(ani_speed<0) {
						ani_speed=0;
					}
					break;
				}
				if(w==vcr.wslow) {
					ani_speed=ani_speed+ani_speed_inc;
					if(ani_speed>100) {
						ani_speed=100;
					}
					break;
				}
				break;
			}
		}
		/* Okay no events  so lets go! */
		/* first set all the variables */
		XSetForeground(display,ani_gc,WhitePixel(display,screen));
		XFillRectangle(display,ani_pixmap,ani_gc,0,0,vcr.wid,vcr.hgt);
		XSetForeground(display,ani_gc,BlackPixel(display,screen));
		row=vcr.pos;
		t=(double)ss[0][row];
		for(i=0;i<NODE+NMarkov;i++) {
			y[i]=(double)ss[i+1][row];
		}
		set_fix_rhs(t,y);

		/* now draw the stuff  */
		render_ani();

		/*  done drawing   */
		XCopyArea(display,ani_pixmap,vcr.view,ani_gc,0,0,vcr.wid,vcr.hgt,0,0);
		XFlush(display);
		waitasec(ani_speed);
		if(mpeg.aviflag==1 || mpeg.flag>0) {
			waitasec(5*ani_speed);
		}
		vcr.pos=vcr.pos+vcr.inc;
		if(vcr.pos>=my_browser.maxrow) {
			done=1;
			vcr.pos=0;
			reset_comets();
		}

		/* now check mpeg stuff */
		if(mpeg.flag>0 && ((mpeg_frame%mpeg.skip)==0)) {
			sprintf(fname,"%s_%d.ppm",mpeg.root,mpeg_write);
			mpeg_write++;
			writeframe(fname,ani_pixmap,vcr.wid,vcr.hgt);
		}
		mpeg_frame++;
		/* now check AVI stuff */

		if(mpeg.aviflag==1) {
			add_ani_gif(vcr.view,angiffile,count);
		}
		count++;
	}
	/* always stop mpeg writing */
	mpeg.flag=0;
	if(mpeg.aviflag==1) {
		end_ani_gif(angiffile);
		fclose(angiffile);
		set_global_map(0);
	}
}


static void ani_disk_warn(void) {
	unsigned int total=(my_browser.maxrow*vcr.wid*vcr.hgt*3)/(mpeg.skip*vcr.inc);
	char junk[DEFAULT_STRING_LENGTH];
	char ans;
	total=total/(1024*1024);
	if(total>10) {
		sprintf(junk," %d Mb disk space needed! Continue?",total);
		ans=(char)TwoChoice("YES","NO",junk,"yn");
		if(ans!='y') {
			mpeg.flag=0;
		}
	}
}


static int ani_new_file(char *filename) {
	FILE *fp;
	char bob[100];
	fp=fopen(filename,"r");
	if(fp==NULL) {
		err_msg("Couldn't open ani-file");
		return -1;
	}
	if(n_anicom>0) {
		free_ani();
	}
	if(load_ani_file(fp)==0) {
		sprintf(bob,"Bad ani-file at line %d",ani_line);
		err_msg(bob);
		return -1;
	}
	return 0;
}


static int load_ani_file(FILE *fp) {
	char old[DEFAULT_STRING_LENGTH],new_char[DEFAULT_STRING_LENGTH],big[DEFAULT_STRING_LENGTH];
	int notdone=1,jj1,jj2,jj;
	int ans=0,flag;
	ani_line=1;
	while(notdone) {
		read_ani_line(fp,old);
		search_array(old,new_char,&jj1,&jj2,&flag);
		for(jj=jj1;jj<=jj2;jj++) {
			subsk(new_char,big,jj,flag);
			ans=parse_ani_string(big,fp);
		}

		if(ans==0 || feof(fp)) {
			break;
		}
		if(ans<0) { /* error occurred !! */
			plintf(" error at line %d\n",ani_line);
			free_ani();
			return 0;
		}
		ani_line++;
	}
	return 1;
}


/*  This has changed to add the FILE fp to the arguments since the GRAB
*	command requires that you load in two additional lines
*/
static int parse_ani_string(char *s, FILE *fp) {
	char x1[DEFAULT_STRING_LENGTH],x2[DEFAULT_STRING_LENGTH],
		 x3[DEFAULT_STRING_LENGTH],x4[DEFAULT_STRING_LENGTH];
	char col[DEFAULT_STRING_LENGTH],thick[DEFAULT_STRING_LENGTH];
	char *ptr,*nxt;
	char *command, *toksave;
	int type=-1;
	int anss;

	x1[0]=0;
	x2[0]=0;
	x3[0]=0;
	x4[0]=0;
	col[0]=0;
	thick[0]=0;
	ptr=s;
	type=COMNT;
	command=strtok_r(ptr,"; ", &toksave);
	if (command == NULL) {
		return -1;
	}
	strupr(command);
	/************** GRAB STUFF *****************/
	if(strprefix("GR",command)) {
		type=GRAB;
	}
	if(strprefix("LI",command)) {
		type=LINE;
	}
	if(strprefix("RL",command)) {
		type=RLINE;
	}
	if(strprefix("RE",command)) {
		type=RECT;
	}
	if(strprefix("FR",command)) {
		type=FRECT;
	}
	if(strprefix("EL",command)) {
		type=ELLIP;
	}
	if(strprefix("FE",command)) {
		type=FELLIP;
	}
	if(strprefix("CI",command)) {
		type=CIRC;
	}
	if(strprefix("FC",command)) {
		type=FCIRC;
	}
	if(strprefix("VT",command)) {
		type=VTEXT;
	}
	if(strprefix("TE",command)) {
		type=TEXT;
	}
	if(strprefix("SE",command)){
		type=SETTEXT;
	}
	if(strprefix("TR",command)) {
		type=TRANSIENT;
	}
	if(strprefix("PE",command)) {
		type=PERMANENT;
	}
	if(strprefix("DI",command)) {
		type=DIMENSION;
	}
	if(strprefix("EN",command)) {
		type=END;
	}
	if(strprefix("DO",command)) {
		type=END;
	}
	if(strprefix("SP",command)) {
		type=SPEED;
	}
	if(strprefix("CO",command)) {
		type=COMET;
	}
	if(strprefix("XN",command)) {
		type=AXNULL;
	}
	if(strprefix("YN",command)) {
		type=AYNULL;
	}
	switch(type) {
	case GRAB:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x3,nxt);
		anss=add_grab_command(x1,x2,x3,fp);
		return(anss);
	case AXNULL:
	case AYNULL:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x3,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x4,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(col,nxt);
		nxt=strtok_r(NULL, "\n", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(thick,nxt);
		break;
	case LINE:
	case RECT:
	case ELLIP:
	case FELLIP:
	case FRECT:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x3,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x4,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if((nxt==NULL) || strlen(nxt)==0) {
			break;
		}
		strcpy(col,nxt);
		nxt=strtok_r(NULL, "\n", &toksave);
		if((nxt==NULL) || strlen(nxt)==0) {
			break;
		}
		strcpy(thick,nxt);
		break;
	case RLINE:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if((nxt==NULL) || strlen(nxt)==0) {
			break;
		}
		strcpy(col,nxt);
		nxt=strtok_r(NULL, "\n", &toksave);
		if((nxt==NULL) || strlen(nxt)==0) {
			break;
		}
		strcpy(thick,nxt);
		break;
	case COMET:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(thick,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x3,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if((nxt==NULL) || strlen(nxt)==0) {
			break;
		}
		strcpy(col,nxt);
		break;
	case CIRC:
	case FCIRC:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x3,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if((nxt==NULL) || strlen(nxt)==0) {
			break;
		}
		strcpy(col,nxt);
		break;
	case SETTEXT:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(col,nxt);
		break;
	case TEXT:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x4,nxt);
		break;
	case VTEXT:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x4,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x3,nxt);
		break;
	case SPEED:
		nxt=strtok_r(NULL, " \n", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		ani_speed=atoi(nxt);
		if(ani_speed<0) {
			ani_speed=0;
		}
		if(ani_speed>1000) {
			ani_speed=1000;
		}
		return 1;
	case DIMENSION:
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x1,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x2,nxt);
		nxt=strtok_r(NULL, ";", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x3,nxt);
		nxt=strtok_r(NULL, ";\n", &toksave);
		if(nxt==NULL) {
			return -1;
		}
		strcpy(x4,nxt);
		break;
	}
	if(type==END) {
		return 0;
	}
	if(type==TRANSIENT) {
		aniflag=TRANSIENT;
		return 1;
	}
	if(type==COMNT) {
		return 1;
	}
	if(type==PERMANENT) {
		aniflag=PERMANENT;
		return 1;
	}

	if(type==DIMENSION) {
		set_ani_dimension(x1,x2,x3,x4);
		return 1;
	}
	return(add_ani_com(type,x1,x2,x3,x4,col,thick));
}


static void set_ani_dimension(char *x1, char *y1, char *x2, char *y2) {
	double xx1,yy1,xx2,yy2;
	xx1=atof(x1);
	xx2=atof(x2);
	yy1=atof(y1);
	yy2=atof(y2);

	if((xx1<xx2) && (yy1<yy2)) {
		ani_xlo=xx1;
		ani_xhi=xx2;
		ani_ylo=yy1;
		ani_yhi=yy2;
	}
}


static int add_ani_com(int type, char *x1, char *y1, char *x2, char *y2,
					   char *col, char *thick) {
	int err=0;
	if(type==COMNT		||
	   type==DIMENSION	||
	   type==PERMANENT	||
	   type==TRANSIENT	||
	   type==END		||
	   type==SPEED) {
		return 1;
	}
	my_ani[n_anicom].type=type;
	my_ani[n_anicom].flag=aniflag;
	my_ani[n_anicom].x1=(int *)malloc(256*sizeof(int));
	my_ani[n_anicom].y1=(int *)malloc(256*sizeof(int));
	my_ani[n_anicom].x2=(int *)malloc(256*sizeof(int));
	my_ani[n_anicom].y2=(int *)malloc(256*sizeof(int));
	my_ani[n_anicom].col=(int *)malloc(256*sizeof(int));
	my_ani[n_anicom].who=(int *)malloc(256*sizeof(int));
	switch(type) {
	case AXNULL:
	case AYNULL:
		err=add_ani_null(&my_ani[n_anicom],x1,y1,x2,y2,col,thick);
		break;
	case COMET:
		err=add_ani_comet(&my_ani[n_anicom],x1,y1,x2,y2,col,thick);
		break;
	case LINE:
		err=add_ani_line(&my_ani[n_anicom],x1,y1,x2,y2,col,thick);
		break;
	case RLINE:
		err=add_ani_rline(&my_ani[n_anicom],x1,y1,col,thick);
		break;
	case RECT:
		err=add_ani_rect(&my_ani[n_anicom],x1,y1,x2,y2,col,thick);
		break;
	case FRECT:
		err=add_ani_frect(&my_ani[n_anicom],x1,y1,x2,y2,col,thick);
		break;
	case ELLIP:
		err=add_ani_ellip(&my_ani[n_anicom],x1,y1,x2,y2,col,thick);
		break;
	case FELLIP:
		err=add_ani_fellip(&my_ani[n_anicom],x1,y1,x2,y2,col,thick);
		break;
	case CIRC:
		err=add_ani_circle(&my_ani[n_anicom],x1,y1,x2,col,thick);
		break;
	case FCIRC:
		err=add_ani_circle(&my_ani[n_anicom],x1,y1,x2,col,thick);
		break;
	case TEXT:
		err=add_ani_text(&my_ani[n_anicom],x1,y1,y2);
		break;
	case VTEXT:
		err=add_ani_vtext(&my_ani[n_anicom],x1,y1,x2,y2);
		break;
	case SETTEXT:
		err=add_ani_settext(&my_ani[n_anicom],x1,y1,col);
		break;
	}
	if(err<0) {
		free_ani();
		return -1;
	}
	n_anicom++;
	return 1;
}


static void init_ani_stuff(void) {
	ani_text_size=1;
	ani_text_font=0;
	ani_text_color=0;
	ani_xlo=0.0;
	ani_ylo=0.0;
	ani_xhi=1.0;
	ani_yhi=1.0;
	aniflag=TRANSIENT;
	n_anicom=0;
	ani_lastx=0.0;
	ani_lasty=0.0;
	vcr.pos=0;
	ani_grab_flag=0; /*********** GRABBER *******************/
	n_ani_grab=0;
}


static void free_ani(void) {
	int i;
	for(i=0;i<n_anicom;i++) {
		free(my_ani[i].x1);
		free(my_ani[i].y1);
		free(my_ani[i].x2);
		free(my_ani[i].y2);
		free(my_ani[i].who);
		free(my_ani[i].col);
		if(my_ani[i].type==COMET) {
			free(my_ani[i].c.x);
			free(my_ani[i].c.y);
			free(my_ani[i].c.col);
		}
	}
	n_anicom=0;
	free_grabber();
	init_ani_stuff();
}


static int chk_ani_color(char *s, int *index) {
	int j;
	char *s2;

	*index=-1;
	de_space(s);
	strupr(s);
	if(strlen(s)==0) {
		*index=0;
		return 1;
	}
	if(s[0]=='$') {
		s2=&s[1];
		for(j=0;j<12;j++) {
			if(strcmp(s2,color_names[j])==0) {
				*index=colorline[j];
				return 1;
			}
		}
	}
	return 0;
}


static int add_ani_expr(char *x, int *c) {
	int i,n;
	int com[DEFAULT_STRING_LENGTH];
	int err;

	err=add_expr(x,com,&n);
	if(err==1) {
		return 1;
	}
	for(i=0;i<n;i++) {
		c[i]=com[i];
	}
	return 0;
}


/*  the commands  */
static int add_ani_rline(ANI_COM *a, char *x1, char *y1, char *col, char *thick) {
	int err,index;
	err=chk_ani_color(col,&index);
	if(err==1) {
		a->col[0]=index;
	} else {
		err=add_ani_expr(col,a->col);
		if(err==1) {
			return -1;
		}
	}
	a->zthick=atoi(thick);
	if(a->zthick <0) {
		a->zthick=0;
	}
	err=add_ani_expr(x1,a->x1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(y1,a->y1);
	if(err) {
		return -1;
	}
	return 0;
}


static void reset_comets(void) {
	int i;
	for(i=0;i<n_anicom;i++) {
		if(my_ani[i].type==COMET) {
			my_ani[i].c.i=0;
		}
	}
}


static void roll_comet(ANI_COM *a, int xn, int yn, int col) {
	int i;

	int n=a->c.n;
	int ii=a->c.i;
	if(ii<n) { /* not loaded yet */
		a->c.x[ii]=xn;
		a->c.y[ii]=yn;
		a->c.col[ii]=col;
		a->c.i=a->c.i+1;
		return;
	}
	/* its full so push down eliminating last */
	for(i=1;i<n;i++) {
		a->c.x[i-1]=a->c.x[i];
		a->c.y[i-1]=a->c.y[i];
		a->c.col[i-1]=a->c.col[i];
	}
	a->c.x[n-1]=xn;
	a->c.y[n-1]=yn;
	a->c.col[n-1]=col;
}


static int add_ani_comet(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						 char *col, char *thick) {
	int err,n,index;
	/* plintf("<<%s>>\n",col); */
	err=chk_ani_color(col,&index);
	if(err==1) {
		a->col[0]=-index;
	} else {
		err=add_ani_expr(col,a->col);
		if(err==1) {
			return -1;
		}
	}
	a->zthick=atoi(thick);
	n=atoi(x2);
	if(n<=0) {
		plintf("4th argument of comet must be positive integer!\n");
		return(-1);
	}
	err=add_ani_expr(x1,a->x1);
	if(err)
		return -1;

	err=add_ani_expr(y1,a->y1);
	if(err) {
		return -1;
	}
	a->c.n=n;
	a->c.x=(int *)malloc(n*sizeof(int));
	a->c.y=(int *)malloc(n*sizeof(int));
	a->c.col=(int *)malloc(n*sizeof(int));
	a->c.i=0;
	return 1;
}


static int add_ani_line(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						char *col, char *thick) {
	int err,index;
	err=chk_ani_color(col,&index);
	if(err==1) {
		a->col[0]=-index;
	} else {
		err=add_ani_expr(col,a->col);
		if(err==1) {
			return -1;
		}
	}
	a->zthick=atoi(thick);
	if(a->zthick <0) {
		a->zthick=0;
	}
	err=add_ani_expr(x1,a->x1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(y1,a->y1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(x2,a->x2);
	if(err) {
		return -1;
	}
	err=add_ani_expr(y2,a->y2);
	if(err) {
		return -1;
	}
	return 0;
}


static int add_ani_null(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						char *col, char *who) {
	int err,index;
	err=chk_ani_color(col,&index);
	if(err==1) {
		a->col[0]=-index;
	} else {
		err=add_ani_expr(col,a->col);
		if(err==1) {
			return -1;
		}
	}
	err=add_ani_expr(who,a->who);
	if(err) {
		return -1;
	}
	err=add_ani_expr(x1,a->x1);
	if(err){
		return -1;
	}
	err=add_ani_expr(y1,a->y1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(x2,a->x2);
	if(err) {
		return -1;
	}
	err=add_ani_expr(y2,a->y2);
	if(err) {
		return -1;
	}
	return 0;
}


static int add_ani_rect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						char *col, char *thick) {
	return(add_ani_line(a,x1,y1,x2,y2,col,thick));
}


static int add_ani_frect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						 char *col, char *thick) {
	return(add_ani_line(a,x1,y1,x2,y2,col,thick));
}


static int add_ani_ellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						 char *col, char *thick) {
	return(add_ani_line(a,x1,y1,x2,y2,col,thick));
}


static int add_ani_fellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2,
						  char *col, char *thick) {
	return(add_ani_line(a,x1,y1,x2,y2,col,thick));
}


static int add_ani_circle(ANI_COM *a, char *x1, char *y1, char *x2,
						  char *col, char *thick) {
	int err,index;
	err=chk_ani_color(col,&index);
	if(err==1) {
		a->col[0]=-index;
	} else {
		err=add_ani_expr(col,a->col);
		if(err==1) {
			return -1;
		}
	}
	a->zthick=atoi(thick);
	if(a->zthick <0) {
		a->zthick=0;
	}
	err=add_ani_expr(x1,a->x1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(y1,a->y1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(x2,a->x2);
	if(err) {
		return -1;
	}
	return 0;
}


static int add_ani_text(ANI_COM *a, char *x1, char *y1, char *y2) {
	int err;
	char *s;
	err=add_ani_expr(x1,a->x1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(y1,a->y1);
	if(err) {
		return -1;
	}
	s=(char *)(a->y2);
	strcpy(s,y2);
	return 0;
}


static int add_ani_vtext(ANI_COM *a, char *x1, char *y1, char *x2, char *y2) {
	int err;
	char *s;
	err=add_ani_expr(x1,a->x1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(y1,a->y1);
	if(err) {
		return -1;
	}
	err=add_ani_expr(x2,a->x2);
	if(err) {
		return -1;
	}
	s=(char *)(a->y2);
	strcpy(s,y2);
	return 0;
}


static int add_ani_settext(ANI_COM *a, char *x1, char *y1, char *col) {
	int size=atoi(x1);
	int font=0;
	int index=0,err;
	de_space(y1);
	if(y1[0]=='s' || y1[0]=='S') {
		font=1;
	}
	err=chk_ani_color(col,&index);
	if(err!=1) {
		index=0;
	}
	if(size<0) {
		size=0;
	}
	if(size>4) {
		size=4;
	}
	a->tsize=size;
	a->tfont=font;
	a->tcolor=index;
	return 0;
}


static void render_ani(void) {
	int i;
	int type,flag;
	redraw_ani_slider();
	for(i=0;i<n_anicom;i++) {
		type=my_ani[i].type;
		flag=my_ani[i].flag;
		if(type==LINE || type==RLINE || type==RECT || type==FRECT || type==CIRC ||
		   type==FCIRC || type==ELLIP || type==FELLIP || type==COMET || type==AXNULL ||
		   type==AYNULL) {
			eval_ani_color(i);
		}
		switch(type) {
		case AXNULL:
		case AYNULL:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_null(i,type-AXNULL);
			break;
		case SETTEXT:
			set_ani_font_stuff(my_ani[i].tsize,my_ani[i].tfont,my_ani[i].tcolor);
			break;
		case TEXT:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_text(i);
			break;
		case VTEXT:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_vtext(i);
			break;
		case LINE:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_line(i);
			break;
		case COMET:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_comet(i);
			break;
		case RLINE:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_rline(i);
			break;
		case RECT:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_rect(i);
			break;
		case FRECT:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_frect(i);
			break;
		case ELLIP:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_ellip(i);
			break;
		case FELLIP:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_fellip(i);
			break;
		case CIRC:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_circ(i);
			break;
		case FCIRC:
			if(flag==TRANSIENT) {
				eval_ani_com(i);
			}
			draw_ani_fcirc(i);
			break;
		}
	}
	if(show_grab_points==1) {
		draw_grab_points();
	}
}


static void set_ani_perm(void) {
	int i,type;
	set_from_init_data();
	for(i=0;i<n_anicom;i++) {
		type=my_ani[i].type;
		if(my_ani[i].flag==PERMANENT) {
			if(my_ani[i].type!=SETTEXT) {
				eval_ani_com(i);
			}
			if(type==LINE || type==RLINE || type==RECT || type==FRECT || type==CIRC ||
			   type==FCIRC || type==ELLIP || type==FELLIP) {
				eval_ani_color(i);
			}
		}
	}
}


static void eval_ani_color(int j) {
	double z;

	if(my_ani[j].col[0]>0) {
		z=evaluate(my_ani[j].col);
		if(z>1) {
			z=1.0;
		}
		if(z<0) {
			z=0.0;
		}
		my_ani[j].zcol=z;
	}
}


static void eval_ani_com(int j) {
	my_ani[j].zx1=evaluate(my_ani[j].x1);
	my_ani[j].zy1=evaluate(my_ani[j].y1);

	switch(my_ani[j].type) {
	case LINE:
	case RECT:
	case FRECT:
	case ELLIP:
	case FELLIP:
	case AXNULL:
	case AYNULL:
		my_ani[j].zx2=evaluate(my_ani[j].x2);
		my_ani[j].zy2=evaluate(my_ani[j].y2);
		break;
	case CIRC:
	case FCIRC:
		my_ani[j].zrad=evaluate(my_ani[j].x2);
		break;
	case VTEXT:
		my_ani[j].zval=evaluate(my_ani[j].x2);
		break;
	}
	if(my_ani[j].type==AXNULL || my_ani[j].type==AYNULL) {
		my_ani[j].zval=evaluate(my_ani[j].who);
	}
}


static void set_ani_thick(int t) {
	if(t<0) {
		t=0;
	}
	XSetLineAttributes(display,ani_gc,t,LineSolid,CapButt,JoinRound);
}

static void set_ani_font_stuff(int size, int font, int color) {
	if(color==0) {
		XSetForeground(display,ani_gc,BlackPixel(display,screen));
	} else {
		XSetForeground(display,ani_gc,ColorMap(color));
	}
	if(font==0) {
		XSetFont(display,ani_gc,romfonts[size]->fid);
	} else {
		XSetFont(display,ani_gc,symfonts[size]->fid);
	}
}


static void set_ani_col(int j) {
	int c=my_ani[j].col[0];
	int icol;

	if(c<=0) {
		icol=-c;
	} else {
		icol=(int)(color_total*my_ani[j].zcol)+FIRSTCOLOR;
	}
	/* plintf(" t=%d j=%d col=%d \n",vcr.pos,j,icol); */
	if(icol==0) {
		XSetForeground(display,ani_gc,BlackPixel(display,screen));
	} else {
		XSetForeground(display,ani_gc,ColorMap(icol));
	}
	LastAniColor=icol;
}


static void xset_ani_col(int icol) {
	if(icol==0) {
		XSetForeground(display,ani_gc,BlackPixel(display,screen));
	} else {
		XSetForeground(display,ani_gc,ColorMap(icol));
	}
}


/**************   DRAWING ROUTINES   *******************/
static void ani_rad2scale(double rx, double ry,int *ix, int *iy) {
	double dx=(double)vcr.wid/(ani_xhi-ani_xlo),
			dy=(double)vcr.hgt/(ani_yhi-ani_ylo);
	double r1=rx*dx,r2=ry*dy;
	*ix=(int)r1;
	*iy=(int)r2;
}


static void ani_radscale(double rad, int *ix, int *iy) {
	double dx=(double)vcr.wid/(ani_xhi-ani_xlo),
			dy=(double)vcr.hgt/(ani_yhi-ani_ylo);
	double r1=rad*dx,r2=rad*dy;
	*ix=(int)r1;
	*iy=(int)r2;
}


static void ani_ij_to_xy(int ix,int iy,double *x,double *y) {
	double dx=(ani_xhi-ani_xlo)/(double)vcr.wid;
	double dy=(ani_yhi-ani_ylo)/(double)vcr.hgt;
	*x=ani_xlo+(double)ix*dx;
	*y=ani_ylo+(double)(vcr.hgt-iy)*dy;


}


static void ani_xyscale(double x, double y, int *ix, int *iy) {
	double dx=(double)vcr.wid/(ani_xhi-ani_xlo),
			dy=(double)vcr.hgt/(ani_yhi-ani_ylo);
	double xx=(x-ani_xlo)*dx;
	double yy=vcr.hgt-dy*(y-ani_ylo);
	*ix=(int)xx;
	*iy=(int)yy;
	if(*ix<0) {
		*ix=0;
	}
	if(*ix>=vcr.wid) {
		*ix=vcr.wid-1;
	}
	if(*iy<0) {
		*iy=0;
	}
	if(*iy>=vcr.hgt) {
		*iy=vcr.hgt-1;
	}
}


static void draw_ani_comet(int j) {
	double x1=my_ani[j].zx1,y1=my_ani[j].zy1;
	int i1,j1,i2,j2;
	int k,nn,ir;
	set_ani_thick(my_ani[j].zthick);
	set_ani_col(j);
	ani_xyscale(x1,y1,&i1,&j1);
	/* we now have the latest x1,y1 */
	roll_comet(&my_ani[j],i1,j1,LastAniColor);
	/* now we draw this thing  */
	nn=my_ani[j].c.i;
	if(my_ani[j].zthick<0) {
		ir=-my_ani[j].zthick;
		for(k=0;k<nn;k++) {
			i1=my_ani[j].c.x[k];
			j1=my_ani[j].c.y[k];
			xset_ani_col(my_ani[j].c.col[k]);
			XFillArc(display,ani_pixmap,ani_gc,i1-ir,j1-ir,2*ir,2*ir,0,360*64);
		}
	} else {
		if(nn>2) {
			for(k=1;k<nn;k++) {
				i1=my_ani[j].c.x[k-1];
				j1=my_ani[j].c.y[k-1];
				i2=my_ani[j].c.x[k];
				j2=my_ani[j].c.y[k];
				xset_ani_col(my_ani[j].c.col[k]);
				XDrawLine(display,ani_pixmap,ani_gc,i1,j1,i2,j2);
			}
		}
	}
}


static void draw_ani_null(int j, int id) {
	double xl=my_ani[j].zx1,xh=my_ani[j].zx2,yl=my_ani[j].zy1,yh=my_ani[j].zy2;
	double z=my_ani[j].zval;
	float *v;
	int n,i,i4,who,i1,j1,i2,j2;
	float x1,y1,x2,y2,dx=xh-xl,dy=yh-yl;
	int err;
	if(dx==0.0 || dy==0.0) {
		return;
	}
	set_ani_col(j);
	who = (int)z;  /* the nullcline that you want  -1 is the default cline */
	err=get_nullcline_floats(&v,&n,who,id);
	if(err==1) {
		return;
	}
	for(i=0;i<n;i++) {
		i4=4*i;
		x1=(v[i4]-xl)/dx;
		y1=(v[i4+1]-yl)/dy;
		x2=(v[i4+2]-xl)/dx;
		y2=(v[i4+3]-yl)/dy;
		ani_xyscale(x1,y1,&i1,&j1);
		ani_xyscale(x2,y2,&i2,&j2);
		XDrawLine(display,ani_pixmap,ani_gc,i1,j1,i2,j2);
	}
}


static void draw_ani_line(int j) {
	double x1=my_ani[j].zx1,x2=my_ani[j].zx2,y1=my_ani[j].zy1,y2=my_ani[j].zy2;
	int i1,j1,i2,j2;

	set_ani_thick(my_ani[j].zthick);
	set_ani_col(j);
	ani_xyscale(x1,y1,&i1,&j1);
	ani_xyscale(x2,y2,&i2,&j2);
	XDrawLine(display,ani_pixmap,ani_gc,i1,j1,i2,j2);
	ani_lastx=x2;
	ani_lasty=y2;
}


static void draw_ani_rline(int j) {
	double x1=ani_lastx+my_ani[j].zx1,y1=ani_lasty+my_ani[j].zy1;
	int i1,j1,i2,j2;

	set_ani_thick(my_ani[j].zthick);
	set_ani_col(j);
	ani_xyscale(ani_lastx,ani_lasty,&i1,&j1);
	ani_xyscale(x1,y1,&i2,&j2);
	XDrawLine(display,ani_pixmap,ani_gc,i1,j1,i2,j2);
	ani_lastx=x1;
	ani_lasty=y1;
}


static void draw_ani_circ(int j) {
	double x1=my_ani[j].zx1,y1=my_ani[j].zy1,rad=my_ani[j].zrad;
	int i1,j1,i2,j2,ir;

	set_ani_col(j);
	set_ani_thick(my_ani[j].zthick);
	ani_xyscale(x1,y1,&i1,&j1);
	ani_radscale(rad,&i2,&j2);
	ir=(i2+j2)/2;
	XDrawArc(display,ani_pixmap,ani_gc,i1-ir,j1-ir,2*ir,2*ir,0,360*64);
}


static void draw_ani_fcirc(int j) {
	double x1=my_ani[j].zx1,y1=my_ani[j].zy1,rad=my_ani[j].zrad;
	int i1,j1,i2,j2,ir;

	set_ani_col(j);
	set_ani_thick(my_ani[j].zthick);
	ani_xyscale(x1,y1,&i1,&j1);
	ani_radscale(rad,&i2,&j2);
	ir=(i2+j2)/2;
	XFillArc(display,ani_pixmap,ani_gc,i1-ir,j1-ir,2*ir,2*ir,0,360*64);
}


static void draw_ani_rect(int j) {
	double x1=my_ani[j].zx1,x2=my_ani[j].zx2,y1=my_ani[j].zy1,y2=my_ani[j].zy2;
	int i1,j1,i2,j2;
	int h,w;
	set_ani_thick(my_ani[j].zthick);
	set_ani_col(j);
	ani_xyscale(x1,y1,&i1,&j1);
	ani_xyscale(x2,y2,&i2,&j2);
	h=abs(j2-j1);
	w=abs(i2-i1);
	if(i1>i2) {
		i1=i2;
	}
	if(j1>j2) {
		j1=j2;
	}
	XDrawRectangle(display,ani_pixmap,ani_gc,i1,j1,w,h);
}


static void draw_ani_frect(int j) {
	double x1=my_ani[j].zx1,x2=my_ani[j].zx2,y1=my_ani[j].zy1,y2=my_ani[j].zy2;
	int i1,j1,i2,j2;
	int h,w;
	set_ani_thick(my_ani[j].zthick);
	set_ani_col(j);
	ani_xyscale(x1,y1,&i1,&j1);
	ani_xyscale(x2,y2,&i2,&j2);

	h=abs(j2-j1);
	w=abs(i2-i1);
	if(i1>i2) {
		i1=i2;
	}
	if(j1>j2) {
		j1=j2;
	}
	XFillRectangle(display,ani_pixmap,ani_gc,i1,j1,w,h);
}


static void draw_ani_ellip(int j) {
	double x1=my_ani[j].zx1,x2=my_ani[j].zx2,y1=my_ani[j].zy1,y2=my_ani[j].zy2;
	int i1,j1,i2,j2;

	set_ani_thick(my_ani[j].zthick);
	set_ani_col(j);
	ani_xyscale(x1,y1,&i1,&j1);
	ani_rad2scale(x2,y2,&i2,&j2);
	XDrawArc(display,ani_pixmap,ani_gc,i1-i2,j1-j2,2*i2,2*j2,0,360*64);
}


static void draw_ani_fellip(int j) {
	double x1=my_ani[j].zx1,x2=my_ani[j].zx2,y1=my_ani[j].zy1,y2=my_ani[j].zy2;
	int i1,j1,i2,j2;

	set_ani_thick(my_ani[j].zthick);
	set_ani_col(j);
	ani_xyscale(x1,y1,&i1,&j1);
	ani_rad2scale(x2,y2,&i2,&j2);
	XFillArc(display,ani_pixmap,ani_gc,i1-i2,j1-j2,2*i2,2*j2,0,360*64);
}


static void draw_ani_text(int j) {
	int n;
	char *s;
	double x1=my_ani[j].zx1,y1=my_ani[j].zy1;
	int i1,j1;
	ani_xyscale(x1,y1,&i1,&j1);
	s=(char *)my_ani[j].y2;
	n=strlen(s);
	XDrawString(display,ani_pixmap,ani_gc,i1,j1,s,n);
}


static void draw_ani_vtext(int j) {
	char s2[DEFAULT_STRING_LENGTH];
	int n;
	char *s;
	double x1=my_ani[j].zx1,y1=my_ani[j].zy1;
	int i1,j1;
	s=(char *)my_ani[j].y2;
	sprintf(s2,"%s%g",s,my_ani[j].zval);
	n=strlen(s2);
	ani_xyscale(x1,y1,&i1,&j1);
	XDrawString(display,ani_pixmap,ani_gc,i1,j1,s2,n);
}


static void tst_pix_draw(void) {
	int i;
	XSetForeground(display,ani_gc,BlackPixel(display,screen));
	XDrawLine(display,ani_pixmap,ani_gc,0,2,vcr.wid,2);
	for(i=1;i<11;i++) {
		XSetForeground(display,ani_gc,ColorMap(colorline[i]));
		XDrawLine(display,ani_pixmap,ani_gc,0,2+i,vcr.wid,2+i);
	}
	for(i=0;i<=color_total;i++) {
		XSetForeground(display,ani_gc,ColorMap(i+FIRSTCOLOR));
		XDrawLine(display,ani_pixmap,ani_gc,0,14+i,vcr.wid,14+i);
	}
	XSetForeground(display,ani_gc,BlackPixel(display,screen));
	XDrawString(display,ani_pixmap,ani_gc,10,vcr.hgt-(DCURYs+6),"THIS SPACE FOR RENT",20);
}


static void read_ani_line(FILE *fp, char *s) {
	char temp[DEFAULT_STRING_LENGTH];
	int i,n,ok,ihat=0;
	/*int nn; Not used anywhere?*/
	s[0]=0;
	ok=1;
	while(ok) {
		ok=0;
		if(fgets(temp,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Error reading file %s", fp);
		}
		n=strlen(temp);
		for(i=n-1;i>=0;i--) {
			if(temp[i]=='\\') {
				ok=1;
				ihat=i;
			}
		}
		if(ok==1) {
			temp[ihat]=0;
		}
		strcat(s,temp);
	}
	n=strlen(s);
	if(s[n-1]=='\n') {
		s[n-1]=' ';
	}
	s[n]=' ';
	s[n+1]=0;
}


/*************************  GRABBER CODE *****************************/
static int add_grab_command(char *xs,char *ys,char *ts, FILE *fp) {
	char start[DEFAULT_STRING_LENGTH],end[DEFAULT_STRING_LENGTH];
	int com[DEFAULT_STRING_LENGTH];
	int nc,j,k,ans;
	double z;
	read_ani_line(fp,start);
	read_ani_line(fp,end);

	if(n_ani_grab>=MAX_ANI_GRAB) {
		plintf("Too many grabbables! \n");
		return(-1);
	}
	j=n_ani_grab;
	z=atof(ts);
	if(z<=0.0) {
		z=.02;
	}
	ani_grab[j].tol=z;
	if(add_expr(xs,com,&nc)) {
		plintf("Bad grab x %s \n",xs);
		return(-1);
	}
	ani_grab[j].x= (int *)malloc(sizeof(int)*(nc+1));
	for(k=0;k<=nc;k++) {
		ani_grab[j].x[k]=com[k];
	}
	if(add_expr(ys,com,&nc)) {
		plintf("Bad grab y %s \n",ys);
		return(-1);
	}
	ani_grab[j].y= (int *)malloc(sizeof(int)*(nc+1));
	for(k=0;k<=nc;k++) {
		ani_grab[j].y[k]=com[k];
	}
	ans=ani_grab_tasks(start,j,1);
	if(ans<0) {
		return(-1);
	}
	if(ani_grab_tasks(end,j,2)==(-1)) {
		return(-1);
	}
	n_ani_grab++;
	return(1);
}


static int ani_grab_tasks(char *line, int igrab,int which) {
	int i,k;
	int n=strlen(line);
	char form[DEFAULT_STRING_LENGTH],c;
	char rhs[DEFAULT_STRING_LENGTH],lhs[20];
	k=0;
	for(i=0;i<n;i++) {
		c=line[i];
		if(c=='{' || c==' ') {
			continue;
		}
		if(c==';' || c=='}') {
			form[k]=0;
			strcpy(rhs,form);
			if(add_grab_task(lhs,rhs,igrab,which)<0) {
				return(-1);
			}
			k=0;
			continue;
		}
		if(c=='=') {
			form[k]=0;
			strcpy(lhs,form);
			k=0;
			continue;
		}
		form[k]=c;
		k++;
	}
	return(1);

}


static int run_now_grab(void) {
	if(who_was_grabbed<0) {
		return(0);
	}
	return(ani_grab[who_was_grabbed].end.runnow);
}


static int search_for_grab(double x,double y) {
	int i;
	double d,u,v;
	double dmin=100000000;
	int imin=-1;
	for(i=0;i<n_ani_grab;i++) 	{
		u=ani_grab[i].zx;
		v=ani_grab[i].zy;
		d=sqrt((x-u)*(x-u)+(y-v)*(y-v));
		if((d<dmin)  &&  (d<ani_grab[i].tol)) {
			dmin=d;
			imin=i;
		}
	}
	return(imin);
}


/* which=1 for start, 2 for end */
static void do_grab_tasks(int which) {
	int i=who_was_grabbed;
	int j,n;
	double z;
	if(i<0) {
		return; /*  no legal grab point */
	}
	if(which==1) {
		n=ani_grab[i].start.n;
		for(j=0;j<n;j++) {
			z=evaluate(ani_grab[i].start.comrhs[j]);
			set_val(ani_grab[i].start.lhsname[j],z);
		}
		return;
	}
	if(which==2) {
		n=ani_grab[i].end.n;
		for(j=0;j<n;j++) {
			z=evaluate(ani_grab[i].end.comrhs[j]);
			set_val(ani_grab[i].end.lhsname[j],z);
		}
		return;
	}
}


static int add_grab_task(char *lhs,char *rhs, int igrab,int which) {
	int com[DEFAULT_STRING_LENGTH];
	int i,nc,k;
	int rn;
	if(which==1) {
		i=ani_grab[igrab].start.n;
		if(i>=MAX_GEVENTS) {
			return(-1); /* too many events */
		}
		strcpy(ani_grab[igrab].start.lhsname[i],lhs);
		if(add_expr(rhs,com,&nc)) {
			plintf("Bad right-hand side for grab event %s\n",rhs);
			return(-1);
		}
		ani_grab[igrab].start.comrhs[i] = (int *)malloc(sizeof(int)*(nc+1));
		for(k=0;k<=nc;k++) {
			ani_grab[igrab].start.comrhs[i][k]=com[k];
		}
		ani_grab[igrab].start.n=ani_grab[igrab].start.n+1;
		return(1);
	}
	if(which==2) {
		if(strncmp("runnow",lhs,6)==0) {
			rn=atoi(rhs);
			ani_grab[igrab].end.runnow=rn;
			return(1);
		}
		i=ani_grab[igrab].end.n;
		if(i>=MAX_GEVENTS) {
			return(-1); /* too many events */
		}
		strcpy(ani_grab[igrab].end.lhsname[i],lhs);
		if(add_expr(rhs,com,&nc)) {
			plintf("Bad right-hand side for grab event %s\n",rhs);
			plintf("should return -1\n");
			return(-1);
		}
		ani_grab[igrab].end.comrhs[i] = (int *)malloc(sizeof(int)*(nc+1));
		for(k=0;k<=nc;k++) {
			ani_grab[igrab].end.comrhs[i][k]=com[k];
		}
		ani_grab[igrab].end.n=ani_grab[igrab].end.n+1;
		return(1);
	}
	return(-1);
}


/* Draw little black x's where the grab points are */
static void draw_grab_points(void) {
	double xc,yc;
	double x1,y1,x2,y2,z;
	int i1,j1,i2,j2,ic,jc;
	int i;
	XSetForeground(display,ani_gc,BlackPixel(display,screen));
	for(i=0;i<n_ani_grab;i++) {
		xc=evaluate(ani_grab[i].x);
		yc=evaluate(ani_grab[i].y);
		ani_grab[i].zx=xc;
		ani_grab[i].zy=yc;
		z=ani_grab[i].tol;
		x1=xc+z;
		x2=xc-z;
		y1=yc+z;
		y2=yc-z;
		ani_xyscale(xc,yc,&ic,&jc);
		ani_xyscale(x1,y1,&i1,&j1);
		ani_xyscale(x2,y2,&i2,&j2);
		XDrawLine(display,ani_pixmap,ani_gc,i1,j1,i2,j2);
		XDrawLine(display,ani_pixmap,ani_gc,i1,j2,i2,j1);
	}
	show_grab_points=0;
}


static void free_grabber(void) {
	int i,j,m;
	for(i=0;i<n_ani_grab;i++) {
		free(ani_grab[i].x);
		free(ani_grab[i].y);
		m=ani_grab[i].start.n;
		for(j=0;j<m;j++) {
			free(ani_grab[i].start.comrhs[j]);
		}
		m=ani_grab[i].end.n;
		for(j=0;j<m;j++) {
			free(ani_grab[i].end.comrhs[j]);
		}
		ani_grab[i].start.n=0;
		ani_grab[i].end.n=0;
	}
}
