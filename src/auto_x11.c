#include "auto_x11.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "auto_nox.h"
#include "browse.h"
#include "color.h"
#include "derived.h"
#include "diagram.h"
#include "init_conds.h"
#include "integrate.h"
#include "form_ode.h"
#include "ggets.h"
#include "main.h"
#include "many_pops.h"
#include "menus.h"
#include "mykeydef.h"
#include "numerics.h"
#include "parserslow.h"
#include "pop_list.h"
#include "rubber.h"
#include "util/timeutil.h"
#include "xAuto.h"
#include "xpplim.h"

#include "bitmap/auto.bitmap"

/* --- Macros --- */
#define STD_WID 460	  /* golden mean  */
#define STD_HGT 284

#define SBW XSetWindowBorderWidth(display,w,1)


/* --- Types --- */
typedef struct  {
	Window canvas, axes,numerics,grab,next,run,clear,redraw,base,per;
	Window info,param,file,abort,stab,hint,kill;
} AUTOWIN;


/* --- Forward declarations --- */
static void a_msg(int i, int v);
static void auto_kill(void);
static void auto_update_view(float xlo,float xhi, float ylo, float yhi);
static void clear_msg(void);
static void find_point(int ibr,int pt);
static Window lil_button(Window root, int x, int y, char *name);
static void MarkAuto(int x,int y);
static int query_special(char* title,char *nsymb);
static void RedrawMark(void);
static void traverse_out(DIAGRAM *d, int *ix, int *iy, int dodraw);
static void XORCross(int x, int y);


/* --- Data --- */
int AutoRedrawFlag=1;

/* stuff for marking a branch  */
int mark_flag=0;
int mark_ibrs,mark_ibre;
int mark_ipts,mark_ipte;
int mark_ixs,mark_ixe,mark_iys,mark_iye;

static int STD_HGT_var =0;
static int STD_WID_var =0;
static int Auto_extra_wid,Auto_extra_hgt;
static int Auto_x0,Auto_y0;

static AUTOWIN AutoW;
static DIAGRAM *CUR_DIAGRAM;


/* --- Functions --- */
void ALINE(int a, int b, int c, int d) {
	XDrawLine(display,AutoW.canvas,small_gc,(a),(b),(c),(d));
}


void auto_button(XEvent ev) {
	Window w=ev.xbutton.window;
	if(Auto.exist==0) {
		return;
	}
	if(w==AutoW.axes) {
		SBW;
		auto_plot_par();
		return;
	}
	if(w==AutoW.numerics) {
		SBW;
		auto_num_par();
		return;
	}
	if(w==AutoW.grab) {
		SBW;
		auto_grab();
		return;
	}
	if(w==AutoW.run) {
		SBW;
		auto_run();
		return;
	}
	if(w==AutoW.redraw) {
		SBW;
		redraw_diagram();
		return;
	}
	if(w==AutoW.clear) {
		SBW;
		draw_bif_axes();
		return;
	}
	if(w==AutoW.per) {
		SBW;
		auto_per_par();
		return;
	}
	if(w==AutoW.param) {
		SBW;
		auto_params();
		return;
	}
	if(w==AutoW.kill) {
		SBW;
		auto_kill();
		return;
	}
	if(w==AutoW.file) {
		SBW;
		auto_file();
		return;
	}
}


/*  Auto event handlers   */
void auto_enter(Window w, int v) {
	if(Auto.exist==0) {
		return;
	}

	if(w==AutoW.param) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(0,v);
		return;
	}
	if(w==AutoW.axes) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(1,v);
		return;
	}
	if(w==AutoW.numerics) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(2,v);
		return;
	}
	if(w==AutoW.run) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(3,v);
		return;
	}
	if(w==AutoW.grab) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(4,v);
		return;
	}
	if(w==AutoW.per) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(5,v);
		return;
	}
	if(w==AutoW.clear) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(6,v);
		return;
	}
	if(w==AutoW.redraw) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(7,v);
		return;
	}
	if(w==AutoW.file) {
		XSetWindowBorderWidth(display,w,v);
		a_msg(8,v);
		return;
	}
	if(w==AutoW.kill) {
		XSetWindowBorderWidth(display,w,v);
		return;
	}
}


void  auto_get_info( int *n, char *pname) {
	int i1,i2,ibr;
	DIAGRAM *d,*dnew;

	if(mark_flag==2) {
		i1=abs(mark_ipts);
		ibr=mark_ibrs;
		i2=abs(mark_ipte);
		*n=abs(i2-i1);
		d=bifd;
		while(1) {
			if(d->ibr==ibr && ((d->ntot==i1)||(d->ntot==(-i1)))) {
				strcpy(pname,upar_names[AutoPar[d->icp1]]);
				break;
			}
			dnew=d->next;
			if(dnew==NULL) {
				break;
			}
			d=dnew;
		}
	}
}


void auto_keypress(XEvent ev, int *used) {
	Window w=ev.xkey.window;
	char ks;
	Window w2;
	int rev;

	*used=0;
	if(Auto.exist==0)
		return;
	XGetInputFocus(display,&w2,&rev);

	if(w==AutoW.base||w==AutoW.canvas||w2==AutoW.base) {
		*used=1;
		ks=(char)get_key_press(&ev);
		/* XLookupString(&ev,buf,maxlen,&ks,&comp); */

		if(ks=='a'||ks=='A') {
			auto_plot_par();
			return;
		}
		if(ks=='n'||ks=='N') {
			auto_num_par();
			return;
		}
		if(ks=='G'||ks=='g') {
			auto_grab();
			return;
		}
		if(ks=='R'||ks=='r') {
			auto_run();
			return;
		}
		if(ks=='D'||ks=='d') {
			redraw_diagram();
			return;
		}
		if(ks=='C'||ks=='c') {
			draw_bif_axes();
			return;
		}
		if(ks=='U'||ks=='u') {
			auto_per_par();
			return;
		}
		if(ks=='P'||ks=='p') {
			auto_params();
			return;
		}
		if(ks=='F'||ks=='f') {
			auto_file();
			return;
		}
		if(ks==ESC) {
			XSetInputFocus(display,command_pop, RevertToParent,CurrentTime);
			return;
		}
	}
}


void auto_motion(XEvent ev) {
	int i=ev.xmotion.x;
	int j=ev.xmotion.y;
	double x,y;
	Window w=ev.xmotion.window;
	if(Auto.exist==0) {
		return;
	}
	if(w==AutoW.canvas) {
		x=Auto.xmin+(double)(i-Auto.x0)*(Auto.xmax-Auto.xmin)/(double)Auto.wid;
		y=Auto.ymin+(double)(Auto.y0-j+Auto.hgt)*(Auto.ymax-Auto.ymin)/(double)Auto.hgt;
		sprintf(Auto.hinttxt,"x=%g,y=%g",x,y);
		storeautopoint(x,y);
		display_auto(AutoW.hint);
	}
}


int auto_pop_up_list(char *title, char **list, char *key, int n, int max, int def,
					 int x, int y, char **hints, char *httxt) {
	Window temp=AutoW.base;
	return pop_up_list(&temp,title,list,key,n,max,def,x,y,hints,AutoW.hint,httxt);
}


int auto_rubber(int *i1, int *j1, int *i2, int *j2, int flag) {
	return(rubber(i1,j1,i2,j2,AutoW.canvas,flag));
}


void auto_scroll_window(void) {
	XEvent ev;
	int i=0,j=0;
	int i0=0,j0=0;
	int state=0;
	float xlo=Auto.xmin;
	float ylo=Auto.ymin;
	float xhi=Auto.xmax;
	float yhi=Auto.ymax;
	float dx=0.0,dy=0.0;
	int alldone=0;
	XSelectInput(display,AutoW.canvas,
				 KeyPressMask|ButtonPressMask|ButtonReleaseMask|
				 PointerMotionMask|ButtonMotionMask|ExposureMask);
	while(!alldone) {
		XNextEvent(display,&ev);
		switch(ev.type) {
		case KeyPress:
			alldone=1;
			break;
		case Expose:
			do_expose(ev);
			break;
		case ButtonPress:
			if(state==0) {
				i0=ev.xkey.x;
				j0=ev.xkey.y;
				state=1;
			}
			break;
		case MotionNotify:
			if(state==1) {
				i0=ev.xmotion.x;
				j0=ev.xmotion.y;
				dx=0.0;
				dy=0.0;
				auto_update_view(xlo+dx,xhi+dx,ylo+dy,yhi+dy);
				state=2;
				break;
			}
			if(state==2) {
				i=ev.xmotion.x;
				j=ev.xmotion.y;
				dx=(float)(i0-i)*(Auto.xmax-Auto.xmin)/(float)Auto.wid;
				dy=(float)(j-j0)*(Auto.ymax-Auto.ymin)/(float)Auto.hgt;
				auto_update_view(xlo+dx,xhi+dx,ylo+dy,yhi+dy);
			}
			break;
		case ButtonRelease:
			state=0;
			xlo=xlo+dx;
			xhi=xhi+dx;
			ylo=ylo+dy;
			yhi=yhi+dy;
			break;
		}
	}
}


void auto_set_mark(int i) {
	int pt,ibr;
	if(mark_flag==2) {
		ibr=mark_ibrs;
		if(abs(mark_ipts)<abs(mark_ipte)) {
			pt=abs(mark_ipts)+i;
		} else {
			pt=abs(mark_ipte)+i;
			find_point(ibr,pt);
		}
	}
}


void auto_stab_line(int x,int y,int xp, int yp) {
	XDrawLine(display,AutoW.stab,small_gc,x,y,xp,yp);
}


void autobw(void) {
	XSetBackground(display,small_gc,MyBackColor);
	XSetForeground(display,small_gc,MyForeColor);
}


void autocol(int col) {
	set_scolor(col);
}


void ATEXT(int a, int b, char *c) {
	XDrawString(display,AutoW.canvas,small_gc,(a),(b),(c),strlen(c));
}


int byeauto_(int *iflag) {
	XEvent event;
	Window w;
	char ch;
	if(Auto.exist==0) {
		return(1);
	}
	*iflag=0;
	while(XPending(display)>0) {
		XNextEvent(display,&event);
		switch(event.type) {
		case Expose:
			do_expose(event);
			break;
		case ButtonPress:
			w=event.xbutton.window;
			if(w==AutoW.abort) {
				SBW;
				*iflag=1;
				return(1);
			}
			break;
		case KeyPress:
			ch=get_key_press(&event);
			if(ch==ESC) {
				*iflag=1;
				return(0);
			}
			break;
		}
	}
	return(0);
}


void Circle(int x, int y, int r) {
	XDrawArc(display,AutoW.canvas,small_gc,x-r,y-r,r<<1,r<<1,0,360*64);
}


void clear_auto_info(void) {
	XClearWindow(display,AutoW.info);
}


void clear_auto_plot(void) {
	XClearWindow(display,AutoW.canvas);
	redraw_auto_menus();
}


void clr_stab(void) {
	int r=Auto.st_wid/4;
	XClearWindow(display,AutoW.stab);
	XDrawArc(display,AutoW.stab,small_gc,r,r,2*r,2*r,0,360*64);
}


void display_auto(Window w) {

	int ix,iy;
	if(Auto.exist==0) {
		return;
	}
	if(w==AutoW.canvas) {
		if(AutoRedrawFlag==1) {
			redraw_diagram();
		}
	}
	if(w==AutoW.stab) {
		XFlush(display);
		int r=Auto.st_wid/4;
		XDrawArc(display,AutoW.stab,small_gc,r,r,2*r,2*r,0,360*64);
		if(CUR_DIAGRAM != NULL) {
			traverse_out(CUR_DIAGRAM,&ix,&iy,1);/*clr_stab();*/
		}
		XFlush(display);
	}
	if(w==AutoW.axes) {
		xds("Axes");
	}
	if(w==AutoW.numerics) {
		xds("Numerics");
	}
	if(w==AutoW.grab) {
		xds("Grab");
	}
	if(w==AutoW.run) {
		xds("Run");
	}
	if(w==AutoW.redraw) {
		xds("reDraw");
	}
	if(w==AutoW.clear) {
		xds("Clear");
	}
	if(w==AutoW.per) {
		xds("Usr period");
	}
	if(w==AutoW.kill) {
		xds("Close");
	}
	if(w==AutoW.param) {
		xds("Parameter");
	}
	if(w==AutoW.file) {
		xds("File");
	}
	if(w==AutoW.abort) {
		xds("ABORT");
	}
	if(w==AutoW.hint) {
		XClearWindow(display,w);
		XDrawString(display,w,gc,8,CURY_OFF,Auto.hinttxt,strlen(Auto.hinttxt));
		return;
	}
}


void DLINE(double a, double b, double c, double d) {
	ALINE(IXVal(a),IYVal(b),IXVal(c),IYVal(d));
}


void do_auto_range(void) {
	double t=TEND;

	if(mark_flag==2) {
		do_auto_range_go();
	}
	TEND=t;
}


void draw_auto_info(char *bob,int x,int y) {
	XDrawString(display,AutoW.info,small_gc,x,y,bob,strlen(bob));
}


void FillCircle(int x, int y, int r) {
	int  r2 = (int) (r / 1.41421356 + 0.5);
	int wh = 2 * r2;
	XFillArc(display, AutoW.canvas, small_gc, x - r2, y - r2, wh, wh, 0, 360*64);
}


void LineWidth(int wid) {
	XSetLineAttributes(display,small_gc, wid, LineSolid, CapButt, JoinRound);
}/* this makes the auto window  */
void make_auto(char *wname, char *iname) {
	int x,y,wid,hgt,addwid=16*DCURX,addhgt=3.0*DCURY,hinthgt=DCURY+6;
	Window base;
	int dely=DCURY+5;
	STD_HGT_var =20*DCURY;
	/*STD_WID_var =1.62*STD_HGT_var;*/
	STD_WID_var = 67*DCURX;
	int ymargin=4*DCURYs,xmargin=12*DCURXs;
	XTextProperty winname,iconname;
	XSizeHints size_hints;
	Auto_extra_wid=10+addwid;
	Auto_extra_hgt=addhgt+2*DCURY+hinthgt;
	wid=10+addwid+STD_WID_var+xmargin;
	hgt=addhgt+2*DCURY+STD_HGT_var+ymargin+hinthgt;
	x=addwid+5;
	y=DCURY;
	Auto_x0=x;
	Auto_y0=y;
	base=make_plain_window(RootWindow(display,screen),0,0,wid,hgt,4);
	XSetWindowBackground(display,base,MyMainWinColor);
	AutoW.base=base;
	strcpy(Auto.hinttxt,"hint");
	XSelectInput(display,base,ExposureMask|KeyPressMask|ButtonPressMask|
				 StructureNotifyMask);
	XStringListToTextProperty(&wname,1,&winname);
	XStringListToTextProperty(&iname,1,&iconname);

	size_hints.flags=PPosition|PSize|PMinSize;
	size_hints.x=0;
	size_hints.y=0;
	size_hints.min_width=wid;
	size_hints.min_height=hgt;

	XClassHint class_hints;
	class_hints.res_name="";
	class_hints.res_class="";

	XSetWMProperties(display,base,&winname,&iconname,NULL,0,
					 &size_hints,NULL,&class_hints);
	make_icon((char*)auto_bits,auto_width,auto_height,base);
	AutoW.canvas=make_plain_window(base,x,y,STD_WID_var+xmargin,STD_HGT_var+ymargin,1);
	XSetWindowBackground(display,AutoW.canvas,MyDrawWinColor);
	XSelectInput(display,AutoW.canvas,MYMASK);
	x=DCURX;
	y=DCURY+STD_HGT_var+ymargin-8*DCURX;
	AutoW.stab=make_plain_window(base,x,y,12*DCURX,12*DCURX,2);
	Auto.st_wid=12*DCURX;
	x=DCURX+2;
	y=2*DCURY;
	Auto.hgt=STD_HGT_var;
	Auto.wid=STD_WID_var;
	Auto.x0=10*DCURXs;
	Auto.y0=2*DCURYs;
	AutoW.kill=lil_button(base,2,2,"Close");
	AutoW.param=lil_button(base,x,y,"Parameter");
	y+=dely;
	AutoW.axes=lil_button(base,x,y,"Axes");
	y+=dely;
	AutoW.numerics=lil_button(base,x,y,"Numerics");
	y+=dely;
	AutoW.run=lil_button(base,x,y,"Run");
	y+=dely;
	AutoW.grab=lil_button(base,x,y,"Grab");
	y+=dely;
	AutoW.per=lil_button(base,x,y,"Usr Function");
	y+=dely;
	AutoW.clear=lil_button(base,x,y,"Clear");
	y+=dely;
	AutoW.redraw=lil_button(base,x,y,"reDraw");
	y+=dely;
	AutoW.file=lil_button(base,x,y,"File");

	y+=3*dely;
	AutoW.abort=lil_button(base,x,y,"ABORT");

	y=DCURY+STD_HGT_var+ymargin+5;
	x=addwid+5;
	AutoW.info=make_plain_window(base,x,y,STD_WID_var+xmargin,addhgt,2);
	AutoW.hint=make_plain_window(base,x,y+addhgt+6,STD_WID_var+xmargin,DCURY+2,2);
	draw_bif_axes();
}


void redraw_auto_menus(void) {
	display_auto(AutoW.axes);
	display_auto(AutoW.numerics);
	display_auto(AutoW.grab);
	display_auto(AutoW.run);
	display_auto(AutoW.redraw);
	display_auto(AutoW.clear);
	display_auto(AutoW.per);
	display_auto(AutoW.param);
	display_auto(AutoW.kill);
	display_auto(AutoW.file);
	display_auto(AutoW.abort);
}


void refreshdisplay(void) {
	XFlush(display);
}


void resize_auto_window(XEvent ev) {

	int wid,hgt,addhgt=3.5*DCURY;

	STD_HGT_var =20*DCURY;
	/*STD_WID_var =1.62*STD_HGT_var;*/
	STD_WID_var = 50*DCURX;
	int ymargin=4*DCURYs,xmargin=12*DCURXs;
	if(ev.xconfigure.window==AutoW.base) {
		wid=ev.xconfigure.width-Auto_extra_wid;
		hgt=ev.xconfigure.height-Auto_extra_hgt;

		addhgt=3.0*DCURY;

		XResizeWindow(display,AutoW.canvas,wid,hgt);
		Window root;
		int xloc;
		int yloc;
		unsigned int cwid;
		unsigned int chgt;
		unsigned int cbwid;
		unsigned int cdepth;

		XGetGeometry(display,AutoW.canvas,&root,&xloc,&yloc,&cwid,&chgt,&cbwid,&cdepth);
		Auto.hgt=chgt-ymargin;
		Auto.wid=cwid-xmargin;

		XMoveResizeWindow(display,AutoW.info,xloc,yloc+chgt+4,wid,addhgt);
		XMoveResizeWindow(display,AutoW.hint,xloc,yloc+chgt+addhgt+10,wid,DCURY+2);
		int ix,iy;

		if(NBifs<2) {
			return;
		}
		traverse_out(CUR_DIAGRAM,&ix,&iy,1);
	}
}



void traverse_diagram(void) {
	DIAGRAM *d,*dnew,*dold;
	int done=0;
	int ix,iy,i;
	int lalo;
	XEvent ev;
	int kp;
	mark_flag=0;
	if(NBifs<2) {
		return;
	}
	d=bifd;
	DONT_XORCross=0;
	traverse_out(d,&ix,&iy,1);

	while(done==0) {
		XNextEvent(display,&ev);
		if(ev.type==ButtonPress) {
			int xm=ev.xmotion.x;
			int ym=ev.xmotion.y;
			Window w=ev.xmotion.window;

			if(w==AutoW.canvas)	{
				clear_msg();
				/* GO HOME */
				XORCross(ix,iy);
				DONT_XORCross = 1;
				while (1) {
					dnew=d->prev;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					/*bifd = dnew;*/
					d=dnew;
				}
				d=dnew;
				CUR_DIAGRAM=d;
				traverse_out(d,&ix,&iy,0);
				/* END GO HOME */
				/* GO END */
				int mindex=0;
				double dist;
				double ndist = Auto.wid*Auto.hgt;
				XORCross(ix,iy);
				lalo=load_all_labeled_orbits;
				load_all_labeled_orbits=0;
				while (1) {
					dist = sqrt(((double)(xm-ix))*((double)(xm-ix)) + ((double)(ym-iy))*((double)(ym-iy)));
					if(dist<ndist)	{
						ndist = dist;
						mindex=d->index;
					}
					dnew=d->next;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					d=dnew;
					traverse_out(d,&ix,&iy,0);/*Need this each time to update the distance calc*/
				}
				d=dnew;
				CUR_DIAGRAM=d;
				load_all_labeled_orbits=lalo;
				traverse_out(d,&ix,&iy,0);
				/* END GO END */
				/* GO HOME */
				XORCross(ix,iy);
				while (1) {
					if(d->index == mindex) {
						dnew=d;
						break;
					}
					dnew=d->prev;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					/*bifd = dnew;*/
					d=dnew;
				}
				d=dnew;
				CUR_DIAGRAM=d;
				DONT_XORCross = 0;
				traverse_out(d,&ix,&iy,1);
				/* END GO HOME */
			}
		} else if(ev.type==KeyPress) {
			clear_msg();
			kp=get_key_press(&ev);
			char symb[3],nsymb[3];

			int found=0;

			switch(kp) {
			case RIGHT:
				dnew=d->next;
				if(dnew==NULL) {
					dnew=bifd;
				}
				XORCross(ix,iy);
				d=dnew;
				CUR_DIAGRAM=dnew;
				traverse_out(d,&ix,&iy,1);
				break;

			case LEFT:
				dnew=d->prev;
				if(dnew==NULL) {
					dnew=bifd;
				}
				XORCross(ix,iy);
				d=dnew;
				CUR_DIAGRAM=dnew;
				traverse_out(d,&ix,&iy,1);
				break;
			case UP:
				if(!query_special("Next...",nsymb)) {
					break;
				}
				XORCross(ix,iy);
				found=0;
				dold=d;
				while(1) {
					dnew=d->next;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					get_bif_sym(symb,dnew->itp);
					if(strcmp(symb,nsymb)==0) {
						d=dnew;
						found=1;
						break;
					}
					d=dnew;
					/*if(d->lab==0)break;*/
				}
				if(found) {
					d=dnew;
				} else {
					snprintf(Auto.hinttxt,255,"  Higher %s not found",nsymb);
					display_auto(AutoW.hint);
					d=dold;
				}
				CUR_DIAGRAM=d;
				traverse_out(d,&ix,&iy,1);
				break;
			case DOWN:
				if(!query_special("Previous...",nsymb)) {
					break;
				}
				XORCross(ix,iy);
				found=0;
				dold=d;
				while(1) {
					dnew=d->prev;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					get_bif_sym(symb,dnew->itp);
					if(strcmp(symb,nsymb)==0) {
						d=dnew;
						found=1;
						break;
					}
					d=dnew;
				}
				if(found) {
					d=dnew;
				} else {
					snprintf(Auto.hinttxt,255,"  Lower %s not found",nsymb);
					display_auto(AutoW.hint);
					d=dold;
				}
				CUR_DIAGRAM=d;
				traverse_out(d,&ix,&iy,1);
				break;
			case TAB:
				XORCross(ix,iy);
				while(1) {
					dnew=d->next;
					if(dnew==NULL) {
						dnew=bifd;
						break;
					} /*TAB wraps*/
					d=dnew;
					if(d->lab!=0) {
						break;
					}
				}
				d=dnew;
				CUR_DIAGRAM=d;
				traverse_out(d,&ix,&iy,1);
				break;
				/* New code */
			case 's': /* mark the start of a branch */
				if(mark_flag==0) {
					MarkAuto(ix,iy);
					mark_ibrs=d->ibr;
					mark_ipts=d->ntot;
					mark_flag=1;
					mark_ixs=ix;
					mark_iys=iy;
				}
				break;
			case 'e': /* mark end of branch */
				if(mark_flag==1) {
					MarkAuto(ix,iy);
					mark_ibre=d->ibr;
					mark_ipte=d->ntot;
					mark_flag=2;
					mark_ixe=ix;
					mark_iye=iy;
				}
				break;
			case END:/*All the way to end*/
				XORCross(ix,iy);
				while (1) {
					dnew=d->next;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					/*bifd = dnew;*/
					d=dnew;
				}
				d=dnew;
				CUR_DIAGRAM=d;
				traverse_out(d,&ix,&iy,1);
				break;
			case HOME:/*All the way to beginning*/
				XORCross(ix,iy);
				while (1) {
					dnew=d->prev;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					/*bifd = dnew;*/
					d=dnew;
				}
				d=dnew;
				CUR_DIAGRAM=d;
				traverse_out(d,&ix,&iy,1);
				break;
			case PGUP: /*Same as TAB except we don't wrap*/
				XORCross(ix,iy);
				while(1) {
					dnew=d->next;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					d=dnew;
					if(d->lab!=0) {
						break;
					}
				}
				d=dnew;
				CUR_DIAGRAM=d;
				traverse_out(d,&ix,&iy,1);
				break;
			case PGDOWN: /*REVERSE TAB*/
				XORCross(ix,iy);
				while(1) {
					dnew=d->prev;
					if(dnew==NULL) {
						dnew=d;
						break;
					}
					d=dnew;
					if(d->lab!=0) {
						break;
					}
				}
				d=dnew;
				CUR_DIAGRAM=d;
				traverse_out(d,&ix,&iy,1);
				break;

			case FINE:
				done=1;
				XORCross(ix,iy);
				/*Cross should be erased now that we have made our selection.*/
				/*Seems XORing it with new draw can tend to bring it back randomly
	depending on the order of window expose events.  Best not
	to do the XORCross function at all.*/
				DONT_XORCross = 1;
				redraw_diagram();
				RedrawMark();
				break;
			case ESC:
				done=-1;
				break;
			}
		}

	}
	/*XORCross(ix,iy); */
	/* check mark_flag branch similarity */
	if(mark_flag==2) {
		if(mark_ibrs!=mark_ibre) {
			mark_flag=0;
		}
	}
	if(done==1) {
		grabpt.ibr=d->ibr;
		grabpt.lab=d->lab;
		for(i=0;i<8;i++)
			grabpt.par[i]=d->par[i];

		grabpt.icp1=d->icp1;
		grabpt.icp2=d->icp2;
		grabpt.per=d->per;
		grabpt.torper=d->torper;
		for(i=0;i<NODE;i++) {
			grabpt.uhi[i]=d->uhi[i];
			grabpt.ulo[i]=d->ulo[i];
			grabpt.u0[i]=d->u0[i];
			grabpt.ubar[i]=d->ubar[i];
			set_ivar(i+1,grabpt.u0[i]);
		}
		get_ic(0,grabpt.u0);
		grabpt.flag=1;
		grabpt.itp=d->itp;
		grabpt.ntot=d->ntot;
		grabpt.nfpar=d->nfpar;
		grabpt.index=d->index;
		for(i=0;i<NAutoPar;i++)
			constants[Auto_index_to_array[i]]=grabpt.par[i];
	}
	evaluate_derived();
	redraw_params();
	redraw_ics();
}


/* --- Static functions --- */
static void a_msg(int i, int v) {
	if(v==0||TipsFlag==0) {
		return;
	}
	snprintf(Auto.hinttxt,DEFAULT_STRING_LENGTH-1,"%s",auto_hint[i]);
	display_auto(AutoW.hint);
}


static void auto_kill(void) {
	Auto.exist=0;
	waitasec(ClickTime);
	XDestroySubwindows(display,AutoW.base);
	XDestroyWindow(display,AutoW.base);
}


static void auto_update_view(float xlo,float xhi, float ylo, float yhi) {
	Auto.xmin=xlo;
	Auto.ymin=ylo;
	Auto.xmax=xhi;
	Auto.ymax=yhi;
	redraw_diagram();
}


static void clear_msg(void) {
	Auto.hinttxt[0]='\0';
	display_auto(AutoW.hint);
}


static void find_point(int ibr, int pt) {
	int i;
	DIAGRAM *d,*dnew;
	if(NBifs<2) {
		return;
	}
	d=bifd;
	while(1) {
		if(d->ibr==ibr && ((d->ntot==pt)||(d->ntot==(-pt)))) {
			/* need to look at both signs to ignore stability */
			/* now we use this info to set parameters and init data */
			for(i=0;i<NODE;i++)
				set_ivar(i+1,d->u0[i]);

			get_ic(0,d->u0);
			for(i=0;i<NAutoPar;i++)
				constants[Auto_index_to_array[i]]=d->par[i];

			evaluate_derived();
			redraw_params();
			redraw_ics();
			if((d->per)>0) {
				set_total(d->per);
			}
			break;
		}
		dnew=d->next;
		if(dnew==NULL) {
			break;
		}
		d=dnew;
	}
}


static Window lil_button(Window root, int x, int y, char *name) {
	Window win;
	/*int width=strlen(name)*DCURX+5;
  */
	int width=12*DCURX;
	win=make_window(root,x,y,width,DCURY+1,1);
	XSelectInput(display,win,MYMASK);
	return(win);
}


static void MarkAuto(int x, int y) {
	LineWidth(2);
	ALINE(x-8,y-8,x+8,y+8);
	ALINE(x+8,y-8,x-8,y+8);
	LineWidth(1);
}


static int query_special(char* title,char *nsymb) {
	int status=1;
	static char *m[]={"BP","EP","HB","LP","MX","PD","TR","UZ"};
	static  char key[]="behlmptu";
	int ch=(char)auto_pop_up_list(title,m,key,8,11,1,10,10,
								  aspecial_hint,Auto.hinttxt);
	if(ch=='b') {
		sprintf(nsymb,"BP");
	} else if(ch=='e') {
		sprintf(nsymb,"EP");
	} else if(ch=='h') {
		sprintf(nsymb,"HB");
	} else if(ch=='l') {
		sprintf(nsymb,"LP");
	} else if(ch=='m') {
		sprintf(nsymb,"MX");
	} else if(ch=='p') {
		sprintf(nsymb,"PD");
	} else if(ch=='t') {
		sprintf(nsymb,"TR");
	} else if(ch=='u') {
		sprintf(nsymb,"UZ");
	} else {
		status=0;
		sprintf(nsymb,"  ");
	}
	redraw_auto_menus();
	return(status);
}


static void RedrawMark(void) {
	if(mark_flag==2) {
		MarkAuto(mark_ixs,mark_iys);
		MarkAuto(mark_ixe,mark_iye);
	}
}


static void traverse_out(DIAGRAM *d, int *ix, int *iy, int dodraw) {
	double norm,per,*par,par1,par2=0,*evr,*evi;
	int pt,itp,ibr,lab,icp1,icp2,flag2;
	double x,y1,y2;
	char symb[3];
	if(d==NULL) {
		return;
	}
	norm=d->norm;
	par=d->par;

	per=d->per;
	lab=d->lab;
	itp=d->itp;
	ibr=d->ibr;
	icp1=d->icp1;
	icp2=d->icp2;
	flag2=d->flag2;
	pt=d->ntot;

	evr=d->evr;
	evi=d->evi;

	get_bif_sym(symb,itp);
	par1=par[icp1];
	if(icp2<NAutoPar) {
		par2=par[icp2];
	}
	auto_xy_plot(&x,&y1,&y2,par1,par2,per,d->uhi,d->ulo,d->ubar,norm);

	*ix=IXVal(x);
	*iy=IYVal(y1);
	if(dodraw==1) {
		XORCross(*ix,*iy);
		plot_stab(evr,evi,NODE);
		new_info(ibr,pt,symb,lab,par,norm,d->u0[Auto.var],per,flag2,icp1,icp2);
	}
	if(lab>0 && load_all_labeled_orbits>0) {
		load_auto_orbitx(ibr,1,lab,per);
	}
}


static void XORCross(int x, int y) {
	if(DONT_XORCross) {
		return;
	}

	if(xorfix) {
		XSetForeground(display,small_gc,MyDrawWinColor);
		XSetBackground(display,small_gc,MyForeColor);
	}

	XSetFunction(display,small_gc,GXxor);
	LineWidth(2);
	ALINE(x-8,y,x+8,y);
	ALINE(x,y+8,x,y-8);
	XSetFunction(display,small_gc,GXcopy);
	LineWidth(1);
	if(xorfix) {
		XSetForeground(display,small_gc,MyForeColor);
		XSetBackground(display,small_gc,MyDrawWinColor);
	}
	XFlush(display);
}
