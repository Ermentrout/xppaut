#ifndef XPPAUT_BROWSE_H
#define XPPAUT_BROWSE_H

#include <stdio.h>
#include <X11/Xlib.h>

#include "xpplim.h"

/* --- Macros --- */
#define xds(a) { XDrawString(display,w,small_gc,5,CURY_OFFs,a,strlen(a));\
	return;}

#define BMAXCOL 20

#define MYMASK	   (ButtonPressMask 	|\
	KeyPressMask		|\
	ExposureMask		|\
	StructureNotifyMask	|\
	LeaveWindowMask		|\
	EnterWindowMask)

#define SIMPMASK   (ButtonPressMask		|\
	KeyPressMask		|\
	ExposureMask		|\
	StructureNotifyMask)

/* --- Types --- */
typedef struct {
	Window base,upper;
	Window find,up,down,pgup,pgdn,home,end,left,right;
	Window first,last,restore,write,get,close;
	Window load,repl,unrepl,table,addcol,delcol;
	Window main;
	Window label[BMAXCOL];
	Window time;
	Window hint;
	char hinttxt[DEFAULT_STRING_LENGTH];
	int dataflag,xflag;
	int col0,row0,ncol,nrow;
	int maxrow,maxcol;
	float **data;
	int istart,iend;
} BROWSER;

/* --- Data --- */
extern BROWSER my_browser;

/* --- Functions --- */
Window br_button(Window root, int row, int col, char *name, int iflag);
void data_get_mybrowser(int row);
void expose_my_browser(XEvent ev);
void enter_my_browser(XEvent ev, int yn);
void find_variable(char *s, int *col);
float **get_browser_data(void);
float *get_data_col(int c);
void get_data_xyz(float *x, float *y, float *z, int i1, int i2, int i3, int off);
int get_maxrow_browser(void);
void init_browser(void);
void make_new_browser(void);
void my_browse_button(XEvent ev);
void my_browse_keypress(XEvent ev, int *used);
void open_write_file(FILE **fp, char *fil, int *ok);
void refresh_browser(int length);
void reset_browser(void);
void resize_my_browser(Window win);
void set_browser_data(float **data, int col0);
void wipe_rep(void);
void write_mybrowser_data(FILE *fp);

#endif /* XPPAUT_BROWSE_H */











