/*
	Copyright (C) 2002-2015  Bard Ermentrout & Daniel Dougherty

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

	The author can be contacted at
	bard@pitt.edu
*/
#include "main.h"

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <X11/cursorfont.h>
#include <X11/Xatom.h>
#include <X11/Xlib.h>
#include <X11/Xos.h>
#include <X11/Xproto.h>
#include <X11/Xutil.h>

#include "adj2.h"
#include "aniparse.h"
#include "arrayplot.h"
#include "auto_nox.h"
#include "auto_x11.h"
#include "axes2.h"
#include "browse.h"
#include "calc.h"
#include "color.h"
#include "comline.h"
#include "dae_fun.h"
#include "do_fit.h"
#include "edit_rhs.h"
#include "eig_list.h"
#include "extra.h"
#include "form_ode.h"
#include "ggets.h"
#include "graf_par.h"
#include "graphics.h"
#include "help_defs.h"
#include "init_conds.h"
#include "integrate.h"
#include "lunch-new.h"
#include "many_pops.h"
#include "menu.h"
#include "menudrive.h"
#include "myfonts.h"
#include "my_rhs.h"
#include "nullcline.h"
#include "numerics.h"
#include "pop_list.h"
#include "read_dir.h"
#include "simplenet.h"
#include "storage.h"
#include "struct.h"
#include "txtread.h"
#include "userbut.h"

#include "bitmap/pp.bitmap"

/* --- Macros --- */
#define lowbit(x) ((x) & (~(x) + 1))
#define TOPBUTTONCOLOR 27
#define TOO_SMALL 0
#define BIG_ENOUGH 1
#define cstringmaj MYSTR1
#define cstringmin MYSTR2

/* --- Forward declarations --- */
static XKeyEvent createKeyEvent(Window win, Window winRoot, int press, int keycode, int modifiers);
static void do_events(unsigned int min_wid, unsigned int min_hgt);
static Window init_win(unsigned int bw, char *icon_name, char *win_name, int x, int y, unsigned int min_wid, unsigned int min_hgt, int argc, char **argv);
static void getGC(GC *gc);
static int getxcolors(XWindowAttributes *win_info, XColor **colors);
static void init_X(void);
static void load_fonts(void);
static void make_pops(void);
static void make_top_buttons(void);
static int script_make(char *s,int *k);
static void set_big_font(void);
static void test_color_info(void);
static void top_button_cross(Window w, int b);
static void top_button_events(XEvent report);
static void top_button_press(Window w);
static void xpp_events(XEvent report, int min_wid, int min_hgt);

/* --- Data --- */
static int ALREADY_SWAPPED=0;
static unsigned int Black,White;
static Window TopButton[6];

char anifile[XPP_MAX_NAME];
char batchout[XPP_MAX_NAME];
char big_font_name[MAXVEC],small_font_name[MAXVEC];
char PlotFormat[MAXVEC];
char UserBlack[8];
char UserWhite[8];
char UserMainWinColor[8];
char UserDrawWinColor[8];
char UserBGBitmap[XPP_MAX_NAME];
char UserOUTFILE[XPP_MAX_NAME];

int (*rhs)();

int allwinvis=0;
int DCURYb,DCURXb,CURY_OFFb;
int DCURYs,DCURXs,CURY_OFFs;
int DCURY,DCURX,CURY_OFF;
int DisplayHeight,DisplayWidth;
/*Set this to 1 if you want the tutorial to come up at start-up as default behavior */
int DoTutorial=0;
int PaperWhite=-1;
int SCALEX,SCALEY;
int screen;
int tfBell;
int TrueColorFlag;
int use_ani_file=0;
int UserGradients=-1;
int UserMinWidth=0,UserMinHeight=0;
int Xup,TipsFlag=1;
int XPPBatch=0,batch_range=0;
int XPPVERBOSE=1;

float xppvermaj,xppvermin;

unsigned int MyBackColor,MyForeColor,MyMainWinColor,MyDrawWinColor;
unsigned int GrFore,GrBack;

Atom deleteWindowAtom=0;
Display *display;
FILE *logfile;
GC gc, gc_graph,small_gc, font_gc;

Window draw_win;
Window main_win;
Window command_pop,info_pop;
XFontStruct *big_font,*small_font;


/* --- Functions --- */
void bye_bye(void) {
	int i;
	yes_reset_auto();
	XUnloadFont(display,big_font->fid);
	XUnloadFont(display,small_font->fid);
	for(i=0;i<5;i++) {
		if(avsymfonts[i]) {
			XUnloadFont(display,symfonts[i]->fid);
		}
		if(avromfonts[i]) {
			XUnloadFont(display,romfonts[i]->fid);
		}
	}
	XFreeGC(display,gc);
	XCloseDisplay(display);
	exit(1);
}


void clr_scrn(void) {
	blank_screen(draw_win);
	restore_off();
	do_axes();
}


void commander(int ch) {
	switch(help_menu) {
	case MAIN_HELP:
	{
		switch(ch) {
		case 'i':
			/*  initial data  */
			ini_data_menu();
			break;
		case 'c':
			/* continue */
			cont_integ();
			break;
		case 'n':
			/*nullclines */
			new_clines();
			break;
		case 'd':
			/*dir fields */
			direct_field();
			break;
		case 'w':
			/* window */
			window_zoom();
			break;
		case 'a':
			/*phase-space */
			do_torus();
			break;
		case 'k':
			/*kinescope */
			do_movie();
			break;
		case 'g' :
			add_a_curve();
			break;
		case 'u':
			help_num();
			break;
		case 'f':
			/* files */
			help_file();
			break;
		case 'p':
			/*parameters */
			/* change_par(-1); */
			new_param();
			break;
		case 'e':
			/*erase */
			clear_screens();
			break;
		case 'h':
		case 'm':
			do_windows();
			/*half windows */
			break;
		case 't':
			/*text */
			do_gr_objs();
			break;
		case 's':
			/*sing pts */
			find_equilibrium();
			break;
		case 'v':
			/*view_axes */
			change_view();
			break;
		case 'b':
			find_bvp();
			break;
		case 'x':
			/*x vs t */
			x_vs_t();
			break;
		case 'r':
			/*restore*/
			redraw_them_all();
			break;
		case '3': get_3d_par();
			break;
		case 'y':
			/* scripty();  */
			break;
		} /* End main switch  */
	} /* MAIN HELP ENDS  */
		break;
	case NUM_HELP:
	{
		get_num_par(ch);
	} /* end num case   */ break;
	case FILE_HELP:
	{
		switch(ch) {
		case 't':
			do_transpose();
			break;
		case 'g':
			get_intern_set();
			break;
		case 'i':
			TipsFlag=1-TipsFlag;
			break;
		case 'p':
			/* file stuff */
			/* do_info(stdout); */
			make_txtview();
			break;
		case 'w':
			/* write set */
			do_lunch(WRITEM);
			break;
		case 's':
			/* make eqn */
			file_inf();
			break;
		case 'a':
			/* AUTO !! */
#ifdef AUTO
			do_auto_win();
#endif
			break;
		case 'c':
			/* calculator */
			q_calc();
			break;
		case 'r':
			/* read set */
			do_lunch(READEM);
			break;
		case 'e':
			/* script */
			edit_menu();
			/*  Insert generic code here ...  */
			break;
		case 'b':
			tfBell=1-tfBell;
			break;
		case 'h':
			xpp_hlp();
			/*	   make_key_stroke(); */
			break;
		case 'q':
			if(yes_no_box()) {
				bye_bye();
			}
			/* quit */
			break;
			/* CLONE ! */
		case 'l':
			clone_ode();
			break;
		case 'x':
			edit_xpprc();
			break;
		case 'u':
			do_tutorial();
			break;
		} /* end file switch  */
		help();
	} /*  end file case   */ break;
	}  /* end help_menu switch  */
	/* redraw_menu(); */
}


void FixWindowSize(Window w, int width, int height, int flag) {
	XSizeHints size_hints;
	switch(flag) {
	case FIX_SIZE:
		size_hints.flags=PSize|PMinSize|PMaxSize;
		size_hints.width=width;
		size_hints.min_width=width;
		size_hints.max_width=width;
		size_hints.height=height;
		size_hints.min_height=height;
		size_hints.max_height=height;
		break;
	case FIX_MIN_SIZE:
		size_hints.flags=PMinSize;
		size_hints.min_width=width;
		size_hints.min_height=height;
		break;
	case FIX_MAX_SIZE:
		size_hints.flags=PMaxSize;
		size_hints.max_width=width;
		size_hints.max_height=height;
		break;
	}
	XSetWMProperties(display,w,NULL,NULL,NULL,0,&size_hints,NULL,NULL);
}


int main(int argc, char **argv) {
	char myfile[XPP_MAX_NAME];
	char pptitle[80];
	unsigned int min_wid=450,min_hgt=360;

	get_directory(myfile);

	SCALEX=640;
	SCALEY=480;

	Xup=0;
	sprintf(batchout,"output.dat");
	sprintf(PlotFormat,"ps");

	logfile=stdout;

	loadeqn_init_options();
	loadeqn_load_xpprc();
	do_comline(argc, argv);

	if(!XPPBatch) {
		init_X();
	}
	/* Load the ODE file (possibly by showing a file selector). */
	load_eqn();

	/* Let command line options override ODE file options. */
	do_comline(argc, argv);

	/* Set up subsystems based on options. */
	init_alloc_info();
	loadeqn_setup_all();
	set_init_guess();
	update_all_ffts();

#ifdef AUTO
	init_auto_win();
#endif

	if(disc(this_file)) {
		METHOD=METHOD_DISCRETE;
	}
	xppvermaj=(float)cstringmaj;
	xppvermin=(float)cstringmin;
	if(strlen(this_file)<60) {
		sprintf(pptitle,"XPP Ver %g.%g >> %s",xppvermaj,xppvermin,this_file);
	} else {
		sprintf(pptitle,"XPP Version %g.%g",xppvermaj,xppvermin);
	}
	do_meth();

	set_delay();
	rhs=my_rhs;
	init_fit_info();
	strip_saveqn();
	create_plot_list();
	auto_load_dll();

	if(XPPBatch) {
		MakeColormap();
		init_browser();
		init_all_graph();
		if_needed_load_set();
		if_needed_load_par();
		if_needed_load_ic();
		if_needed_select_sets();
		if_needed_load_ext_options();
		set_extra_graphs();
		set_colorization_stuff();

		batch_integrate();
		silent_nullclines();
		silent_dfields();
		return 0;
	}

	gtitle_text(pptitle, main_win);
	Xup=1;
	MakeColormap();
	XMapWindow(display,main_win);
	make_pops();
	make_top_buttons();
	initialize_box();
	init_browser();
	if(allwinvis==1) {
		make_new_ic_box();
		make_new_bc_box();
		make_new_delay_box();
		make_new_param_box();
		make_new_browser();
		create_eq_list();
	}

	Xup=1;
	ani_zero();
	set_extra_graphs();
	set_colorization_stuff();
	make_scrbox_lists();

	/*          MAIN LOOP             */
	test_color_info();
	if_needed_load_set();
	if_needed_load_par();
	if_needed_load_ic();
	if_needed_load_ext_options();
	if(use_ani_file) {
		new_vcr();
		get_ani_file(anifile);
	}

	if(DoTutorial==1) {
		do_tutorial();
	}
	default_window();
	do_events(min_wid,min_hgt);

	return 0;
}


void redraw_all(void) {
	redraw_dfield();
	restore(0,my_browser.maxrow);
	draw_label(draw_win);
	draw_freeze(draw_win);
	restore_on();
}


void scripty(void) {
	char scr[MAXVEC]="edf#b#b8#r";
	int bob[200];
	XKeyEvent ev;
	int i;
	int k;
	k=script_make(scr,bob);
	for(i=0;i<k;i++) {
		ev=createKeyEvent(main_win,RootWindow(display,screen),1,bob[i],0);
		XSendEvent(display,main_win,1,KeyPressMask,(XEvent *)&ev);
	}
}


/* not sure what to do with this - but it works pretty well!
   it allows you to create a KB script and send it as
   fake presses to the X11 event handler

   special keypresses are
   #t tab
   #e escape
   #b backspace
   #d delete
   #r return
  so for example to change a parameter to some numbert and then
   run the integration, you would script
   "piapp#r#b#b#b#b.12#r#rig"

   p calls parameter prompt
   iapp#r  types in iapp with a return
   #b#b#b#b deletes the current value (assuming no more than 4 numbers)
   .12# types in the number
   #r gets out of the parameter picker
   ig  runs XPP
*/
void top_button_draw(Window w) {
	if(w==TopButton[0]) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"ICs  ",5);
	}
	if(w==TopButton[1]) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"BCs  ",5);
	}
	if(w==TopButton[2]) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Delay",5);
	}
	if(w==TopButton[3]) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Param",5);
	}
	if(w==TopButton[4]) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Eqns ",5);
	}
	if(w==TopButton[5]) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Data ",5);
	}
}


/* --- Static functions -- */
static XKeyEvent createKeyEvent(Window win,
								Window winRoot, int press,
								int keycode, int modifiers) {
	XKeyEvent event;

	event.display     = display;
	event.window      = win;
	event.root        = winRoot;
	event.subwindow   = None;
	event.time        = CurrentTime;
	event.x           = 1;
	event.y           = 1;
	event.x_root      = 1;
	event.y_root      = 1;
	event.same_screen = True;
	event.keycode     = XKeysymToKeycode(display, keycode);
	event.state       = modifiers;

	if(press==1) {
		event.type = KeyPress;
	} else {
		event.type = KeyRelease;
	}
	return event;
}


static void do_events(unsigned int min_wid, unsigned int min_hgt) {
	XEvent report;

	blank_screen(main_win);
	help();
	if(RunImmediately==1) {
		run_the_commands(4);
		RunImmediately=0;
	}
	while(1) {
		XNextEvent(display,&report);
		xpp_events(report,min_wid,min_hgt);
	} /* end while */
}


static Window init_win(unsigned int bw, char *icon_name, char *win_name,
					   int x, int y, unsigned int min_wid,
					   unsigned int min_hgt, int argc, char **argv) {
	Window wine;
	int count;
	unsigned dp_h,dp_w;
	Pixmap icon_map;
	XIconSize *size_list;
	XSizeHints size_hints;
	char *display_name=NULL;

	if((display=XOpenDisplay(display_name))==NULL) {
		plintf(" Failed to open X-Display \n");
		exit(-1);
	}
	screen=DefaultScreen(display);
	if(!deleteWindowAtom) {
		deleteWindowAtom = XInternAtom(display,"WM_DELETE_WINDOW", 0);
	}
	dp_w=DisplayWidth(display,screen);
	dp_h=DisplayHeight(display,screen);
	DisplayWidth=dp_w;
	DisplayHeight=dp_h;
	if(SCALEX>dp_w) {
		SCALEX=dp_w;
	}
	if(SCALEY>dp_h) {
		SCALEY=dp_h;
	}
	wine=XCreateSimpleWindow(display,RootWindow(display,screen),
							 x,y,SCALEX,SCALEY,bw,MyForeColor,
							 MyBackColor);
	XGetIconSizes(display,RootWindow(display,screen),&size_list,&count);
	icon_map=XCreateBitmapFromData(display,wine,
								   (char*)pp_bits,pp_width,pp_height);
#ifdef X11R3
	size_hints.flags=PPosition|PSize|PMinsize;
	size_hints.x=x;
	size_hints.y=y;
	size_hints.width=width;
	size_hints.height=height;
	size_hints.min_width=min_wid;
	size_hints.min_height=min_hgt;
#else
	size_hints.flags=PPosition|PSize|PMinSize;
	size_hints.min_width=min_wid;
	size_hints.min_height=min_hgt;
#endif

#ifdef X11R3
	XSetStandardProperties(display,wine,win_name, icon_name,icon_map,argv,
						   argc,&size_hints);
#else
	{
		XWMHints wm_hints;
		XClassHint class_hints;
		XTextProperty winname,iconname;
		if(XStringListToTextProperty(&icon_name,1,&iconname)==0) {
			plintf("X error: failure for iconname\n");
			exit(-1);
		}
		if(XStringListToTextProperty(&win_name,1,&winname)==0) {
			plintf("X error: failure for winname\n");
			exit(-1);
		}

		wm_hints.initial_state=NormalState;
		wm_hints.input=True;
		wm_hints.icon_pixmap=icon_map;
		wm_hints.flags=StateHint|IconPixmapHint|InputHint;
		class_hints.res_name="base";
		class_hints.res_class=win_name;

		XSetWMProperties(display,wine,&winname,&iconname,argv,argc,&size_hints,&wm_hints,&class_hints);
		XSetWMProtocols(display, wine, &deleteWindowAtom, 1);
	}
#endif
	return(wine);
}


static void getGC(GC *gc) {
	unsigned int valuemask=0;
	XGCValues values;
	*gc=XCreateGC(display,main_win,valuemask,&values);
	XSetForeground(display,*gc,MyForeColor);
}


static int getxcolors(XWindowAttributes *win_info, XColor **colors) {
	int i, ncolors;
	*colors = (XColor *) NULL;
	TrueColorFlag=0;
	if(win_info->visual->class == TrueColor) {
		TrueColorFlag=1;
		plintf("TrueColor visual:  no colormap needed\n");
		return 0;
	} else if(!win_info->colormap) {
		plintf("no colormap associated with window\n");
		return 0;
	}

	ncolors = win_info->visual->map_entries;
	plintf("%d entries in colormap\n", ncolors);
	*colors = (XColor *) malloc (sizeof(XColor) * ncolors);

	if(win_info->visual->class == DirectColor) {
		int red, green, blue, red1, green1, blue1;

		plintf("DirectColor visual\n");

		red = green = blue = 0;
		red1   = lowbit(win_info->visual->red_mask);
		green1 = lowbit(win_info->visual->green_mask);
		blue1  = lowbit(win_info->visual->blue_mask);
		for(i=0; i<ncolors; i++) {
			(*colors)[i].pixel = red|green|blue;
			(*colors)[i].pad = 0;
			red += red1;
			if(red > win_info->visual->red_mask) {
				red = 0;
			}
			green += green1;
			if(green > win_info->visual->green_mask) {
				green = 0;
			}
			blue += blue1;
			if(blue > win_info->visual->blue_mask)   {
				blue = 0;
			}
		}
	} else {
		for(i=0; i<ncolors; i++) {
			(*colors)[i].pixel = i;
			(*colors)[i].pad = 0;
		}
	}
	XQueryColors(display, win_info->colormap, *colors, ncolors);
	return(ncolors);
}


static void init_X (void) {
	char *icon_name = "xpp";
	char *win_name = "XPPAUT";
	unsigned int x = 0, y = 0;
	unsigned int min_wid = 450, min_hgt = 360;
	char *getenv();

	char teststr[] = "The Quick Brown Fox Jumped Over The Lazy Dog?";

	if(UserMinWidth > 0) {
		min_wid = UserMinWidth;
		SCALEX = min_wid;
	}
	if(UserMinHeight> 0) {
		min_hgt = UserMinHeight;
		SCALEY = min_hgt;
	}
	if(PaperWhite == 0) {
		GrFore = White;
		GrBack = Black;
	}
	main_win = init_win (4, icon_name, win_name,
						 x, y, min_wid, min_hgt, 0, NULL);

	/*Set up foreground and background colors*/
	Black = BlackPixel (display, screen);
	White = WhitePixel (display, screen);

	if(strlen(UserBlack) != 0) {
		XColor user_col;
		XParseColor(display, DefaultColormap(display,screen), UserBlack, &user_col);
		XAllocColor(display, DefaultColormap(display,screen), &user_col);
		MyForeColor = GrFore = user_col.pixel;
		Black = MyForeColor;
	}

	if(strlen(UserWhite) != 0) {
		XColor user_col;
		XParseColor(display, DefaultColormap(display,screen), UserWhite, &user_col);
		XAllocColor(display, DefaultColormap(display,screen), &user_col);
		MyBackColor = GrBack = user_col.pixel;
		White = MyBackColor;
	}

	/*  Switch for reversed video  */
	MyForeColor = GrFore = Black;
	MyBackColor = GrBack = White;
	if(PaperWhite == 1) {
		printf("Doing swap!\n");
		char swapcol[8];
		strcpy(swapcol,UserWhite);
		strcpy(UserWhite,UserBlack);
		strcpy(UserBlack,swapcol);
		MyForeColor = GrFore = White;
		MyBackColor = GrBack = Black;
		ALREADY_SWAPPED=1;
	}

	if(strlen(UserMainWinColor) != 0) {
		XColor main_win_col;
		XParseColor(display, DefaultColormap(display,screen), UserMainWinColor, &main_win_col);
		XAllocColor(display, DefaultColormap(display,screen), &main_win_col);
		MyMainWinColor = main_win_col.pixel;
	} else {
		MyMainWinColor = MyBackColor;
	}
	XSetWindowBorder(display,main_win,MyForeColor);
	XSetWindowBackground(display,main_win,MyMainWinColor);

	if(strlen(UserDrawWinColor) != 0) {
		XColor draw_win_col;
		XParseColor(display, DefaultColormap(display,screen), UserDrawWinColor, &draw_win_col);
		XAllocColor(display, DefaultColormap(display,screen), &draw_win_col);
		MyDrawWinColor = draw_win_col.pixel;
	} else {
		MyDrawWinColor = MyBackColor;
	}
	FixWindowSize (main_win, SCALEX, SCALEY, FIX_MIN_SIZE);
	periodic = 1;

	if(DefaultDepth (display, screen) >= 8) {
		COLOR = 1;
	} else {
		COLOR = 0;
	}

	XSelectInput (display, main_win,
				  ExposureMask | KeyPressMask | ButtonPressMask |
				  StructureNotifyMask | ButtonReleaseMask | ButtonMotionMask);

	load_fonts();

	/*BETTER SUPPORT FOR VARIABLE WIDTH FONTS
  Use a statistical average to get average spacing. Some fonts don't
  or are not able to report this accurately so this is reliable way to
  get the information. If person only has variable width font on their
  system they can get by.
  The average spacing will be too small for some short strings having
  capital letters (for example "GO"). Thus, we divide by the string
  length of our test string minus 2 for a little more wiggle room.
  */

	/*   DCURXb = XTextWidth (big_font, "#", 1);
   */
	DCURXb = XTextWidth(big_font,teststr,strlen(teststr))/(strlen(teststr)-2);
	DCURYb = big_font->ascent + big_font->descent;
	CURY_OFFb = big_font->ascent - 1;

	DCURXs = XTextWidth (small_font,teststr,strlen(teststr))/(strlen(teststr)-2);
	DCURYs = small_font->ascent + small_font->descent;
	CURY_OFFs = small_font->ascent - 1;

	getGC (&gc);
	getGC (&gc_graph);
	getGC (&small_gc);
	getGC (&font_gc);

	/*User supplied */
	if(strlen(UserBGBitmap) != 0) {
		unsigned int width_return,height_return;
		int x_hot,y_hot;
		unsigned char* pixdata;
		int success = XReadBitmapFileData(UserBGBitmap,&width_return, &height_return,&pixdata,&x_hot,&y_hot);

		if(success != BitmapSuccess) {
			if(success == BitmapOpenFailed) {
				plintf("Problem reading bitmap file %s -> BitmapOpenFailed\n",UserBGBitmap);
			} else if(success == BitmapFileInvalid) {
				plintf("Problem reading bitmap file %s -> BitmapFileInvalid\n",UserBGBitmap);
			} else if(success == BitmapNoMemory) {
				plintf("Problem reading bitmap file %s -> BitmapNoMemory\n",UserBGBitmap);
			}
		} else {
			Pixmap pmap=XCreatePixmapFromBitmapData(display,main_win,(char *)pixdata,width_return,height_return,MyForeColor,MyMainWinColor,DefaultDepth(display, DefaultScreen(display)));
			XSetWindowBackgroundPixmap(display,main_win,pmap);
			XFreePixmap(display,pmap);
			XFree(pixdata);
		}
	}
	if(COLOR) {
		MakeColormap ();
	}
	set_big_font ();

	XSetFont (display, small_gc, small_font->fid);

	/*
  If the user didn't specify specifically heights and widths
  we try to set the initial size to fit everything nicely especially
  if they are using wacky fonts...
  */
	if(UserMinWidth <= 0) {
		SCALEX = 10+36*2*DCURXs + 32*DCURXs;
	}
	if(UserMinHeight <= 0) {
		SCALEY = 25*DCURYb+7*DCURYs;
	}

	XResizeWindow(display,main_win,SCALEX,SCALEY);
}


static void load_fonts(void) {
	int i;
	if((big_font=XLoadQueryFont(display,big_font_name))==NULL) {
		plintf("X Error: Failed to load big font: %s\n",big_font_name);
		exit(-1);
	}
	if((small_font=XLoadQueryFont(display,small_font_name))==NULL) {
		plintf("X Error: Failed to load small font: %s\n",small_font_name);
		exit(-1);
	}

	for(i=0;i<5;i++) {
		if((symfonts[i]=XLoadQueryFont(display,symbolfonts[i]))==NULL) {
			if(i==0||i==1) {
				symfonts[i]=small_font;
			} else {
				symfonts[i]=big_font;
			}
			avsymfonts[i]=1;
		} else {
			avsymfonts[i]=1;
			plintf(" sym %d loaded ..",i);
		}

		if((romfonts[i]=XLoadQueryFont(display,timesfonts[i]))==NULL) {
			if(i==0||i==1) {
				romfonts[i]=small_font;
			} else {
				romfonts[i]=big_font;
			}
			avromfonts[i]=1;
		} else{
			avromfonts[i]=1;
			plintf( " times %d loaded ..",i);
		}
	}
	plintf("\n");
}


static void make_pops(void) {
	int x,y;
	unsigned int h,w,bw,d;
	Window wn;
	int var_ind[3] = {
		!strlen(SLIDERVAR[0]) ? -1 : 0,
		!strlen(SLIDERVAR[1]) ? -1 : 1,
		!strlen(SLIDERVAR[2]) ? -1 : 2,
	};

	XGetGeometry(display,main_win,&wn,&x,&y,&w,&h,&bw,&d);
	create_the_menus(main_win);

	command_pop=XCreateSimpleWindow(display,main_win,0,DCURYs+4,w-2,DCURY+4,2,
									MyForeColor,
									MyBackColor);
	info_pop=XCreateSimpleWindow(display,main_win,0,h-DCURY-4,w-2,DCURY,2,
								 MyForeColor,
								 MyBackColor);
	XCreateFontCursor(display,XC_hand2);
	XSelectInput(display,command_pop,KeyPressMask|ButtonPressMask|ExposureMask);
	XSelectInput(display,info_pop,ExposureMask);
	XMapWindow(display,info_pop);
	XMapWindow(display,command_pop);
	init_grafs(16*DCURX+6,DCURYs+DCURYb+6,w-16-16*DCURX,h-6*DCURY-16);
	create_par_sliders(main_win,0,h-5*DCURY+8, var_ind);
	get_draw_area();
}


static void make_top_buttons(void) {
	int x1=2,x2=6*DCURXs+5,dx=DCURXs;
	TopButton[0]=make_fancy_window(main_win,x1,1,x2,DCURYs,1,ColorMap(20),ColorMap(TOPBUTTONCOLOR));
	x1=x1+x2+dx;
	TopButton[1]=make_fancy_window(main_win,x1,1,x2,DCURYs,1,ColorMap(20),ColorMap(TOPBUTTONCOLOR));
	x1=x1+x2+dx;

	TopButton[2]=make_fancy_window(main_win,x1,1,x2,DCURYs,1,ColorMap(20),ColorMap(TOPBUTTONCOLOR));
	x1=x1+x2+dx;

	TopButton[3]=make_fancy_window(main_win,x1,1,x2,DCURYs,1,ColorMap(20),ColorMap(TOPBUTTONCOLOR));
	x1=x1+x2+dx;

	TopButton[4]=make_fancy_window(main_win,x1,1,x2,DCURYs,1,ColorMap(20),ColorMap(TOPBUTTONCOLOR));
	x1=x1+x2+dx;

	TopButton[5]=make_fancy_window(main_win,x1,1,x2,DCURYs,1,ColorMap(20),ColorMap(TOPBUTTONCOLOR));
	x1=x1+x2+dx;
	create_user_buttons(x1,1,main_win);
}


static int script_make(char *s,int *k) {
	int l=strlen(s);
	int i=0;
	char c;
	int j=0;
	while(i<l) {
		c=s[i];
		if(c=='#') {
			i++;
			c=s[i];
			switch(c) {
			case 't':
				k[j]=65289;
				break;
			case 'r':
				k[j]=65293;
				break;
			case 'd':
				k[j]=65288;
				break;
			case 'b':
				k[j]=0xff08;
				break;
			case 'e':
				k[j]=65307;
				break;
			case 'l':
				k[j]=0xff0a;
				break;
			default:
				k[j]=32;
			}
			j++;
		} else {
			k[j]=(int)c;
			j++;
		}
		i++;
	}
	return j;
}


static void set_big_font(void) {
	DCURX=DCURXb;
	DCURY=DCURYb;
	CURY_OFF=CURY_OFFb;
	XSetFont(display,gc,big_font->fid);
}


static void test_color_info(void) {
	XColor *colors;
	XWindowAttributes xwa;
	TrueColorFlag=0;

	XGetWindowAttributes(display,main_win, &xwa);
	getxcolors(&xwa,&colors);
	if(colors) free((char *) colors);
}

static void top_button_cross(Window w,int b) {
	int i;
	for(i=0;i<6;i++) {
		if(w==TopButton[i]) {
			XSetWindowBorderWidth(display,w,b);
			return;
		}
	}
}


static void top_button_events(XEvent report) {
	switch(report.type) {
	case Expose:
	case MapNotify:
		top_button_draw(report.xany.window);
		break;
	case EnterNotify:
		top_button_cross(report.xcrossing.window,2);
		break;
	case LeaveNotify:
		top_button_cross(report.xcrossing.window,1);
		break;
	case ButtonPress:
		top_button_press(report.xbutton.window);
		break;
	}
	user_button_events(report);
}


static void top_button_press(Window w) {
	if(w==TopButton[0]) {
		make_new_ic_box();
	}
	if(w==TopButton[1]) {
		make_new_bc_box();
	}
	if(w==TopButton[2]) {
		make_new_delay_box();
	}
	if(w==TopButton[3]) {
		make_new_param_box();
	}
	if(w==TopButton[4]) {
		create_eq_list();
	}
	if(w==TopButton[5]) {
		make_new_browser();
	}
}


static void xpp_events(XEvent report,int min_wid,int min_hgt) {
	int com;
	char ch;
	int used=0;

	do_array_plot_events(report);
	txt_view_events(report);
	do_ani_events(report);
	top_button_events(report);
	switch(report.type) {
	case ConfigureNotify: /* this needs to be fixed!!! */
		/*    printf("CN %ld \n",report.xany.window); */
		resize_par_box(report.xany.window);
		resize_my_browser(report.xany.window);
		resize_eq_list(report.xany.window);
		resize_auto_window(report);
		if(report.xconfigure.window==main_win) {
			SCALEX=report.xconfigure.width;
			SCALEY=report.xconfigure.height;
			if((SCALEX<min_wid) || (SCALEY<min_hgt)) {
				/*window_size=TOO_SMALL;*/
				SCALEX=min_wid;
				SCALEY=min_hgt;
			} else {
				/*window_size=BIG_ENOUGH;*/
				XResizeWindow(display,command_pop,SCALEX-4,DCURY+1);
				XMoveResizeWindow(display,info_pop,0,SCALEY-DCURY-4,SCALEX-4,DCURY);
				resize_par_slides(SCALEY-3*DCURYs-1*DCURYb-13);
				resize_all_pops(SCALEX,SCALEY);
				redraw_all();
			}
		}
		break;
	case Expose:
	case MapNotify:
		/*  printf("E %ld \n",report.xany.window); */
		if(report.xany.window==command_pop) {
			put_command("Command:");
		}
		do_expose(report);
		break;
	case KeyPress:
		used=0;
		box_keypress(report,&used);
		if(used) {
			break;
		}
		eq_list_keypress(report,&used);
		if(used) {
			break;
		}
		my_browse_keypress(report,&used);
		if(used) {
			break;
		}
#ifdef AUTO
		auto_keypress(report,&used);
		if(used) {
			break;
		}
#endif
		ch=(char)get_key_press(&report);
		commander(ch);
		break;
	case EnterNotify:
		enter_eq_stuff(report.xcrossing.window,2);
		enter_my_browser(report,1);
		enter_slides(report.xcrossing.window,1);
		box_enter_events(report.xcrossing.window,1);
		menu_crossing(report.xcrossing.window,1);
#ifdef AUTO
		auto_enter(report.xcrossing.window,2);
#endif
		break;
	case LeaveNotify:
		enter_eq_stuff(report.xcrossing.window,1);
		enter_my_browser(report,0);
		enter_slides(report.xcrossing.window,0);
		box_enter_events(report.xcrossing.window,0);
		menu_crossing(report.xcrossing.window,0);
#ifdef AUTO
		auto_enter(report.xcrossing.window,1);
#endif
		break;
	case MotionNotify:
		do_motion_events(report);
		break;
	case ButtonRelease:
		slide_release(report.xbutton.window);
		break;
	case ButtonPress:
		if(!rotate3dcheck(report)) {
			menu_button(report.xbutton.window);
			box_buttons(report.xbutton.window);

			slide_button_press(report.xbutton.window);
			eq_list_button(report);
			my_browse_button(report);
#ifdef AUTO
			auto_button(report);
#endif

			show_position(report,&com);
		}
		break;
	} /* end switch */
}
