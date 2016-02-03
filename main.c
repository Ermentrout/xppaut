#include "main.h"

#include "aniparse.h"
#include "adj2.h"
#include "storage.h"
#include "load_eqn.h"
#include "ggets.h"
#include "many_pops.h"
#include "read_dir.h"
#include "comline.h"
#include "simplenet.h"
#include "dae_fun.h"
#include "auto_x11.h"
#include "auto_nox.h"
#include "extra.h"
#include "menudrive.h"
#include "init_conds.h"
#include "graphics.h"
#include "integrate.h"
#include "numerics.h"
#include "form_ode.h"
#include "pop_list.h"
#include "arrayplot.h"
#include "menu.h"
#include "userbut.h"
#include "color.h"
#include "eig_list.h"
#include "txtread.h"
#include "edit_rhs.h"
#include "axes2.h"
#include "do_fit.h"
#include "graf_par.h"
#include "auto_x11.h"


#include "nullcline.h"
#include "lunch-new.h"

#include "calc.h"
#include <stdlib.h> 




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



#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xproto.h>
#include <X11/Xatom.h>
/* #include <X11/bitmaps/icon> */
#include <math.h>
#include "pp.bitmap"
#include <stdio.h>
/* #include <errno.h> */
#include "help_defs.h"
#include "browse.h"
#include "struct.h"
#include <dirent.h>

#define FIX_SIZE 3
#define FIX_MIN_SIZE 2
#define FIX_MAX_SIZE 1
#define lowbit(x) ((x) & (~(x) + 1))

#define TOPBUTTONCOLOR 27

#include <X11/cursorfont.h>
#define BITMAPDEPTH 1
#define TOO_SMALL 0
#define BIG_ENOUGH 1

#include "myfonts.h"

#define cstringmaj MYSTR1
#define cstringmin MYSTR2

#ifdef NOERRNO
int errno;
#endif

int allwinvis=0;
int use_intern_sets=1;
int use_ani_file=0;
/*char anifile[256]; */
char anifile[XPP_MAX_NAME];
extern int ani_grab_flag;

float xppvermaj,xppvermin;

/*extern char this_file[100];*/
extern char this_file[XPP_MAX_NAME];
extern int METHOD,storind;
extern int (*rhs)();
extern XFontStruct *symfonts[5],*romfonts[5];
extern int avsymfonts[5],avromfonts[5];
extern int RunImmediately;
int Xup,TipsFlag=1;
Atom deleteWindowAtom=0;
int XPPBatch=0,batch_range=0;
char batchout[256];
char UserOUTFILE[256];
XKeyEvent createKeyEvent(Window w,Window wr,int p,int kc,int m);
void scripty();
int my_rhs(); 

 int DisplayHeight,DisplayWidth;
int TrueColorFlag;
char big_font_name[100],small_font_name[100];
char PlotFormat[100];

int PaperWhite=-1;
extern int DF_FLAG;
char mycommand[100];

Window TopButton[6];
Window init_win();
Window make_fancy_window();
/* Window win; */
Window draw_win;
Window make_input_strip();
Window main_win;
Window command_pop,info_pop;
GC gc, gc_graph,small_gc, font_gc;
extern int help_menu,current_pop;
unsigned int Black,White;
char UserBlack[8];
char UserWhite[8];
char UserMainWinColor[8];
char UserDrawWinColor[8];
/*char UserBGBitmap[100];*/
char UserBGBitmap[XPP_MAX_NAME];

int UserGradients=-1;
int UserMinWidth=0,UserMinHeight=0;
unsigned int MyBackColor,MyForeColor,MyMainWinColor,MyDrawWinColor;
unsigned int GrFore,GrBack;
int SCALEX,SCALEY;
extern int COLOR;
Display *display;
int screen;
extern int periodic;
int DCURYb,DCURXb,CURY_OFFb;
int DCURYs,DCURXs,CURY_OFFs;
int DCURY,DCURX,CURY_OFF;
FILE *logfile;
int XPPVERBOSE=1;
int OVERRIDE_QUIET=0;
int OVERRIDE_LOGFILE=0;
extern BROWSER my_browser;
int tfBell;

int SLIDER1=-1;
int SLIDER2=-1;
int SLIDER3=-1;
char SLIDER1VAR[20];
char SLIDER2VAR[20];
char SLIDER3VAR[20];
double SLIDER1LO=0.0;
double SLIDER2LO=0.0;
double SLIDER3LO=0.0;
double SLIDER1HI=1.0;
double SLIDER2HI=1.0;
double SLIDER3HI=1.0;
double SLIDER1INIT=0.5;
double SLIDER2INIT=0.5;
double SLIDER3INIT=0.5;

int ALREADY_SWAPPED=0;

/*Set this to 1 if you want the tutorial to come up at start-up
as default behavior
*/
int DoTutorial=0;

OptionsSet notAlreadySet;

XFontStruct *big_font,*small_font;

 int popped=0;

void do_main(argc,argv)
char **argv;
int argc;
{
/* Moved to init_X() */
/*  char *icon_name="xpp"; */
  /*char myfile[256];*/
  char myfile[XPP_MAX_NAME];
/*  char *win_name; */
  char pptitle[80];
  
  /*Track which options have not been set already*/
  notAlreadySet.BIG_FONT_NAME=1;
  notAlreadySet.SMALL_FONT_NAME=1;
  notAlreadySet.BACKGROUND=1;
  notAlreadySet.IXPLT=1;
  notAlreadySet.IYPLT=1;
  notAlreadySet.IZPLT=1;
  notAlreadySet.AXES=1;
  notAlreadySet.NMESH=1;
  notAlreadySet.METHOD=1;
  notAlreadySet.TIMEPLOT=1;
  notAlreadySet.MAXSTOR=1;
  notAlreadySet.TEND=1;
  notAlreadySet.DT=1;
  notAlreadySet.T0=1;
  notAlreadySet.TRANS=1;
  notAlreadySet.BOUND=1;
  notAlreadySet.TOLER=1;
  notAlreadySet.DELAY=1;
  notAlreadySet.XLO=1;
  notAlreadySet.XHI=1;
  notAlreadySet.YLO=1;
  notAlreadySet.YHI=1;  
  notAlreadySet.UserBlack=1;
  notAlreadySet.UserWhite=1;
  notAlreadySet.UserMainWinColor=1;
  notAlreadySet.UserDrawWinColor=1;
  notAlreadySet.UserGradients=1;
  notAlreadySet.UserBGBitmap=1;
  notAlreadySet.UserMinWidth=1;
  notAlreadySet.UserMinHeight=1;
  notAlreadySet.YNullColor=1;
  notAlreadySet.XNullColor=1;
  notAlreadySet.StableManifoldColor=1;
  notAlreadySet.UnstableManifoldColor=1;
  notAlreadySet.START_LINE_TYPE=1;
  notAlreadySet.RandSeed=1;
  notAlreadySet.PaperWhite=1;
  notAlreadySet.COLORMAP=1;
  notAlreadySet.NPLOT=1;
  notAlreadySet.DLL_LIB=1;
  notAlreadySet.DLL_FUN=1;
  notAlreadySet.XP=1;
  notAlreadySet.YP=1;
  notAlreadySet.ZP=1;
  notAlreadySet.NOUT=1;
  notAlreadySet.VMAXPTS=1;
  notAlreadySet.TOR_PER=1;
  notAlreadySet.JAC_EPS=1;
  notAlreadySet.NEWT_TOL=1;
  notAlreadySet.NEWT_ITER=1;
  notAlreadySet.FOLD=1;
  notAlreadySet.DTMIN=1;
  notAlreadySet.DTMAX=1;
  notAlreadySet.ATOL=1;
  notAlreadySet.TOL=1;
  notAlreadySet.BANDUP=1;
  notAlreadySet.BANDLO=1;
  notAlreadySet.PHI=1;
  notAlreadySet.THETA=1;
  notAlreadySet.XMIN=1;
  notAlreadySet.XMAX=1;
  notAlreadySet.YMIN=1;
  notAlreadySet.YMAX=1;
  notAlreadySet.ZMIN=1;
  notAlreadySet.ZMAX=1;
  notAlreadySet.POIVAR=1;
  notAlreadySet.OUTPUT=1;
  notAlreadySet.POISGN=1;
  notAlreadySet.POISTOP=1;
  notAlreadySet.STOCH=1;
  notAlreadySet.POIPLN=1;
  notAlreadySet.POIMAP=1;
  notAlreadySet.RANGEOVER=1;
  notAlreadySet.RANGESTEP=1;
  notAlreadySet.RANGELOW=1;
  notAlreadySet.RANGEHIGH=1;
  notAlreadySet.RANGERESET=1;
  notAlreadySet.RANGEOLDIC=1;
  notAlreadySet.RANGE=1;
  notAlreadySet.NTST=1;
  notAlreadySet.NMAX=1;
  notAlreadySet.NPR=1;
  notAlreadySet.NCOL=1;
  notAlreadySet.DSMIN=1;
  notAlreadySet.DSMAX=1;
  notAlreadySet.DS=1;
  notAlreadySet.PARMAX=1;
  notAlreadySet.NORMMIN=1;
  notAlreadySet.NORMMAX=1;
  notAlreadySet.EPSL=1;
  notAlreadySet.EPSU=1;
  notAlreadySet.EPSS=1;
  notAlreadySet.RUNNOW=1;
  notAlreadySet.SEC=1;
  notAlreadySet.UEC=1;
  notAlreadySet.SPC=1;
  notAlreadySet.UPC=1;
  notAlreadySet.AUTOEVAL=1;
  notAlreadySet.AUTOXMAX=1;
  notAlreadySet.AUTOYMAX=1;
  notAlreadySet.AUTOXMIN=1;
  notAlreadySet.AUTOYMIN=1;
  notAlreadySet.AUTOVAR=1;
  notAlreadySet.PS_FONT=1;
  notAlreadySet.PS_LW=1;   
  notAlreadySet.PS_FSIZE=1;
  notAlreadySet.PS_COLOR=1;
  notAlreadySet.FOREVER=1;
  notAlreadySet.BVP_TOL=1;
  notAlreadySet.BVP_EPS=1;
  notAlreadySet.BVP_MAXIT=1;
  notAlreadySet.BVP_FLAG=1;
  notAlreadySet.SOS=1;
  notAlreadySet.FFT=1;
  notAlreadySet.HIST=1;
  notAlreadySet.PltFmtFlag=1;
  notAlreadySet.ATOLER=1;
  notAlreadySet.MaxEulIter=1;
  notAlreadySet.EulTol=1;
  notAlreadySet.EVEC_ITER=1;
  notAlreadySet.EVEC_ERR=1;
  notAlreadySet.NULL_ERR=1;
  notAlreadySet.NEWT_ERR=1;
  notAlreadySet.NULL_HERE=1;
  notAlreadySet.TUTORIAL=1;
  notAlreadySet.SLIDER1=1;
  notAlreadySet.SLIDER2=1;
  notAlreadySet.SLIDER3=1;
  notAlreadySet.SLIDER1LO=1;
  notAlreadySet.SLIDER2LO=1;
  notAlreadySet.SLIDER3LO=1;
  notAlreadySet.SLIDER1HI=1;
  notAlreadySet.SLIDER2HI=1;
  notAlreadySet.SLIDER3HI=1;
  notAlreadySet.POSTPROCESS=1;
  notAlreadySet.HISTCOL=1;
  notAlreadySet.HISTLO=1;
  notAlreadySet.HISTHI=1;
  notAlreadySet.HISTBINS=1;
  notAlreadySet.SPECCOL=1;
  notAlreadySet.SPECCOL2=1;
  notAlreadySet.SPECWIDTH=1;
  notAlreadySet.SPECWIN=1;
  notAlreadySet.PLOTFORMAT=1;
  notAlreadySet.DFGRID=1;
  notAlreadySet.DFBATCH=1;
  notAlreadySet.NCBATCH=1;
  notAlreadySet.COLORVIA=1;
  notAlreadySet.COLORIZE=1;
    notAlreadySet.COLORHI=1;
  notAlreadySet.COLORLO=1;

  unsigned int min_wid=450,min_hgt=360;
  
/*  unsigned int x=0,y=0; */
  /*
  sprintf(myfile,"lecar.ode");
  */
  get_directory(myfile);
  
  SCALEX=640;
  SCALEY=480;
  
  Xup=0;
  sprintf(batchout,"output.dat");
  sprintf(PlotFormat,"ps");
  

  /*Read visualization environement variables here
  since some may be overridden by command line
  */
  logfile=stdout;
  check_for_quiet(argc, argv);
  /*do_vis_env();*/
  
  
   /* Old code did it this way...
   do_vis_env();
  do_comline(argc, argv);
   if(!XPPBatch)
   {
    	init_X();
   }
   
  load_eqn(); 
  */
  
  do_comline(argc, argv);

  

  /*We need to init_X here if there is no file on command line
  so that a file browser can be opened.
  */
  if(!XPPBatch)
  {    
       /*Swap out the current options for a temporary place holder
       */
       OptionsSet *tempNS = (OptionsSet*)malloc(sizeof(OptionsSet));
       *tempNS = notAlreadySet;
       /*Initialize what's needed to open a browser based on 
       the current options.
       */
       do_vis_env(); 
	 set_all_vals(); 
       init_X();
       /*       XSynchronize(display,1); */
       /*
       Now swap back the options for proper precedence ordering of options.
       */
       notAlreadySet = *tempNS;
       free(tempNS);
  }
  

  load_eqn();

  OptionsSet *tempNS = (OptionsSet*)malloc(sizeof(OptionsSet));
  *tempNS = notAlreadySet;
  set_internopts(tempNS);
  free(tempNS);
  
  init_alloc_info();
      do_vis_env();
  set_all_vals();
  
 
  /*if(!XPPBatch)
  {
       init_X();
  }*/
  
  init_alloc_info();
  set_init_guess();
  update_all_ffts(); 
 
  #ifdef AUTO 
    init_auto_win();
  #endif
  
  /* if(make_kernels()==0){
    plintf("Illegal kernel -- aborting \n");
     exit(0);
     } */
     
 
 
  if(disc(this_file))METHOD=0;
  xppvermaj=(float)cstringmaj;
  xppvermin=(float)cstringmin;
  if(strlen(this_file)<60)
    sprintf(pptitle,"XPP Ver %g.%g >> %s",xppvermaj,xppvermin,this_file);
  else sprintf(pptitle,"XPP Version %g.%g",xppvermaj,xppvermin);
/*  win_name=pptitle; */
  do_meth();
  
  
  set_delay();
  rhs=my_rhs;
  init_fit_info();
  strip_saveqn();
  create_plot_list();
  auto_load_dll();
  
if(XPPBatch){
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
    exit(0);
  }

  gtitle_text(pptitle, main_win);
  Xup=1;
  MakeColormap();


XMapWindow(display,main_win);

make_pops();

make_top_buttons();



initialize_box();

init_browser();
 if(allwinvis==1){
   make_new_ic_box();
   make_new_bc_box();
   make_new_delay_box();
   make_new_param_box();
   make_new_browser();
   create_eq_list();
 }




/* create_eq_list(); */
/* make_my_aplot("z"); */
Xup=1;
ani_zero();
 set_extra_graphs();
  set_colorization_stuff();

make_scrbox_lists();

/*   This is for testing widgets  ---    */



/*          MAIN LOOP             */
test_color_info();
  if_needed_load_set();
  if_needed_load_par();
  if_needed_load_ic();
  if_needed_load_ext_options();
if(use_ani_file)
{
	new_vcr();
	get_ani_file(anifile);
}

if(DoTutorial==1)
{
	do_tutorial();
}

default_window();

do_events(min_wid,min_hgt);


} 

void check_for_quiet(argc,argv)
char **argv;
int argc;
{
	/*First scan, check for any QUIET option set...*/
	int i = 0;
	/*Allow for multiple calls to the QUIET and LOGFILE options 
	on the command line. The last setting is the one that will stick.
	Settings of logfile and quiet in the xpprc file will be ignored
	if they are set on the command line.
	*/ 
	int quiet_specified_once=0;
	int logfile_specified_once=0;
	
	for(i=1;i<argc;i++)
	{
 	       if (strcmp(argv[i],"-quiet")==0)
	       {
	       	       set_option("QUIET",argv[i+1],1,NULL);
		       quiet_specified_once=1;
     		       i++;
	       }
	       else if (strcmp(argv[i],"-logfile")==0)
	       {
		       set_option("LOGFILE",argv[i+1],1,NULL);
		       logfile_specified_once = 1;
     		       i++;
	       }
	}
	/*If -quiet or -logfile were specified at least once on the command line
	we lock those in now...
	*/
	if (quiet_specified_once == 1)
	{
		OVERRIDE_QUIET=1;
	}
	if (logfile_specified_once == 1)
	{
		OVERRIDE_LOGFILE=1;
	}
}


void do_vis_env()
{
  set_X_vals();
  check_for_xpprc();
  set_internopts_xpprc_and_comline();
  
}


void init_X ()
{
  char *icon_name = "xpp";
  char *win_name = "XPPAUT";
  unsigned int x = 0, y = 0;
  unsigned int min_wid = 450, min_hgt = 360; 
  char *getenv();
  
  char teststr[] = "The Quick Brown Fox Jumped Over The Lazy Dog?";
  
   
  if (UserMinWidth > 0)
  {
  	min_wid = UserMinWidth;
	SCALEX = min_wid;
  }
  
  if (UserMinHeight> 0)
  {
  	min_hgt = UserMinHeight;
	SCALEY = min_hgt;
  }

  if (PaperWhite == 0)
    {
      GrFore = White;
      GrBack = Black;
    }


  main_win = init_win (4, icon_name, win_name,
                       x, y, min_wid, min_hgt, 0, NULL);
		       
   /*Set up foreground and background colors*/
  
 
  Black = BlackPixel (display, screen);
  White = WhitePixel (display, screen);
  
  



  if (strlen(UserBlack) != 0)
  {       
 	XColor user_col;

  	XParseColor(display, DefaultColormap(display,screen), UserBlack, &user_col);
  	XAllocColor(display, DefaultColormap(display,screen), &user_col);
  
  	MyForeColor = GrFore = user_col.pixel;
	Black = MyForeColor;
  }
  
  if (strlen(UserWhite) != 0)
  {
  	XColor user_col;

  	XParseColor(display, DefaultColormap(display,screen), UserWhite, &user_col);
  	XAllocColor(display, DefaultColormap(display,screen), &user_col);
  
  	MyBackColor = GrBack = user_col.pixel;
	White = MyBackColor;
	     
  }
  
  
  
  /*  Switch for reversed video  */
  MyForeColor = GrFore = Black;
  MyBackColor = GrBack = White;
 /*  This is reversed  
  MyForeColor=White;
  MyBackColor=Black; */
  
  if (PaperWhite == 1) /*Respect the swapping implied by the -white option.*/
  {
  	  printf("Doing swap!\n");
	  char swapcol[8];
	  strcpy(swapcol,UserWhite);
	  strcpy(UserWhite,UserBlack);
	  strcpy(UserBlack,swapcol);


	  MyForeColor = GrFore = White;
	  MyBackColor = GrBack = Black;
	  
	  ALREADY_SWAPPED=1;
  }	
	
	
  if (strlen(UserMainWinColor) != 0)
  {
  	XColor main_win_col;

  	XParseColor(display, DefaultColormap(display,screen), UserMainWinColor, &main_win_col);
  	XAllocColor(display, DefaultColormap(display,screen), &main_win_col);
  
    	MyMainWinColor = main_win_col.pixel;
  }
  else
  {
  	MyMainWinColor = MyBackColor;
  }
  
  XSetWindowBorder(display,main_win,MyForeColor);
  XSetWindowBackground(display,main_win,MyMainWinColor); 
  
  if (strlen(UserDrawWinColor) != 0)
  {
  	XColor draw_win_col;
  	XParseColor(display, DefaultColormap(display,screen), UserDrawWinColor, &draw_win_col);
  	XAllocColor(display, DefaultColormap(display,screen), &draw_win_col);
  
    	MyDrawWinColor = draw_win_col.pixel;
  }
  else
  { 
  	MyDrawWinColor = MyBackColor;
  }
  
  /* win = main_win; */
  FixWindowSize (main_win, SCALEX, SCALEY, FIX_MIN_SIZE);
  periodic = 1;
  if (DefaultDepth (display, screen) >= 8)
  {
    COLOR = 1;
  }
  else
  {
    COLOR = 0;
  }
  /* if(DefaultDepth(display,screen)==4){
     map16();
     COLOR=1;
     } */


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
  
  /*  DCURXs = XTextWidth (small_font, "#", 1);
   */
 
  DCURXs = XTextWidth (small_font,teststr,strlen(teststr))/(strlen(teststr)-2);
  

  DCURYs = small_font->ascent + small_font->descent;
  CURY_OFFs = small_font->ascent - 1;

  getGC (&gc);
  getGC (&gc_graph);
  getGC (&small_gc);
  getGC (&font_gc);
  

  if (strlen(UserBGBitmap) != 0) /*User supplied */
  {
  	  /*Pixmap pmap=XCreatePixmapFromBitmapData(display,main_win,lines_bits,lines_width,lines_height,MyForeColor,MyBackColor,DefaultDepth(display, DefaultScreen(display)));
	 */
	 unsigned int width_return,height_return;
	 int x_hot,y_hot;
	 unsigned char* pixdata;
	 
	 int success = XReadBitmapFileData(UserBGBitmap,&width_return, &height_return,&pixdata,&x_hot,&y_hot);
	
	 if (success != BitmapSuccess)
	 {
	 	if (success == BitmapOpenFailed)
		{
		plintf("Problem reading bitmap file %s -> BitmapOpenFailed\n",UserBGBitmap);
		}
		else if (success == BitmapFileInvalid)
		{
		plintf("Problem reading bitmap file %s -> BitmapFileInvalid\n",UserBGBitmap);
	 	}
		else if (success == BitmapNoMemory)
		{
		plintf("Problem reading bitmap file %s -> BitmapNoMemory\n",UserBGBitmap);
	 	}
	 } 
	 else
	 { 
		Pixmap pmap=XCreatePixmapFromBitmapData(display,main_win,(char *)pixdata,width_return,height_return,MyForeColor,MyMainWinColor,DefaultDepth(display, DefaultScreen(display)));
		XSetWindowBackgroundPixmap(display,main_win,pmap); 
		XFreePixmap(display,pmap);
		XFree(pixdata);
	 }    
  }
         
	 
  if (COLOR)
    MakeColormap ();

  set_big_font ();
/* set_small_font(); */

  XSetFont (display, small_gc, small_font->fid);
  /*make_pops();*/
  
  /*
  If the user didn't specify specifically heights and widths
  we try to set the initial size to fit everything nicely especially
  if they are using wacky fonts...
  */
  if (UserMinWidth <= 0)
  {
	SCALEX = 10+36*2*DCURXs + 32*DCURXs;
  }
  
  if (UserMinHeight <= 0)
  {
	SCALEY = 25*DCURYb+7*DCURYs;
  }
  
  XResizeWindow(display,main_win,SCALEX,SCALEY);
  /*FixWindowSize (main_win, SCALEX, SCALEY, FIX_MIN_SIZE);
  */
}


void set_big_font()
{
 DCURX=DCURXb;
 DCURY=DCURYb;
 CURY_OFF=CURY_OFFb;
 XSetFont(display,gc,big_font->fid);
}

void set_small_font()
{
  DCURX=DCURXs;
 DCURY=DCURYs;
 CURY_OFF=CURY_OFFs;
 XSetFont(display,gc,small_font->fid);
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
int script_make(char *s,int *k)
{
  int l=strlen(s);
  int i=0;
  char c;
  int j=0;
  while(i<l){
    c=s[i];
    if(c=='#'){
      i++;
      c=s[i];
      switch(c){
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
    }
    else
      {
	k[j]=(int)c;
	j++;
      }
    i++;
  }
  return j;
}



void scripty()
{
  /*  char scr[100]="piapp#r#b#b#b#b.12#r#rig"; */
  /* char scr[100]="eir#t"; */
  char scr[100]="edf#b#b8#r";
  int bob[200];
  XKeyEvent ev;
  int i;
  int k;
  k=script_make(scr,bob);
  for(i=0;i<k;i++){
    ev=createKeyEvent(main_win,RootWindow(display,screen),1,bob[i],0);
    XSendEvent(display,main_win,1,KeyPressMask,(XEvent *)&ev);
  }
}
    

void xpp_events(XEvent report,int min_wid,int min_hgt)
{
/*int window_size,com;
 */
 
 int com;
 char ch;
 
 int used=0;

/*  put_command("Command:");  */
 
   do_array_plot_events(report);
   txt_view_events(report);
   do_ani_events(report);
   top_button_events(report);
 switch(report.type){
 /* case ClientMessage:
	  if(report.xclient.data.l[0]==deleteWindowAtom){
        
	  break;
	  }
	  break; */
 
 case ConfigureNotify: /* this needs to be fixed!!! */
   /*    printf("CN %ld \n",report.xany.window); */
  resize_par_box(report.xany.window);  
  resize_my_browser(report.xany.window);
   resize_eq_list(report.xany.window);
   resize_auto_window(report); 
  if(report.xconfigure.window==main_win){

  SCALEX=report.xconfigure.width;
  SCALEY=report.xconfigure.height;
  if((SCALEX<min_wid) || (SCALEY<min_hgt)){
   /*window_size=TOO_SMALL;*/
  SCALEX=min_wid;
  SCALEY=min_hgt;}
  else {
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
	if(report.xany.window==command_pop)put_command("Command:");
     do_expose(report);


   break;
 case KeyPress:
   used=0;
                box_keypress(report,&used);
		if(used)break;
                eq_list_keypress(report,&used);
                if(used)break;
                my_browse_keypress(report,&used);
		 if(used)break;
#ifdef AUTO		
	        auto_keypress(report,&used);
	         if(used)break;
#endif
	             ch=(char)get_key_press(&report);
		commander(ch);
	    
	    /* do_key_stuff(report); */


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
   /* check_box_cursor(); */
   if(!rotate3dcheck(report)){
	        menu_button(report.xbutton.window);
	        /* box_select_events(report.xbutton.window,&i1); */
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


void do_events(min_wid,min_hgt)
unsigned int min_wid,min_hgt;  
{
 XEvent report;
 
  blank_screen(main_win);
  help();
  if(RunImmediately==1){
    run_the_commands(4);
    RunImmediately=0;
  }
while(1)
{
/*  put_command("Command:");  */
 XNextEvent(display,&report);
 xpp_events(report,min_wid,min_hgt);
 } /* end while */
}


void bye_bye()
{
   int i;
   yes_reset_auto();
   XUnloadFont(display,big_font->fid);
   XUnloadFont(display,small_font->fid);
   for(i=0;i<5;i++){
     if(avsymfonts[i])XUnloadFont(display,symfonts[i]->fid);
     if(avromfonts[i])XUnloadFont(display,romfonts[i]->fid);
   }
   XFreeGC(display,gc);
   XCloseDisplay(display);
   exit(1);
}


void clr_scrn()
{
	  blank_screen(draw_win);
			 restore_off();
			 do_axes();
}

void redraw_all()
{
    redraw_dfield();
  restore(0,my_browser.maxrow);
  draw_label(draw_win);
  draw_freeze(draw_win);
  restore_on();
}


void commander(ch)
 char ch;
 {
   switch(help_menu)
  {
      case MAIN_HELP:
  	{ 
          switch(ch)
		{
		case 'i':flash(0);
			/*  initial data  */
                        ini_data_menu();

				 flash(0);
			 break;
		case 	'c':flash(1);
			 /* continue */
		  cont_integ(); 

			 
		  flash(1);
			 break;
		case 'n':flash(2);
			 /*nullclines */
		  /* test_color_info(); */
                          new_clines(); 
			 flash(2);
			 break;
		case 'd':flash(3);
			 /*dir fields */
		  direct_field(); 

			 flash(3);
			 break;
		case 'w':flash(4);
			 /* window */
		         window_zoom();
			 flash(4);
			 break;
		case 'a':flash(5);
			 /*phase-space */
		         do_torus();
			 flash(5);
			 break;
		case 'k':flash(6);
			 /*kinescope */
			 do_movie();
			 flash(6);
			 break;
		case 'g' :flash(7);
			  flash(7);
			  add_a_curve();
			  break;
		case 'u':flash(8);
			 flash(8);
			 help_num();
			 break;
		case 'f':flash(9);
			 /* files */
			 flash(9);
			help_file();
			break;
		case 'p':flash(10);
			 /*parameters */
		  /* change_par(-1); */
		         new_param();
			 flash(10);
			 break;
		case 'e':flash(11);
			 /*erase */
		         clear_screens();
			 flash(11);
			 break;
                case 'h':
		case 'm':
                         do_windows();
			 /*half windows */
			 flash(12);
			 break;
		case 't':flash(13);
			 /*text */
			 do_gr_objs();
			 flash(13);
			 break;
		case 's':flash(14);
			 /*sing pts */
		         find_equilibrium();
			 flash(14);
			 break;
		case 'v':flash(15);
			 /*view_axes */
		        change_view();
			 flash(15);
			 break;
		case 'b':
		        find_bvp();
		        break;
                        
		case 'x':flash(16);
			 /*x vs t */
		  x_vs_t();
			 flash(16);
			 break;
		case 'r':flash(17);
			 /*restore*/
			redraw_them_all();
			 flash(17);
			 break;
		case '3': get_3d_par();
			 break;
		case 'y':
		  /* scripty();  */
		  break;


                 } /* End main switch  */
		
   
            } /* MAIN HELP ENDS  */ break;

  case NUM_HELP:{
    get_num_par(ch);
  } /* end num case   */ break;
  
  
  
  case FILE_HELP: {
          switch(ch){
	  case 't': do_transpose();
	            break;
          case 'g': get_intern_set();
	            break;
	  case 'i':  TipsFlag=1-TipsFlag;
	            break;
		case 'p': flash(0);
			/* file stuff */
		  /* do_info(stdout); */
		  make_txtview();
			flash(0);
			break;
		case 'w': flash(1);
			/* write set */
		       do_lunch(0);
			flash(1);
			break;
		case 's': flash(2);
			/* make eqn */
		       file_inf();
			flash(2);
			break;
		case 'a': flash(3);
			/* AUTO !! */
#ifdef AUTO
		        do_auto_win();
#endif
			flash(3);
			break;
		case 'c': flash(4);
			/* calculator */
		       q_calc();  
		  
			flash(4);
			break;
		case 'r': flash(5);
			/* read set */
		        do_lunch(1);
			flash(5);
			break;
		case 'e': flash(6);
			/* script */
                        edit_menu();
	               /*  Insert generic code here ...  */
			flash(6);
			break;
                case 'b': 
		        tfBell=1-tfBell;
		        break;
                 case 'h':
		   /* c_hints(); */ 
		   xpp_hlp();   
		    /*	   make_key_stroke(); */
 

		        break;
		case 'q': flash(7);
			if(yes_no_box())bye_bye();
			 
			/*  if(are_you_sure())bye_bye();  */
			
			/* quit */
			 flash(7);
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

/*
	
Window make_unmapped_window(root,x,y,width,height,bw)
	Window root;
	int x,y,width,height,bw;
	{
	 Window win;
	 win=XCreateSimpleWindow(display,root,x,y,width,height,
		bw,MyForeColor,MyBackColor);
	 XSelectInput(display,win,ExposureMask|KeyPressMask|ButtonPressMask|
              StructureNotifyMask|ButtonReleaseMask|ButtonMotionMask);
         return(win);
         }

	 */

XKeyEvent createKeyEvent(Window win,
                           Window winRoot, int press,
                           int keycode, int modifiers)
{
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

   if(press==1)
      event.type = KeyPress;
   else
      event.type = KeyRelease;

   return event;
}

 Window init_win(bw,icon_name,win_name,
                 x,y,min_wid,min_hgt,argc,argv)
 int argc;
 int x,y;
 char **argv;
 unsigned int min_wid,min_hgt,bw;
 char *icon_name,*win_name;
 {
/*  XSetWindowAttributes xswa;
 XWindowAttributes xwa;
  */
  Window wine;
  int count;
  unsigned dp_h,dp_w;
  Pixmap icon_map;
 XIconSize *size_list;
  XSizeHints size_hints;
 char *display_name=NULL;

 if((display=XOpenDisplay(display_name))==NULL){
  plintf(" Failed to open X-Display \n");
   exit(-1);
 }
 /*  printf("Display=%s\n",display_name); */
/*   Remove after debugging is done */   
 /* XSynchronize(display,1); */
 screen=DefaultScreen(display);
 if (!deleteWindowAtom) {
	deleteWindowAtom = XInternAtom(display,"WM_DELETE_WINDOW", 0);
    }
 dp_w=DisplayWidth(display,screen);
 dp_h=DisplayHeight(display,screen);
 DisplayWidth=dp_w;
 DisplayHeight=dp_h;
 if(SCALEX>dp_w)SCALEX=dp_w;
 if(SCALEY>dp_h)SCALEY=dp_h;
 wine=XCreateSimpleWindow(display,RootWindow(display,screen),
     x,y,SCALEX,SCALEY,bw,MyForeColor,
     MyBackColor);
/*  xswa.override_redirect=1;
 XChangeWindowAttributes(display,wine,CWOverrideRedirect,&xswa); */
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
 if(XStringListToTextProperty(&icon_name,1,&iconname)==0)
 {
plintf("X error: failure for iconname\n");
 exit(-1);
 }
if(XStringListToTextProperty(&win_name,1,&winname)==0)
 {
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

void top_button_draw(Window w)
{
  if(w==TopButton[0])
    XDrawString(display,w,small_gc,5,CURY_OFFs,"ICs  ",5);
if(w==TopButton[1])
    XDrawString(display,w,small_gc,5,CURY_OFFs,"BCs  ",5);
if(w==TopButton[2])
    XDrawString(display,w,small_gc,5,CURY_OFFs,"Delay",5);
if(w==TopButton[3])
    XDrawString(display,w,small_gc,5,CURY_OFFs,"Param",5);
if(w==TopButton[4])
    XDrawString(display,w,small_gc,5,CURY_OFFs,"Eqns ",5);
if(w==TopButton[5])
    XDrawString(display,w,small_gc,5,CURY_OFFs,"Data ",5);
}
 
void top_button_cross(Window w,int b)
{
  int i;
  for(i=0;i<6;i++)
    if(w==TopButton[i]){
      XSetWindowBorderWidth(display,w,b);
      return;
    }
}
void top_button_press(Window w)
{
  if(w==TopButton[0]){
    make_new_ic_box();

  }
  if(w==TopButton[1]){
    make_new_bc_box();

  }
  if(w==TopButton[2]){
    make_new_delay_box();

  }
  if(w==TopButton[3]){
    make_new_param_box();

  }
  if(w==TopButton[4]){
    create_eq_list(); 

  }
  if(w==TopButton[5]){
    make_new_browser();

  }
}
    
void top_button_events(XEvent report)
{
  switch(report.type){
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
void make_top_buttons()
{
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
void getGC(gc)
GC *gc;
{
 unsigned int valuemask=0;
 XGCValues values;
 /* unsigned int lw=6;
  int ls=LineOnOffDash;
 int cs=CapRound;
 int js=JoinRound;
 int dash_off=0;
 static char dash[]={12,24};
 int ll=2; */
 *gc=XCreateGC(display,main_win,valuemask,&values);
 XSetForeground(display,*gc,MyForeColor);
/* XSetLineAttributes(display,*gc,lw,ls,cs,js);
 XSetDashes(display,*gc,dash_off,dash,ll); */
}

void load_fonts()
 {
 
  int i;
  /*printf("\n\nFONTS %s %s \n",big_font_name,small_font_name); 
   */
   if((big_font=XLoadQueryFont(display,big_font_name))==NULL)
 {
 plintf("X Error: Failed to load big font: %s\n",big_font_name);
  exit(-1);
 }
 
  if((small_font=XLoadQueryFont(display,small_font_name))==NULL)
 {
 plintf("X Error: Failed to load small font: %s\n",small_font_name);
  exit(-1);
 }
  

  for(i=0;i<5;i++){
    if((symfonts[i]=XLoadQueryFont(display,symbolfonts[i]))==NULL){
      if(i==0||i==1)
	symfonts[i]=small_font;
      else
	symfonts[i]=big_font;
      avsymfonts[i]=1;
    }
    else{ 
      avsymfonts[i]=1;
     plintf(" sym %d loaded ..",i);
    }
    
    if((romfonts[i]=XLoadQueryFont(display,timesfonts[i]))==NULL){
       if(i==0||i==1)
	romfonts[i]=small_font;
      else
	romfonts[i]=big_font;
      avromfonts[i]=1;
    }
    else{
      avromfonts[i]=1;
     plintf( " times %d loaded ..",i);
    }
  }
 plintf("\n");

 
}


void make_pops()

{  int x,y;
  unsigned int h,w,bw,d;
   Window wn;
   /*Cursor cursor;*/
   XGetGeometry(display,main_win,&wn,&x,&y,&w,&h,&bw,&d);
   /* menu_pop=XCreateSimpleWindow(display,main_win,
     0,DCURY+6,16*DCURX,22*DCURY,2,MyForeColor,
     MyBackColor);
     */
  create_the_menus(main_win);
 command_pop=XCreateSimpleWindow(display,main_win,0,DCURYs+4,w-2,DCURY+4,2,
MyForeColor,
     MyBackColor);
 info_pop=XCreateSimpleWindow(display,main_win,0,h-DCURY-4,w-2,DCURY,2,
MyForeColor,
     MyBackColor);
 /*cursor=XCreateFontCursor(display,XC_hand2);*/
 XCreateFontCursor(display,XC_hand2);
 /* XDefineCursor(display,menu_pop,cursor); */
 /* XSelectInput(display,menu_pop,KeyPressMask|ButtonPressMask|ExposureMask); */
XSelectInput(display,command_pop,KeyPressMask|ButtonPressMask|ExposureMask);
XSelectInput(display,info_pop,ExposureMask);
 XMapWindow(display,info_pop);
 XMapWindow(display,command_pop);
 /* XMapWindow(display,menu_pop); */
 init_grafs(16*DCURX+6,DCURYs+DCURYb+6,w-16-16*DCURX,h-6*DCURY-16);
 create_par_sliders(main_win,0,h-5*DCURY+8);
  get_draw_area(); 
 
}


void FixWindowSize(w,width,height,flag)
     Window w;
     int width,height,flag;
{
 XSizeHints size_hints;
 switch(flag){
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




int getxcolors(win_info, colors)
     XWindowAttributes *win_info;
     XColor **colors;
{

  int i, ncolors;
 

  *colors = (XColor *) NULL;
  TrueColorFlag=0;
  if (win_info->visual->class == TrueColor) {
    TrueColorFlag=1;
   plintf("TrueColor visual:  no colormap needed\n"); 
    return 0;
  }

  else if (!win_info->colormap) {
   plintf("no colormap associated with window\n");
    return 0;
  }

  ncolors = win_info->visual->map_entries;
 plintf("%d entries in colormap\n", ncolors);

  *colors = (XColor *) malloc (sizeof(XColor) * ncolors);


  if (win_info->visual->class == DirectColor) {
    int red, green, blue, red1, green1, blue1;

   plintf("DirectColor visual\n");

    red = green = blue = 0;
    red1   = lowbit(win_info->visual->red_mask);
    green1 = lowbit(win_info->visual->green_mask);
    blue1  = lowbit(win_info->visual->blue_mask);
    for (i=0; i<ncolors; i++) {
      (*colors)[i].pixel = red|green|blue;
      (*colors)[i].pad = 0;
      red += red1;
      if (red > win_info->visual->red_mask)     red = 0;
      green += green1;
      if (green > win_info->visual->green_mask) green = 0;
      blue += blue1;
      if (blue > win_info->visual->blue_mask)   blue = 0;
    }
  }
  else {
    for (i=0; i<ncolors; i++) {
      (*colors)[i].pixel = i;
      (*colors)[i].pad = 0;
    }
  }

  XQueryColors(display, win_info->colormap, *colors, ncolors);

  return(ncolors);
}
    






void test_color_info()
{
 XColor *colors;
 XWindowAttributes xwa;
 /*int n;
 */
 TrueColorFlag=0;
 
 XGetWindowAttributes(display,main_win, &xwa);
 /*n=getxcolors(&xwa,&colors);
 */
 getxcolors(&xwa,&colors);
 
if (colors) free((char *) colors);
}


/*   this is stuff for ndrand48(), nsrand48()  
     uncomment if you dont have it 
     RAND_MAX is either 2^16-1 or 2^32-1
     it is usually in stdlib.h

#ifdef NORAND48
double ndrand48()
{
  return (double)rand()/(double)RAND_MAX;
}

nsrand48(seed)
int seed;
{
  srand(seed);
}

#endif

*/





