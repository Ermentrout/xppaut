#include "arrayplot.h"
#include "array_print.h"

#include <stdlib.h> 
#include <string.h>
/*   routines for plotting arrays as functions of time  

     makes a window 
     of  N X M pixels 
     user specifies   starting variable  x0 and ending variable xn
                      starting time  ending time 
                      max var  min var

                                TITLE

                   [Kill]  [Edit]  [Print]  [Style] [Fit] [Range]
   ________________________________________________________
           1 |  |      tic marks              |  | N
            ---------------------------------------
     T0
          -                                                   MAX
          -                                                   ---
          -                                                   | |
                                                              | |
                                                              | |
                                                              | | 
                                                              | | 
                                                              | |  
          - 
          -                                                   MIN
     TN     ---------------------------------------                   
 

    and it creates a color plot 

*/
#include "ggets.h"
#include "main.h"
#include "many_pops.h"
#include "pop_list.h"
#include "kinescope.h"
#include "scrngif.h"
#include "lunch-new.h"
#include "color.h"
#include "init_conds.h"
#include "load_eqn.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include <math.h>
#include <math.h>
#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif
#include "array.bitmap"
#include "xpplim.h"
#define READEM 1
#include "browse.h"
#define MAX_LEN_SBOX 25
#define FIRSTCOLOR 30
#define FIX_MIN_SIZE 2
extern int COLOR;
extern unsigned int GrFore,GrBack;
/*extern char this_file[100];*/
extern char this_file[XPP_MAX_NAME];
extern Display *display;
extern int DCURX,DCURXs,DCURY,DCURYs,CURY_OFFs,CURY_OFF,color_total,screen;
extern GC gc, small_gc,gc_graph;
double atof();
extern char uvar_names[MAXODE][12];
extern BROWSER my_browser;
int aplot_range;
int aplot_range_count=0;
char aplot_range_stem[256]="rangearray";
int aplot_still=1,aplot_tag=0;
APLOT aplot;
extern Window draw_win;
int plot3d_auto_redraw=0;
FILE *ap_fp;
GC aplot_gc;
int first_aplot_press;
int do_range(double *, int);
extern double MyData[MAXODE];

#define MYMASK  (ButtonPressMask 	|\
                ButtonReleaseMask |\
		KeyPressMask		|\
		ExposureMask		|\
		StructureNotifyMask	|\
		LeaveWindowMask		|\
		EnterWindowMask)



void draw_one_array_plot(char *bob)
{
  char filename[300];
 
  redraw_aplot(aplot);
  if(aplot_tag)tag_aplot(bob);
  XFlush(display);
  sprintf(filename,"%s.%d.gif",aplot_range_stem,aplot_range_count);
  gif_aplot_all(filename,aplot_still);
  aplot_range_count++;

}

void set_up_aplot_range()
{ 
  static char *n[]={"Basename","Still(1/0)","Tag(0/1)"};
  char values[3][MAX_LEN_SBOX];
  int status;
  double *x;
 sprintf(values[0],"%s",aplot_range_stem);
 sprintf(values[1],"%d",aplot_still);
 sprintf(values[2],"%d",aplot_tag);
 status=do_string_box(3,3,1,"Array range saving",n,values,28); 
 if(status!=0){
   sprintf(aplot_range_stem,"%s",values[0]);
   aplot_still=atoi(values[1]);
   aplot_tag=atoi(values[2]);
 aplot_range=1;
 aplot_range_count=0;
 x=&MyData[0];
 do_range(x,0);
 }
}
void fit_aplot()
{
double zmax,zmin;
 scale_aplot(&aplot,&zmax,&zmin);
  aplot.zmin=zmin;
  aplot.zmax=zmax;
  redraw_aplot(aplot);

}
void optimize_aplot(int *plist)
{
  int i0=plist[0]-1;
  int i1=plist[1]-1;
  int nr,ns;
  double zmax,zmin;
  int nrows=my_browser.maxrow;
  int ncol=i1+1-i0;
  if(ncol<2||nrows<2)return;
  make_my_aplot("Array!");

  aplot.index0=i0+1;
  strcpy(aplot.name,uvar_names[i0]);
  aplot.nacross=ncol;
  nr=201;
  if(nrows<nr)
    nr=nrows;
  aplot.ndown=nr;
  ns=nrows/nr;
  aplot.nskip=ns;
  aplot.ncskip=1;
  scale_aplot(&aplot,&zmax,&zmin);
  aplot.zmin=zmin;
  aplot.zmax=zmax;
  aplot.plotdef=1;
  reset_aplot_axes(aplot);
  redraw_aplot(aplot);
}
  
  
  
void make_my_aplot(name)
     char *name;
{
  if(aplot.alive==1)return;
  create_arrayplot(&aplot,name,name);
}

void scale_aplot(ap,zmax,zmin)
APLOT *ap;
double *zmax,*zmin;
{
  int i,j,ib,jb,row0=ap->nstart,col0=ap->index0;
  int nrows=my_browser.maxrow;
  double z;
  ib=col0;
  jb=row0;
  *zmax=my_browser.data[ib][jb];
  *zmin=*zmax;
  for(i=0;i<ap->nacross/ap->ncskip;i++){
      ib=col0+i*ap->ncskip;
      if(ib<=my_browser.maxcol);{
	for(j=0;j<ap->ndown;j++){
	  jb=row0+ap->nskip*j;
	  if(jb<nrows&&jb>=0){
	    z=my_browser.data[ib][jb];
	    if(z<*zmin)*zmin=z;
	    if(z>*zmax)*zmax=z;
	  }
	}
      }
  }
  if(*zmin>=*zmax)
    *zmax=fabs(*zmin)+1+*zmin;
 
}

void init_arrayplot(ap)
APLOT *ap;
{
 ap->height=400;
 ap->width=400;
 ap->zmin=0.0;
 ap->zmax=1.0;
 ap->alive=0;
 ap->plotdef=0;
 ap->index0=1;
 ap->indexn=0;
 ap->nacross=1;
 ap->ndown=50;
 ap->nstart=0;
 ap->nskip=8;
 ap->ncskip=1;
 ap->tstart=0.0;
 ap->tend=20.0;
 strcpy(ap->filename,"output.ps");
 strcpy(ap->xtitle,"index");
 strcpy(ap->ytitle,"time");
 strcpy(ap->bottom,"");
 ap->type=-1;
}

void expose_aplot(w)
     Window w;
{
  if(aplot.alive)
    display_aplot(w,aplot);
}


void do_array_plot_events(ev)
     XEvent ev;
{
  int x,y;
  if(aplot.alive==0)return;
  switch(ev.type){
 /* case Expose:
  case MapNotify:
    display_aplot(ev.xany.window,aplot);
    break;
  */
  case MotionNotify:
    if(ev.xany.window==aplot.wplot){
    /*printf("%d\n",ev.xmotion.y-first_aplot_press); */
    aplot.nstart=aplot.nstart-ev.xmotion.y+first_aplot_press;
    if(aplot.nstart<0)aplot.nstart=0;
    redraw_aplot(aplot);}
    break;
  case ConfigureNotify:
    if(ev.xconfigure.window!=aplot.base)return;
    x=ev.xconfigure.width;
    y=ev.xconfigure.height;
    aplot.width=x;
    aplot.height=y;
    aplot.ploth=y-55;
    aplot.plotw=x-30-10*DCURXs;
    XMoveResizeWindow(display,aplot.wplot,20+10*DCURXs,45,aplot.plotw,
		      aplot.ploth);
    break;
  case EnterNotify:
    wborder(ev.xexpose.window,2,aplot);
    break;
  case LeaveNotify:
    wborder(ev.xexpose.window,1,aplot);
    break;
  case ButtonPress:
    if(ev.xany.window==aplot.wplot)
      first_aplot_press=ev.xbutton.y;
    /*apbutton(ev.xbutton.window,aplot);*/
      apbutton(ev.xbutton.window);
    break;
  }
}

void wborder(w,i,ap)
     Window w;
     int i;
     APLOT ap;
{
 /* if(w==ap.wedit||w==ap.wprint||w==ap.wkill||w==ap.wstyle||w==ap.wredraw) */
  if(w==ap.wedit||w==ap.wprint||w==ap.wclose||w==ap.wredraw||w==ap.wgif||w==ap.wrange||w==ap.wfit)
    XSetWindowBorderWidth(display,w,i);
}

void destroy_aplot()
{
  aplot.alive=0;
  waitasec(ClickTime);
  XDestroySubwindows(display,aplot.base);
    XDestroyWindow(display,aplot.base);
}

void init_my_aplot()
{
 init_arrayplot(&aplot);
}


void create_arrayplot(ap,wname,iname)
     APLOT *ap;
     char *wname,*iname;

{
  Window base;
  int width,height;
  unsigned int valuemask=0;
  XGCValues values;
  XTextProperty winname,iconname;
  XSizeHints size_hints;
  /* init_arrayplot(ap); */ 
  width=ap->width;
  height=ap->height;
  base=make_plain_window(RootWindow(display,screen),0,0,ap->width,ap->height,1);
  ap->base=base;
  XSelectInput(display,base,ExposureMask|KeyPressMask|ButtonPressMask|
	       StructureNotifyMask);
  XStringListToTextProperty(&wname,1,&winname);
  XStringListToTextProperty(&iname,1,&iconname);
  
  size_hints.flags=PPosition|PSize|PMinSize;
  size_hints.x=0;
  size_hints.y=0;
  size_hints.min_width=width;
  size_hints.min_height=height;
   XSetWMProperties(display,base,&winname,&iconname,NULL,0,&size_hints,NULL,NULL);
  FixWindowSize(base,width,height,FIX_MIN_SIZE);
  make_icon((char *)array_bits,array_width,array_height,base);

  ap->wredraw=br_button(base,0,0,"Redraw",0);
  ap->wedit=br_button(base,0,1,"Edit",0);
  ap->wprint=br_button(base,0,2,"Print",0);
  ap->wclose=br_button(base,0,3,"Close",0);
  ap->wfit=br_button(base,0,4,"Fit",0);
  ap->wrange=br_button(base,0,5,"Range",0);
  ap->wgif=br_button(base,1,0,"GIF",0);
  ap->wmax=make_window(base,10,45,10*DCURXs,DCURYs,1);
  ap->wmin=make_window(base,10,51+DCURYs+color_total,10*DCURXs,DCURYs,1);
  ap->wscale=make_window(base,10+4*DCURXs,48+DCURYs,2*DCURXs,color_total,0);
  ap->wtime=make_window(base,20+10*DCURXs,30,20*DCURXs,DCURYs,0);
  ap->wplot=make_plain_window(base, 20+10*DCURXs,45,width-30-10*DCURXs,height-55,2);
  ap->plotw=width-30-10*DCURXs;
  ap->ploth=height-55;
  ap->alive=1;
  aplot_gc=XCreateGC(display,ap->wplot,valuemask,&values);
}
 
void print_aplot(ap)
     APLOT *ap;
{
  double tlo,thi;
  int status,errflag;
  static char *n[]={"Filename","Top label","Side label","Bottom label", 
	       "Render(-1,0,1,2)"};
   char values[5][MAX_LEN_SBOX];
  int nrows=my_browser.maxrow;
  int row0=ap->nstart;
  int col0=ap->index0;
  int jb;
  if(nrows<=2)return;
  if(ap->plotdef==0||ap->nacross<2||ap->ndown<2)return;
  jb=row0;
  tlo=0.0;
  thi=20.0;
  if(jb>0&&jb<nrows)tlo=my_browser.data[0][jb];
  jb=row0+ap->nskip*(ap->ndown-1);
  if(jb>=nrows)jb=nrows-1;
  if(jb>=0)thi=my_browser.data[0][jb];
  sprintf(values[0],"%s",ap->filename);
  sprintf(values[1],"%s",ap->xtitle);
  sprintf(values[2],"%s",ap->ytitle);
    sprintf(values[3],"%s",ap->bottom);
  sprintf(values[4],"%d",ap->type);
  status=do_string_box(5,5,1,"Print arrayplot",n,values,40);
 if(status!=0){
   strcpy(ap->filename,values[0]);
   strcpy(ap->xtitle,values[1]);
   strcpy(ap->ytitle,values[2]);
   strcpy(ap->bottom,values[3]);
   ap->type=atoi(values[4]);
   if(ap->type<-1||ap->type>2)ap->type=-1;
   errflag=array_print(ap->filename,ap->xtitle,ap->ytitle,ap->bottom,
		       ap->nacross,
		       ap->ndown,col0,row0,ap->nskip,ap->ncskip,
		       nrows,my_browser.maxcol,
		      my_browser.data,ap->zmin,ap->zmax,tlo,thi,ap->type);
   if(errflag==-1)err_msg("Couldn't open file");
 }
}

void apbutton(w)
     Window w;
{
  if(w==aplot.wedit){
    editaplot(&aplot);
  }
  if(w==aplot.wfit){
    fit_aplot();
  }
  if(w==aplot.wrange){
    set_up_aplot_range();

  }
  if(w==aplot.wredraw){
    redraw_aplot(aplot);
  }
  if(w==aplot.wprint){
    print_aplot(&aplot);
  }
  if(w==aplot.wclose){
    destroy_aplot();
  }
  if(w==aplot.wgif){
    gif_aplot();
  }
}

void draw_scale(ap)
     APLOT ap;
{
  int i,y;
  Window w=ap.wscale;
  for(i=0;i<color_total;i++){
  y=color_total-i-1;
  set_color(i+FIRSTCOLOR);
  XDrawLine(display,w,gc_graph,0,y,2*DCURXs,y);
  }
}

void draw_aplot(ap)
     APLOT ap;
{
  if(plot3d_auto_redraw!=1)return;
  redraw_aplot(ap);
}

void edit_aplot()
{
  editaplot(&aplot);
}

void get_root(s,sroot,num)
     char *s,*sroot;
     int *num;
{
  int n=strlen(s);
    int i=n-1,j;

  char me[100];
  *num=0;
  while(1){
   
    if(!isdigit(s[i])){
   
      break;
    }
    i--;
    if(i<0)break;
  }
  if(i<0)strcpy(sroot,s);
  else {
    for(j=0;j<=i;j++)
      sroot[j]=s[j];
    sroot[i+1]=0;
  }
  if(i>=0&&i<n){
    for(j=i+1;j<n;j++)
      me[j-i-1]=s[j];
    me[n-i]=0;
   /* plintf(" i=%d me=%s sroot=%s \n",i,me,sroot); */  
    *num=atoi(me);
  }
}
  
void reset_aplot_axes(ap)
     APLOT ap;
{
  char bob[200];
  char sroot[100];
  int num;
  if(ap.alive==0)return;
  get_root(ap.name,sroot,&num);
  sprintf(bob,"%s%d..%d",sroot,num,num+ap.nacross-1);
  XClearWindow(display,ap.wmax);
  XClearWindow(display,ap.wmin);
  display_aplot(ap.wmax,ap);
  display_aplot(ap.wmin,ap);
  gtitle_text(bob,ap.base);
}

void dump_aplot(fp,f)
     FILE *fp;
     int f;
{
  char bob[256];
  if(f==READEM)
    fgets(bob,255,fp);
  else
    fprintf(fp,"# Array plot stuff\n");
  io_string(aplot.name,11,fp,f);
  io_int(&aplot.nacross ,fp,f,"NCols");
    io_int(&aplot.nstart ,fp,f,"Row 1");
  io_int(&aplot.ndown ,fp,f,"NRows");
  io_int(&aplot.nskip ,fp,f,"RowSkip");
  io_double(&aplot.zmin,fp,f,"Zmin");
    io_double(&aplot.zmax,fp,f,"Zmax");

}

int editaplot(ap)
     APLOT *ap;
{
 int i,status;
 double zmax,zmin;
  char *n[]={"*0Column 1","NCols","Row 1","NRows","RowSkip",
  "Zmin","Zmax","Autoplot(0/1)","ColSkip"};
 char values[9][MAX_LEN_SBOX];
 sprintf(values[0],"%s",ap->name);
 sprintf(values[1],"%d",ap->nacross);
 sprintf(values[2],"%d",ap->nstart);
 sprintf(values[3],"%d",ap->ndown);
 sprintf(values[4],"%d",ap->nskip);
 sprintf(values[5],"%g",ap->zmin);
 sprintf(values[6],"%g",ap->zmax);
 sprintf(values[7],"%d",plot3d_auto_redraw);
sprintf(values[8],"%d",ap->ncskip);
 status=do_string_box(9,9,1,"Edit arrayplot",n,values,40);
 if(status!=0){
   find_variable(values[0],&i);
   if(i>-1){
     ap->index0=i;
     strcpy(ap->name,values[0]);
   }
   else
     {
       err_msg("No such columns");
       ap->plotdef=0;
       return 0;
     }
    zmax=atof(values[6]);
    zmin=atof(values[5]);
    if(zmin<zmax){
      ap->zmin=zmin;
      ap->zmax=zmax;
    }
    ap->nacross=atoi(values[1]);
    ap->nstart=atoi(values[2]);
    ap->ndown=atoi(values[3]);
    ap->nskip=atoi(values[4]);
    plot3d_auto_redraw=atoi(values[7]);
    ap->plotdef=1;
    ap->ncskip=atoi(values[8]);
    if(ap->ncskip<1)
      ap->ncskip=1;
    reset_aplot_axes(*ap);
 }
   return 1;
}
void close_aplot_files()
{
  if(aplot_still==0)
    fclose(ap_fp);
}
void gif_aplot_all(char *filename,int still)
{
  Pixmap xi;
 int x,y;
 unsigned int h,w,bw,d;
 Window root;
 /* FILE *fp; */
 if(still==0) 
   {

     if(aplot_range_count==0){
     
       if((ap_fp=fopen(filename,"w"))==NULL){
	 err_msg("Cannot open file ");
	 return;
       }
     }
     XGetGeometry(display,aplot.wplot,&root,&x,&y,&w,&h,&bw,&d);
     xi=XCreatePixmap(display,RootWindow(display,screen),w,h,
		      DefaultDepth(display,screen));
     XCopyArea(display,aplot.wplot,xi,aplot_gc,0,0,w,h,0,0); 
     
     add_ani_gif(xi,ap_fp,aplot_range_count);
     XFreePixmap(display,xi);
     return;
   } 

 if(still==1){
   if((ap_fp=fopen(filename,"w"))==NULL){
     err_msg("Cannot open file ");
     return;
   }

   /* redraw_aplot(aplot); */
   XGetGeometry(display,aplot.wplot,&root,&x,&y,&w,&h,&bw,&d);
  xi=XCreatePixmap(display,RootWindow(display,screen),w,h,
		  DefaultDepth(display,screen));
 XCopyArea(display,aplot.wplot,xi,aplot_gc,0,0,w,h,0,0); 
 /*  XFlush(display); */
  
  /* screen_to_gif(aplot.wplot,fp); */  
  screen_to_gif(xi,ap_fp);  
  fclose(ap_fp);
  XFreePixmap(display,xi);
 }
}
 
void gif_aplot()
{
  /*char filename[256];*/
  char filename[XPP_MAX_NAME];
  sprintf(filename,"%s.gif",this_file);
  if(!file_selector("GIF plot",filename,"*.gif"))return;
  gif_aplot_all(filename,1);
}


void grab_aplot_screen(ap)
     APLOT ap;
{
 
 Window temp=draw_win;
 draw_win=ap.wplot;
 if(film_clip()==0)
   err_msg("Out of film!");
 draw_win=temp;
}

void redraw_aplot(ap)
     APLOT ap;
{
  int i,j,w=ap.wplot;
  double z,dx,dy,x,y,tlo,thi;
  char bob[100];
  int nrows=my_browser.maxrow,colr,cmax=FIRSTCOLOR+color_total;
  int row0=ap.nstart;
  int col0=ap.index0,delx,dely;
  int ib,jb,ix,iy;
  if(nrows<=2)return;
  if(ap.plotdef==0||ap.nacross<2||ap.ndown<2)return;
  XClearWindow(display,ap.wtime);
  XClearWindow(display,w);
  jb=row0;
  tlo=0.0;
  thi=20.0;
  if(jb>0&&jb<nrows)tlo=my_browser.data[0][jb];
  jb=row0+ap.nskip*(ap.ndown-1);
  if(jb>=nrows)jb=nrows-1;
  if(jb>=0)thi=my_browser.data[0][jb];
  sprintf(bob, " %g < t < %g ",tlo,thi);
   XDrawString(display,ap.wtime,small_gc,0,CURY_OFFs,bob,strlen(bob));
  dx=(double)ap.plotw/(double)(ap.nacross/ap.ncskip);
  dy=(double)ap.ploth/(double)ap.ndown;
  delx=(int)dx+1;
  dely=(int)dy+1;
    for(i=0;i<ap.nacross/ap.ncskip;i++){
      ib=col0+i*ap.ncskip;
      x=dx*i;
      ix=(int)x;

      if(ib>=my_browser.maxcol)return;
      for(j=0;j<ap.ndown;j++){
	jb=row0+ap.nskip*j;
	
	if(jb<nrows&&jb>=0){
	  y=j*dy;
	  iy=(int)y;
/*	  if(j==0)
	  plintf(" ib=%d ix=%d iy=%d \n",ib,ix,iy); */
	  z=(double)color_total*(my_browser.data[ib][jb]-ap.zmin)/(ap.zmax-ap.zmin);
	  colr=(int)z+FIRSTCOLOR;
	  if(colr<FIRSTCOLOR)colr=FIRSTCOLOR;
	  if(colr>cmax)colr=cmax;
	  set_acolor(colr);
	  XFillRectangle(display,w,aplot_gc,ix,iy,delx,dely);
	}
      }
    }
    XFlush(display);
}
void tag_aplot(char *bob)
{
  set_color(0);
  XDrawString(display,aplot.wplot,small_gc,0,CURY_OFFs,bob,strlen(bob));
}

void set_acolor(col)
int col;
{
 if(col<0)XSetForeground(display,aplot_gc,GrBack);
 if(col==0)XSetForeground(display,aplot_gc,GrFore);
 else{

   if(COLOR)XSetForeground(display,aplot_gc,ColorMap(col));
   else XSetForeground(display,aplot_gc,GrFore);
}

}


void display_aplot(w,ap)
     APLOT ap;
     Window w;
{
  char bob[200];
  
  if(w==ap.wplot){
    draw_aplot(ap);
    return;
  }
  if(w==ap.wscale){
    draw_scale(ap);
    return;
  }
  if(w==ap.wmin){
    sprintf(bob,"%g",ap.zmin);
      XDrawString(display,w,small_gc,0,CURY_OFFs,bob,strlen(bob));
		return;
  }
 if(w==ap.wmax){
    sprintf(bob,"%g",ap.zmax);
      XDrawString(display,w,small_gc,0,CURY_OFFs,bob,strlen(bob));
		return;
  }
 if(w==ap.wedit){
   XDrawString(display,w,small_gc,0,CURY_OFFs,"Edit",4);
   return;
 }
 if(w==ap.wgif) {
   XDrawString(display,w,small_gc,0,CURY_OFFs,"GIF",3);
   return;
 }
if(w==ap.wredraw){
   XDrawString(display,w,small_gc,0,CURY_OFFs,"Redraw",6);
   return;
 }
if(w==ap.wfit){
   XDrawString(display,w,small_gc,0,CURY_OFFs,"Fit",3);
   return;
 }
if(w==ap.wrange){
   XDrawString(display,w,small_gc,0,CURY_OFFs,"Range",5);
   return;
 }
  if(w==ap.wprint){
   XDrawString(display,w,small_gc,0,CURY_OFFs,"Print",5);
   return;
 }

  /*if(w==ap.wkill){
   XDrawString(display,w,small_gc,0,CURY_OFFs,"Kill",4);
   return;
 }
  */
 
if(w==ap.wclose){
   XDrawString(display,w,small_gc,0,CURY_OFFs,"Close",5);
   return;
 }
}





















