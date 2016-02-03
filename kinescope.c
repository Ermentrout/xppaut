#include "kinescope.h"


#include "scrngif.h"
#include "pop_list.h"
#include "aniparse.h"
#include "browse.h"
#include "load_eqn.h"

#include "ggets.h"
#include <stdlib.h> 
/*    Kinescope for X  windows       */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
/* #include <X11/bitmaps/icon> */
#include <stdio.h>
#include <sys/time.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include "mykeydef.h"
extern Display *display;
extern Window draw_win,main_win,info_pop;
extern int DCURY;
extern GC gc_graph;
#define MAXFILM 250
int ks_ncycle=1;
int ks_speed=50;
extern char *info_message,*kin_hint[];
extern int screen;
int mov_ind;
typedef struct {
		unsigned int h,w;
		Pixmap xi;
		} MOVIE;

MOVIE movie[MAXFILM];

void do_movie_com(int c)

{



/*  XDestroyWindow(display,temp); 
  draw_help();
  XFlush(display); */ 
 switch(c){

	     case 0: if(film_clip()==0)
	 respond_box("Okay","Out of film!");
			break;
	      case 1: reset_film();
	                break;
	      case 2: play_back();
		        break;
	      case 3: auto_play();
			 break;
              case 4: save_kine();
		        break;
              case 5: make_anigif();
		        break;
 case 6: /* test_keys(); */
   break;
	    }
      
  }



void reset_film()
{
 int i;
 if(mov_ind==0)return;
 for(i=0;i<mov_ind;i++)XFreePixmap(display,movie[i].xi);
 mov_ind=0;
 }

 

int film_clip()
{
 int x,y;
 unsigned int h,w,bw,d;
 Window root;
 if(mov_ind>=MAXFILM)return(0);
 XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
 movie[mov_ind].h=h;
 movie[mov_ind].w=w;
 movie[mov_ind].xi=XCreatePixmap(display,RootWindow(display,screen),w,h,
		  DefaultDepth(display,screen));
 XCopyArea(display,draw_win,movie[mov_ind].xi,gc_graph,0,0,w,h,0,0);
 mov_ind++;
 return 1;
}

int show_frame(int i, int h, int w)
{
  if(h<movie[i].h||w<movie[i].w){
    too_small();
    return 1;
  }
  XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
  XFlush(display);

  return 0;
}

void play_back()
{
  int x,y;
  unsigned int h,w,bw,d;

  Window root;
  XEvent ev;
  int i=0;
  XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
  if(mov_ind==0)return;
  if(h<movie[i].h||w<movie[i].w){
    too_small();
    return;
  }
  
  XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
  XFlush(display);
  while(1){
    XNextEvent(display,&ev);
    switch(ev.type)
      {
      case ButtonPress:
	i++;
	if(i>=mov_ind)i=0;
	if(show_frame(i,h,w))return;
	break;
      case KeyPress:
	switch(get_key_press(&ev)){
	  case ESC:
	    return;
	  case RIGHT:
	    i++;
	    if(i>=mov_ind)i=0;
	    if(show_frame(i,h,w))return;
	    break;
	  case  LEFT:
	    i--;
	    if(i<0)i=mov_ind-1;
	    if(show_frame(i,h,w))return;
	    break;
	  case  HOME:
	    i=0;
	    if(show_frame(i,h,w))return;
	    break;
	  case END:
	    i=mov_ind-1;
	    if(show_frame(i,h,w))return;
	    break;
	  }
	  
      }
  }
}
	       		

     
void save_kine()
{
 char base[128];
 int fmat=2;
 sprintf(base,"frame");
 /* #ifdef NOGIF
#else
new_int("format:1-ppm,2-gif",&fmat);
#endif
 */
 new_string("Base file name",base);
 if(strlen(base)>0)
   save_movie(base,fmat);
   
}
  


     
void make_anigif()
{
  int i=0;
  int x,y;
  FILE *fp;
  Window root;
  unsigned int h,w,bw,d;
  XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
  if(mov_ind==0)return;
  if(h<movie[i].h||w<movie[i].w){
    too_small();
    return;
    
  }
 h=movie[0].h;
 w=movie[0].w;
 for(i=0;i<mov_ind;i++){
   if((movie[i].h!=h)||(movie[i].w!=w)){
     err_msg("All clips must be same size");
     return;
   }
 }
 fp=fopen("anim.gif","wb");
 set_global_map(1);
  for(i=0;i<mov_ind;i++){
    XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
    XFlush(display);
    /* add_ani_gif(draw_win,fp,i); */
        add_ani_gif(movie[i].xi,fp,i);
  }
  
 end_ani_gif(fp);
   fclose(fp);
   set_global_map(0);

}


     
void save_movie(basename,fmat)
     char *basename;
     int fmat;
{
  /*char file[256];
  */
  char file[XPP_MAX_NAME];
  int i=0;
  int x,y;
  FILE *fp;
  Window root;
  Pixmap xi;
  unsigned int h,w,bw,d;
  XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
  if(mov_ind==0)return;
  if(h<movie[i].h||w<movie[i].w){
    too_small();
    return;
    
  }
  
  for(i=0;i<mov_ind;i++){
    if(fmat==1)
      sprintf(file,"%s_%d.ppm",basename,i);
    else
      sprintf(file,"%s_%d.gif",basename,i);
    XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
    XFlush(display);
   if(fmat==1)
     writeframe(file,draw_win,w,h); 
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
 

     
void  auto_play()
{
 int x,y;
 unsigned int h,w,bw,d,key;
 Window root;
 
  int dt=20;
  int smax=500;
 XEvent ev;
 int i=0,cycle=0;

 
 new_int("Number of cycles",&ks_ncycle);
 new_int("Msec between frames",&ks_speed);
 if(ks_speed<0)ks_speed=0;
 if(ks_ncycle<=0)return;
 XGetGeometry(display,draw_win,&root,&x,&y,&w,&h,&bw,&d);
 if(mov_ind==0)return;
 if(h<movie[i].h||w<movie[i].w){
	too_small();
	return;

	}

		  XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
   XFlush(display);


 while(1)
{
	
	  
	  	  
	/* check for events    */
	  if(XPending(display)>0)
	  {
	   
          XNextEvent(display,&ev);
          switch(ev.type){
			  case ButtonPress:
		
			    return;
			   break;
		
				
			case  KeyPress:
			     key=get_key_press(&ev);
			      if(key==27)return;
			      if(key==','){
				ks_speed-=dt;
			    if(ks_speed<dt)ks_speed=dt;
			      }
			      if(key=='.'){
				 ks_speed+=dt;
			       if(ks_speed>smax)ks_speed=smax;
			      }
			   
		         break;
                        }

           }  /* done checking  now increment pix   */

		waitasec(ks_speed);
		i++;
		if(i>=mov_ind){cycle++; i=0;}
 		if(h<movie[i].h||w<movie[i].w){
					too_small();
					return;
					}
		XCopyArea(display,movie[i].xi,draw_win,gc_graph,0,0,w,h,0,0);
		XFlush(display);
                if(cycle>=ks_ncycle)return;
 

	}  /*  Big loop   */
  }	  


       

void too_small()
 {
  respond_box("Okay","Window too small for film!");
  }

 
 
  
 
 
 
 

		
