#ifndef _pop_list_h
#define _pop_list_h



#include "phsplan.h"
#include <stdlib.h> 
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <stdio.h>
#include "xpplim.h"
#include "math.h"



#define MAX_N_SBOX 20
#define MAX_LEN_SBOX 25 


#define FORGET_ALL 0
#define DONE_ALL 2
#define FORGET_THIS 3
#define DONE_THIS 1





#define EV_MASK (ButtonPressMask 	|\
		KeyPressMask		|\
		ExposureMask		|\
		StructureNotifyMask)	

#define BUT_MASK (ButtonPressMask 	|\
		KeyPressMask		|\
		ExposureMask		|\
		StructureNotifyMask	|\
		EnterWindowMask		|\
		LeaveWindowMask)	



	
	


extern Display *display;
extern int DisplayWidth,DisplayHeight;
extern int screen;
extern Atom deleteWindowAtom;
extern Window main_win,info_pop,draw_win;
extern int DCURY,DCURX,CURY_OFF,DCURXs,DCURYs,CURY_OFFs,xor_flag;
extern GC gc,small_gc;
extern unsigned int MyBackColor,MyForeColor;

extern int TipsFlag;
char *get_next(),*get_first();
extern char UserBlack[8];
extern char UserWhite[8];
extern int UserGradients;


Window make_window();
Window make_plain_window();

/*  This is a string box widget which handles a list of 
	editable strings  
 */

typedef struct {
		Window base,ok,cancel;
		Window win[MAX_N_SBOX];
		char name[MAX_N_SBOX][MAX_LEN_SBOX],
		     value[MAX_N_SBOX][MAX_LEN_SBOX];
		int n,hot;
                int hgt,wid;
                int hh[MAX_N_SBOX];
		} STRING_BOX;


typedef struct {
               char **list;
               int n;
}  SCRBOX_LIST;

extern int NUPAR,NEQ,NODE,NMarkov;
extern char  upar_names[MAXPAR][11],uvar_names[MAXODE][12];
extern  char *color_names[];
SCRBOX_LIST scrbox_list[10];


/*  This is a new improved pop_up widget */
typedef struct {
		Window base,tit;
		Window *w;
		char *title;
		char **entries;
                char **hints;
		int n,max;
		char *key;
		int hot;
		} POP_UP;
		
   
typedef struct {
               Window base,slide,close,text;
               int i0;
               int exist,len,nlines;
               char **list;
               } TEXTWIN;

typedef struct {
              Window base,slide;
              Window *w;
              int nw,nent,i0;
              int len,exist;
              char **list;
              } SCROLLBOX;

TEXTWIN mytext;
#define SB_PLOTTABLE 0
#define SB_VARIABLE 1
#define SB_PARAMETER 2
#define SB_PARVAR 3
#define SB_COLOR 4
#define SB_MARKER 5
#define SB_METHOD 6





void set_window_title(Window win, char *string);
void make_scrbox_lists(void);
int get_x_coord_win(Window win);
void destroy_scroll_box(SCROLLBOX *sb);
void create_scroll_box(Window root, int x0, int y0, int nent, int nw, char **list, SCROLLBOX *sb);
void expose_scroll_box(Window w, SCROLLBOX sb);
void redraw_scroll_box(SCROLLBOX sb);
void crossing_scroll_box(Window w, int c, SCROLLBOX sb);
int scroll_box_motion(XEvent ev, SCROLLBOX *sb);
int select_scroll_item(Window w, SCROLLBOX sb);
void scroll_popup(STRING_BOX *sb, SCROLLBOX *scrb);
int do_string_box(int n, int row, int col, char *title, char **names, char values[][25], int maxchar);
void expose_sbox(STRING_BOX sb, Window w, int pos, int col);
void do_hilite_text(char *name, char *value, int flag, Window w, int pos, int col);
void reset_hot(int inew, STRING_BOX *sb);
void new_editable(STRING_BOX *sb, int inew, int *pos, int *col, int *done, Window *w);
void set_sbox_item(STRING_BOX *sb, int item);
int s_box_event_loop(STRING_BOX *sb, int *pos, int *col, SCROLLBOX *scrb);
void make_sbox_windows(STRING_BOX *sb, int row, int col, char *title, int maxchar);
Window make_fancy_window(Window root, int x, int y, int width, int height, int bw, int fc, int bc);
Window make_unmapped_window(Window root, int x, int y, int width, int height, int bw);
Window make_plain_unmapped_window(Window root, int x, int y, int width, int height, int bw);
Window make_window(Window root, int x, int y, int width, int height, int bw);
Window make_plain_window(Window root, int x, int y, int width, int height, int bw);
void expose_resp_box(char *button, char *message, Window wb, Window wm, Window w);
void respond_box(char *button, char *message);
void message_box(Window *w, int x, int y, char *message);
void expose_choice(char *choice1, char *choice2, char *msg, Window c1, Window c2, Window wm, Window w);
int two_choice(char *choice1, char *choice2, char *string, char *key, int x, int y, Window w,char *title);
int yes_no_box(void);
int pop_up_list(Window *root, char *title, char **list, char *key, int n, int max, int def, int x, int y, char **hints, Window hwin, char *httxt);
void draw_pop_up(POP_UP p, Window w);
Window make_unmapped_icon_window(Window root,int x,int y,int width,int height,int bw,int icx,int icy,unsigned char* icdata);
Window make_icon_window(Window root,int x,int y,int width,int height,int bw,int icx,int icy,unsigned char* icdata);

#endif
