#include "menu.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/cursorfont.h>
#include <X11/Xatom.h>
#include <X11/Xlib.h>
#include <X11/Xos.h>
#include <X11/Xutil.h>

#include "ggets.h"
#include "main.h"
#include "many_pops.h"
#include "menudrive.h"
#include "menus.h"
#include "pop_list.h"


/* --- Forward declarations --- */
static void add_menu(Window base, int j, int n, char **names, char *key, char **hint);
static void show_menu(int j);
static void unshow_menu(int j);


/* --- Data --- */
static MENUDEF my_menus[3];
int help_menu;


/* --- Functions --- */
void create_the_menus(Window base) {
	char key[30];
	strcpy(key,"icndwakgufpemtsvxr3b");
	add_menu(base,MAIN_MENU,MAIN_ENTRIES,main_menu,key,main_hint);
	strcpy(key,"tsrdniobmechpukva");
	key[17]=27;
	key[18]=0;
	add_menu(base,NUM_MENU,NUM_ENTRIES,num_menu,key,num_hint);
	/* CLONE */
	strcpy(key,"pwracesbhqtiglxu");
	add_menu(base,FILE_MENU,FILE_ENTRIES,fileon_menu,key,file_hint);
	help_menu=-1;
}


void draw_help(void) {
	int i,j=help_menu,n;
	if(j<0) {
		return;
	}
	if(my_menus[j].visible==0) {
		return;
	}
	n=my_menus[j].n;
	menu_expose(my_menus[j].title);
	for(i=0;i<n;i++) {
		menu_expose(my_menus[j].w[i]);
	}
}


void help(void) {
	unshow_menu(help_menu);
	show_menu(MAIN_MENU);
}


void help_file(void) {
	if(tfBell) {
		my_menus[FILE_MENU].names=fileon_menu;
	} else {
		my_menus[FILE_MENU].names=fileoff_menu;
	}
	unshow_menu(help_menu);
	show_menu(FILE_MENU);
}


void help_num(void) {
	unshow_menu(help_menu);
	show_menu(NUM_MENU);
}


void menu_button(Window win) {
	int i,n,j=help_menu;
	if(j<0) {
		return;
	}
	if(my_menus[j].visible==0) {
		return;
	}
	n=my_menus[j].n;
	for(i=0;i<n;i++) {
		if(win==my_menus[j].w[i]) {
			XSetWindowBorderWidth(display,win,0);
			commander(my_menus[j].key[i]);
			return;
		}
	}
}


void menu_crossing(Window win, int yn) {
	int i,n,j=help_menu;
	char **z;
	if(j<0) {
		return;
	}
	if(my_menus[j].visible==0) {
		return;
	}
	n=my_menus[j].n;
	z=my_menus[j].hints;
	for(i=0;i<n;i++) {
		if(win==my_menus[j].w[i]) {
			XSetWindowBorderWidth(display,win,yn);
			if(yn&&TipsFlag) {
				bottom_msg(0,z[i]);
			}
			return;
		}
	}
}


void menu_expose(Window win) {
	int i,n,j=help_menu;
	char **z;
	if(j<0) {
		return;
	}
	if(my_menus[j].visible==0) {
		return;
	}
	n=my_menus[j].n;
	z=my_menus[j].names;
	if(win==my_menus[j].title) {
		set_fore();
		bar(0,0,16*DCURX,DCURY,win);
		set_back();
		XDrawString(display,win,gc,DCURX/2+5,CURY_OFF,z[0],strlen(z[0]));
		set_fore();
		return;
	}
	for(i=0;i<n;i++) {
		if(win==my_menus[j].w[i]) {
			BaseCol();
			XDrawString(display,win,gc,5,CURY_OFF,z[i+1],strlen(z[i+1]));
			return;
		}
	}
}


/* --- Static functions --- */
static void add_menu(Window base, int j, int n, char **names, char *key, char **hint) {
	Window w;
	int i;
	Cursor cursor;
	cursor=XCreateFontCursor(display,XC_hand2);
	w=make_plain_unmapped_window(base,0,DCURYs+DCURYb+10,16*DCURX,21*(DCURY+2)-3,1);
	my_menus[j].base=w;
	XDefineCursor(display,w,cursor);
	my_menus[j].names=names;
	my_menus[j].n=n;
	my_menus[j].hints=hint;
	strcpy(my_menus[j].key,key);
	my_menus[j].title=make_unmapped_window(w,0,0,16*DCURX,DCURY,1);
	for(i=0;i<n;i++) {
		my_menus[j].w[i]=make_unmapped_window(w,0,(i+1)*(DCURY+2),16*DCURX,DCURY,0);
	}
	my_menus[j].visible=0;
	XMapRaised(display,my_menus[j].base);
	XMapSubwindows(display,my_menus[j].base);
}


static void show_menu(int j) {
	XRaiseWindow(display,my_menus[j].base);
	my_menus[j].visible=1;
	help_menu=j;
}


static void unshow_menu(int j) {
	if(j<0) {
		return;
	}
	my_menus[j].visible=0;
}
