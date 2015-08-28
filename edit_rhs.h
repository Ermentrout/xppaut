#ifndef _edit_rhs_h_
#define _edit_rhs_h_


#include <X11/Xlib.h>
#include "xpplim.h"
#include <stdio.h>


#define NEQMAXFOREDIT 20
#define MAXARG 20
#define MAX_N_EBOX MAXODE
#define MAX_LEN_EBOX 86
#define FORGET_ALL 0
#define DONE_ALL 2
#define FORGET_THIS 3
#define DONE_THIS 1
#define RESET_ALL 4

#define MAXUFUN 50



/*typedef struct {
  int narg;
  char args[MAXARG][11];
} UFUN_ARG;
*/


/*  This is a edit box widget which handles a list of 
	editable strings  
 */

typedef struct {
		Window base,ok,cancel,reset;
		Window win[MAX_N_EBOX];
		char name[MAX_N_EBOX][MAX_LEN_EBOX],
		     value[MAX_N_EBOX][MAX_LEN_EBOX],
		     rval[MAX_N_EBOX][MAX_LEN_EBOX];
		int n,hot;
		} EDIT_BOX;


void reset_ebox(EDIT_BOX *sb, int *pos, int *col);
int do_edit_box(int n, char *title, char **names, char **values);
void expose_ebox(EDIT_BOX *sb, Window w, int pos, int col);
void ereset_hot(int inew, EDIT_BOX *sb);
void enew_editable(EDIT_BOX *sb, int inew, int *pos, int *col, int *done, Window *w);
int e_box_event_loop(EDIT_BOX *sb, int *pos, int *col);
void make_ebox_windows(EDIT_BOX *sb, char *title);
void edit_menu(void);
void edit_rhs(void);
void user_fun_info(FILE *fp);
void edit_functions(void);
int save_as(void);

#endif
