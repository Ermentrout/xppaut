#ifndef XPPAUT_EDIT_RHS_H
#define XPPAUT_EDIT_RHS_H


#include <X11/Xlib.h>
#include "xpplim.h"
#include <stdio.h>

/* --- Macros --- */
#define NEQMAXFOREDIT 20
#define MAXARG 20
#define MAX_N_EBOX MAXODE
#define MAX_LEN_EBOX 86
#define RESET_ALL 4
#define MAXUFUN 50

/*  This is a edit box widget which handles a list of
	editable strings
 */

/* --- Types --- */
typedef struct {
		Window base,ok,cancel,reset;
		Window win[MAX_N_EBOX];
		char name[MAX_N_EBOX][MAX_LEN_EBOX],
			 value[MAX_N_EBOX][MAX_LEN_EBOX],
			 rval[MAX_N_EBOX][MAX_LEN_EBOX];
		int n,hot;
} EDIT_BOX;

/* --- Functions --- */
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

#endif /* XPPAUT_EDIT_RHS_H */
