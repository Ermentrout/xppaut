#ifndef _browse_h_
#define _browse_h_

#define BMAXCOL 20

#include <X11/Xlib.h>
#include <stdio.h>

typedef struct {
		Window base,upper;
		Window find,up,down,pgup,pgdn,home,end,left,right;
		Window first,last,restore,write,get,close;
		Window load,repl,unrepl,table,addcol,delcol;
                Window main;
                Window label[BMAXCOL];
                Window time;
                Window hint;
		char hinttxt[256];
		int dataflag,xflag;
		int col0,row0,ncol,nrow;
		int maxrow,maxcol;
                float **data;
		int istart,iend;
                } BROWSER;

/*extern BROWSER my_browser;
*/

float **get_browser_data(void);
void set_browser_data(float **data, int col0);
float *get_data_col(int c);
int gettimenow(void);
void waitasec(int msec);
int get_maxrow_browser(void);
void write_mybrowser_data(FILE *fp);
void write_browser_data(FILE *fp, BROWSER *b);
int check_for_stor(float **data);
void del_stor_col(char *var, BROWSER *b);
void data_del_col(BROWSER *b);
void data_add_col(BROWSER *b);
int add_stor_col(char *name, char *formula, BROWSER *b);
void chk_seq(char *f, int *seq, double *a1, double *a2);
void replace_column(char *var, char *form, float **dat, int n);
void wipe_rep(void);
void unreplace_column(void);
void make_d_table(double xlo, double xhi, int col, char *filename, BROWSER b);
void find_value(int col, double val, int *row, BROWSER b);
void find_variable(char *s, int *col);
void browse_but_on(BROWSER *b, int i, Window w, int yn);
void enter_browser(XEvent ev, BROWSER *b, int yn);
void display_browser(Window w, BROWSER b);
void redraw_browser(BROWSER b);
void new_browse_dat(float **new_dat, int dat_len);
void refresh_browser(int length);
void reset_browser(void);
void draw_data(BROWSER b);
void init_browser(void);
void kill_browser(BROWSER *b);
void make_new_browser(void);
Window br_button(Window root, int row, int col, char *name, int iflag);
Window br_button_data(Window root, int row, int col, char *name, int iflag);
void make_browser(BROWSER *b, char *wname, char *iname, int row, int col);
void expose_my_browser(XEvent ev);
void enter_my_browser(XEvent ev, int yn);
void my_browse_button(XEvent ev);
void my_browse_keypress(XEvent ev, int *used);
void resize_my_browser(Window win);
void expose_browser(XEvent ev, BROWSER b);
void resize_browser(Window win, BROWSER *b);
void browse_button(XEvent ev, BROWSER *b);
void browse_keypress(XEvent ev, int *used, BROWSER *b);
void data_up(BROWSER *b);
void data_down(BROWSER *b);
void data_pgup(BROWSER *b);
void data_pgdn(BROWSER *b);
void data_home(BROWSER *b);
void data_end(BROWSER *b);
void get_data_xyz(float *x, float *y, float *z, int i1, int i2, int i3, int off);
void data_get(BROWSER *b);
void data_replace(BROWSER *b);
void data_unreplace(BROWSER *b);
void data_table(BROWSER *b);
void data_find(BROWSER *b);
void open_write_file(FILE **fp, char *fil, int *ok);
void data_read(BROWSER *b);
void data_write(BROWSER *b);
void data_left(BROWSER *b);
void data_right(BROWSER *b);
void data_first(BROWSER *b);
void data_last(BROWSER *b);
void data_restore(BROWSER *b);
void get_col_list(char *s, int *cl, int *n);

#endif











