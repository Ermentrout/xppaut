#ifndef _aniparse_h_
#define _aniparse_h_


#include <X11/Xlib.h>
#include <stdio.h>

/**************  New stuff for the Grabber ***************************/
#define MAX_GEVENTS 20  /* maximum variables you can change per grabbable */
#define MAX_ANI_GRAB 50   /* max grabbable objects  */


typedef struct {  /* tasks have the form {name1=formula1;name2=formula2;...} */
 
  double vrhs[MAX_GEVENTS];
  char lhsname[MAX_GEVENTS][11]; 
  int lhsivar[MAX_GEVENTS];
  int *comrhs[MAX_GEVENTS];
  int runnow;
  int n; /* number of tasks <= MAX_GEVENTS */
}GRAB_TASK;


typedef struct {
  int ok;
  double zx,zy,tol;
  int *x,*y;
  GRAB_TASK start,end;
} ANI_GRAB;

/***************  End of grabber stuff  in header **************/

typedef struct {
  int flag;
 int skip;
  char root[100];
 char filter[256];
 int aviflag,filflag;
} MPEG_SAVE;


typedef struct {
  int n;
  int *x,*y,*col;
  int i;
} Comet;

typedef struct {
  Comet c;
  int type, flag;
  int *col,*x1,*y1,*x2,*y2,*who;
  double zcol,zx1,zy1,zx2,zy2,zrad,zval;
  int zthick,tfont,tsize,tcolor;  
} ANI_COM;

void new_vcr(void);
void create_vcr(char *name);
void ani_border(Window w, int i);
void destroy_vcr(void);
void do_ani_events(XEvent ev);
void ani_motion_stuff(Window w, int x, int y);
double get_current_time(void);
void update_ani_motion_stuff(int x, int y);
void ani_buttonx(XEvent ev, int flag);
void ani_button(Window w);
void ani_create_mpeg(void);
void ani_expose(Window w);
void ani_resize(int x, int y);
void ani_newskip(void);
void check_on_the_fly(void);
void on_the_fly(int task);
void ani_frame(int task);
void set_to_init_data(void);
void set_from_init_data(void);
void ani_flip1(int n);
void ani_flip(void);
void ani_disk_warn(void);
int getppmbits(Window window, int *wid, int *hgt, unsigned char *out);
int writeframe(char *filename, Window window, int wid, int hgt);
void ani_zero(void);
void get_ani_file(char *fname);
int ani_new_file(char *filename);
int load_ani_file(FILE *fp);
int parse_ani_string(char *s, FILE *fp);
void set_ani_dimension(char *x1, char *y1, char *x2, char *y2);
int add_ani_com(int type, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
void init_ani_stuff(void);
void free_ani(void);
int chk_ani_color(char *s, int *index);
int add_ani_expr(char *x, int *c);
int add_ani_rline(ANI_COM *a, char *x1, char *y1, char *col, char *thick);
void reset_comets(void);
void roll_comet(ANI_COM *a, int xn, int yn, int col);
int add_ani_comet(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_line(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_null(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *who);
int add_ani_rect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_frect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_ellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_fellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_circle(ANI_COM *a, char *x1, char *y1, char *x2, char *col, char *thick);
int add_ani_fcircle(ANI_COM *a, char *x1, char *y1, char *x2, char *col, char *thick);
int add_ani_text(ANI_COM *a, char *x1, char *y1, char *y2);
int add_ani_vtext(ANI_COM *a, char *x1, char *y1, char *x2, char *y2);
int add_ani_settext(ANI_COM *a, char *x1, char *y1, char *col);
void render_ani(void);
void set_ani_perm(void);
void eval_ani_color(int j);
void eval_ani_com(int j);
void set_ani_thick(int t);
void set_ani_font_stuff(int size, int font, int color);
void set_ani_col(int j);
void xset_ani_col(int icol);
void ani_rad2scale(double rx, double ry, int *ix, int *iy);
void ani_radscale(double rad, int *ix, int *iy);
void ani_ij_to_xy(int ix, int iy, double *x, double *y);
void ani_xyscale(double x, double y, int *ix, int *iy);
void draw_ani_comet(int j);
void draw_ani_null(int j, int id);
void draw_ani_line(int j);
void draw_ani_rline(int j);
void draw_ani_circ(int j);
void draw_ani_fcirc(int j);
void draw_ani_rect(int j);
void draw_ani_frect(int j);
void draw_ani_ellip(int j);
void draw_ani_fellip(int j);
void draw_ani_text(int j);
void draw_ani_vtext(int j);
void tst_pix_draw(void);
void read_ani_line(FILE *fp, char *s);
void de_space(char *s);
int add_grab_command(char *xs, char *ys, char *ts, FILE *fp);
void info_grab_stuff(void);
int ani_grab_tasks(char *line, int igrab, int which);
int run_now_grab(void);
int search_for_grab(double x, double y);
void do_grab_tasks(int which);
int add_grab_task(char *lhs, char *rhs, int igrab, int which);
void draw_grab_points(void);
void free_grabber(void);
int check_ani_pause(XEvent ev);
void do_ani_slider_motion(Window w,int x);
void draw_ani_slider(Window w,int x);
void redraw_ani_slider(void);
#endif
