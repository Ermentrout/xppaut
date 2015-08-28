#ifndef _aniparse_avi_h
#define _aniparse_avi_h


typedef struct {
  int nframe,wid,hgt,fps;
  unsigned char *image;
  int cur_frame;
  int task;
} AVI_INFO;

AVI_INFO avi_info;
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
  int *col,*x1,*y1,*x2,*y2;
  double zcol,zx1,zy1,zx2,zy2,zrad,zval;
  int zthick,tfont,tsize,tcolor;  
} ANI_COM;



int new_vcr(void);
int create_vcr(char *name);
int ani_border(Window w, int i);
int do_ani_events(XEvent ev);
int ani_button(Window w);
int ani_create_mpeg(void);
int ani_expose(Window w);
int ani_resize(int x, int y);
int ani_newskip(void);
int ani_flip1(int n);
int ani_flip(void);
int ani_disk_warn(void);
int getppmbits(Window window, int *wid, int *hgt, unsigned char *out);
int writeframe(char *filename, Window window, int wid, int hgt);
int ani_zero(void);
int get_ani_file(void);
int ani_new_file(char *filename);
int load_ani_file(FILE *fp);
int parse_ani_string(char *s);
int set_ani_dimension(char *x1, char *y1, char *x2, char *y2);
int add_ani_com(int type, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int init_ani_stuff(void);
int free_ani(void);
int chk_ani_color(char *s, int *index);
int add_ani_expr(char *x, int *c);
int add_ani_rline(ANI_COM *a, char *x1, char *y1, char *col, char *thick);
int reset_comets(void);
int roll_comet(ANI_COM *a, int xn, int yn, int col);
int add_ani_comet(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_line(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_rect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_frect(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_ellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_fellip(ANI_COM *a, char *x1, char *y1, char *x2, char *y2, char *col, char *thick);
int add_ani_circle(ANI_COM *a, char *x1, char *y1, char *x2, char *col, char *thick);
int add_ani_fcircle(ANI_COM *a, char *x1, char *y1, char *x2, char *col, char *thick);
int add_ani_text(ANI_COM *a, char *x1, char *y1, char *y2);
int add_ani_vtext(ANI_COM *a, char *x1, char *y1, char *x2, char *y2);
int add_ani_settext(ANI_COM *a, char *x1, char *y1, char *col);
int render_ani(void);
int set_ani_perm(void);
int eval_ani_color(int j);
int eval_ani_com(int j);
int set_ani_thick(int t);
int set_ani_font_stuff(int size, int font, int color);
int set_ani_col(int j);
int xset_ani_col(int icol);
int ani_rad2scale(double rx, double ry, int *ix, int *iy);
int ani_radscale(double rad, int *ix, int *iy);
int ani_xyscale(double x, double y, int *ix, int *iy);
int draw_ani_comet(int j);
int draw_ani_line(int j);
int draw_ani_rline(int j);
int draw_ani_circ(int j);
int draw_ani_fcirc(int j);
int draw_ani_rect(int j);
int draw_ani_frect(int j);
int draw_ani_ellip(int j);
int draw_ani_fellip(int j);
int draw_ani_text(int j);
int draw_ani_vtext(int j);
int tst_pix_draw(void);
int read_ani_line(FILE *fp, char *s);
int de_space(char *s);


#endif
