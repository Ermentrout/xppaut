#ifndef _graphics_h
#define _graphics_h


void get_scale(double *x1, double *y1, double *x2, double *y2);
void set_scale(double x1, double y1, double x2, double y2);
void get_draw_area(void);
void change_current_linestyle(int new, int *old);
void set_normal_scale(void);
void point(int x, int y);
void line(int x1, int y1, int x2, int y2);
void bead(int x1, int y1);
void frect(int x1, int y1, int w, int h);
void put_text(int x, int y, char *str);
void init_x11(void);
void init_ps(void);
void init_svg(void);
void point_x11(int xp, int yp);
void set_linestyle(int ls);
void set_line_style_x11(int ls);
void bead_x11(int x, int y);
void rect_x11(int x, int y, int w, int h);
void line_x11(int xp1, int yp1, int xp2, int yp2);
void put_text_x11(int x, int y, char *str);
void special_put_text_x11(int x, int y, char *str, int size);
void fancy_put_text_x11(int x, int y, char *str, int size, int font);
void scale_dxdy(double x, double y, double *i, double *j);
void scale_to_screen(double x, double y, int *i, int *j);
void scale_to_real(int i, int j, float *x, float *y);
void init_all_graph(void);
void set_extra_graphs(void);
void reset_graph(void);
void get_graph(void);
void init_graph(int i);
void copy_graph(int i, int l);
void make_rot(double theta, double phi);
void scale3d(double x, double y, double z, float *xp, float *yp, float *zp);
double proj3d(double theta, double phi, double x, double y, double z, int in);
int threedproj(double x2p, double y2p, double z2p, float *xp, float *yp);
void text3d(double x, double y, double z, char *s);
void text_3d(double x, double y, double z, char *s);
int threed_proj(double x, double y, double z, float *xp, float *yp);
void point_3d(double x, double y, double z);
void line3dn(double xs1, double ys1, double zs1, double xsp1, double ysp1, double zsp1);
void line3d(double x01, double y01, double z01, double x02, double y02, double z02);
void line_3d(double x, double y, double z, double xp, double yp, double zp);
void pers_line(double x, double y, double z, double xp, double yp, double zp);
void rot_3dvec(double x, double y, double z, float *xp, float *yp, float *zp);
void point_abs(double x1, double y1);
void line_nabs(double x1_out, double y1_out, double x2_out, double y2_out);
void bead_abs(double x1, double y1);
void frect_abs(double x1, double y1, double w, double h);
void line_abs(double x1, double y1, double x2, double y2);
void text_abs(double x, double y, char *text);
void fillintext(char *old, char *new);
void fancy_text_abs(double x, double y, char *old, int size, int font);
int clip3d(double x1, double y1, double z1, double x2, double y2, double z2, float *x1p, float *y1p, float *z1p, float *x2p, float *y2p, float *z2p);
int clip(double x1, double x2, double y1, double y2, float *x1_out, float *y1_out, float *x2_out, float *y2_out);
void eq_symb(double *x, int type);
void draw_symbol(double x, double y, double size, int my_symb);
void reset_all_line_type();

#endif
