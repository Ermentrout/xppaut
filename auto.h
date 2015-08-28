#ifndef _auto_h_
#define _auto_h_


#include <X11/Xlib.h>
#include <stdio.h>


int get_auto_str(char *xlabel, char *ylabel);
int draw_ps_axes(void);
int draw_bif_axes(void);
int byeauto_(int *nt, int *iflag);
int IXVal(double x);
int IYVal(double y);
int Circle(int x, int y, int r);
int XORCross(int x, int y);
int FillCircle(int x, int y, int r);
int LineWidth(int wid);
int renamef(char *old, char *new);
int copyf(char *old, char *new);
int appendf(char *old, char *new);
int deletef(char *old);
int close_auto(int flag);
int open_auto(int flag);
int do_auto(int iold, int isave, int itp);
int set_auto(void);
int auto_name_to_index(char *s);
int auto_par_to_name(int index, char *s);
int auto_per_par(void);
int auto_params(void);
int auto_num_par(void);
int auto_plot_par(void);
int auto_fit(void);
int auto_zoom(int i1, int j1, int i2, int j2);
int auto_xy_plot(double *x, double *y1, double *y2, double par1, double par2, double per, double *uhigh, double *ulow, double *ubar, double a);
int plot_point(int flag2, int icp1, int icp2);
int add_ps_point(double *par, double per, double *uhigh, double *ulow, double *ubar, double a, int type, int flag, int lab, int npar, int icp1, int icp2, int flag2, double *evr, double *evi);
int add_point(double *par, double per, double *uhigh, double *ulow, double *ubar, double a, int type, int flag, int lab, int npar, int icp1, int icp2, int flag2, double *evr, double *evi);
int redraw_auto_menus(void);
int get_bif_sym(char *at, int itp);
int info_header(int flag2, int icp1, int icp2);
int new_info(int ibr, int pt, char *ty, int lab, double *par, double norm, double u0, double per, int flag2, int icp1, int icp2);
int traverse_diagram(void);
int clear_auto_plot(void);
int do_auto_win(void);
int load_last_plot(int flag);
int keep_last_plot(int flag);
int init_auto_win(void);
int plot_stab(double *evr, double *evi, int n);
int clr_stab(void);
int auto_motion(XEvent ev);
int display_auto(Window w);
Window lil_button(Window root, int x, int y, char *name);
int make_auto(char *wname, char *iname);
int yes_reset_auto(void);
int reset_auto(void);
int auto_grab(void);
int auto_start_diff_ss(void);
int auto_start_at_bvp(void);
int auto_start_at_per(void);
int get_start_period(double *p);
int get_start_orbit(double *u, double t, double p, int n);
int auto_new_ss(void);
int auto_new_discrete(void);
int auto_extend_ss(void);
int auto_start_choice(void);
int torus_choice(void);
int per_doub_choice(void);
int periodic_choice(void);
int hopf_choice(void);
int auto_new_per(void);
int auto_extend_bvp(void);
int auto_switch_per(void);
int auto_switch_bvp(void);
int auto_switch_ss(void);
int auto_2p_limit(int ips);
int auto_twopar_double(void);
int auto_torus(void);
int auto_2p_branch(void);
int auto_branch_choice(int ibr);
int auto_2p_fixper(void);
int auto_2p_hopf(void);
int auto_period_double(void);
int auto_err(char *s);
int auto_run(void);
int load_auto_orbit(void);
int save_auto(void);
int save_auto_numerics(FILE *fp);
int load_auto_numerics(FILE *fp);
int save_auto_graph(FILE *fp);
int load_auto_graph(FILE *fp);
int save_q_file(FILE *fp);
int make_q_file(FILE *fp);
int noinfo(char *s);
int load_auto(void);
int move_to_label(int mylab, int *nrow, int *ndim, FILE *fp);
int get_a_row(double *u, double *t, int n, FILE *fp);
int auto_file(void);
int a_msg(int i, int v);
int auto_enter(Window w, int v);
int auto_button(XEvent ev);
int auto_kill(void);
int auto_keypress(XEvent ev, int *used);

#endif
