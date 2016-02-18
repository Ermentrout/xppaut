#ifndef XPPAUT_ARRAYPLOT_H
#define XPPAUT_ARRAYPLOT_H

#include <stdio.h>
#include <X11/Xlib.h>

/* --- Types --- */
typedef struct {
  Window base,wclose,wedit,wprint,wstyle,wscale,wmax,wmin,wplot,wredraw,wtime,wgif,wrange,wfit;
  int index0,indexn,alive,nacross,ndown,plotdef;
  int height,width,ploth,plotw;
  int nstart,nskip,ncskip;
  char name[20];
  double tstart,tend,zmin,zmax,dt;
  char xtitle[256],ytitle[256],filename[256],bottom[256];
  int type;
} APLOT;

/* --- Data --- */
extern int aplot_range;

/* --- Functions --- */
void set_acolor(int);
void tag_aplot(char *);
void close_aplot_files(void);
void draw_one_array_plot(char *);
void gif_aplot_all(char *,int);
void optimize_aplot(int *plist);
void make_my_aplot(char *name);
void scale_aplot(APLOT *ap, double *zmax, double *zmin);
void init_arrayplot(APLOT *ap);
void expose_aplot(Window w);
void do_array_plot_events(XEvent ev);
void wborder(Window w, int i, APLOT ap);
void destroy_aplot(void);
void init_my_aplot(void);
void create_arrayplot(APLOT *ap, char *wname, char *iname);
void print_aplot(APLOT *ap);
void apbutton(Window w);
void draw_scale(APLOT ap);
void draw_aplot(APLOT ap);
void edit_aplot(void);
void get_root(char *s, char *sroot, int *num);
void reset_aplot_axes(APLOT ap);
void dump_aplot(FILE *fp, int f);
int editaplot(APLOT *ap);
void gif_aplot(void);
void grab_aplot_screen(APLOT ap);
void redraw_aplot(APLOT ap);
void display_aplot(Window w, APLOT ap);

#endif /* XPPAUT_ARRAYPLOT_H */
