#ifndef XPPAUT_ARRAYPLOT_H
#define XPPAUT_ARRAYPLOT_H

#include <stdio.h>
#include <X11/Xlib.h>
#include "xpplim.h"

/* --- Types --- */
typedef struct {
	Window base,wclose,wedit,wprint,wstyle,wscale,wmax,wmin,wplot,wredraw,wtime,wgif,wrange,wfit;
	int index0,indexn,alive,nacross,ndown,plotdef;
	int height,width,ploth,plotw;
	int nstart,nskip,ncskip;
	char name[20];
	double tstart,tend,zmin,zmax,dt;
	char xtitle[XPP_MAX_NAME],ytitle[XPP_MAX_NAME],filename[XPP_MAX_NAME],bottom[XPP_MAX_NAME];
	int type;
} APLOT;

/* --- Data --- */
extern int aplot_range;
extern APLOT aplot;

/* --- Functions --- */
void close_aplot_files(void);
void do_array_plot_events(XEvent ev);
void draw_one_array_plot(char *);
void edit_aplot(void);
void expose_aplot(Window w);
void init_my_aplot(void);
void make_my_aplot(char *name);
void optimize_aplot(int *plist);

#endif /* XPPAUT_ARRAYPLOT_H */
