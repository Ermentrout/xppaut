#ifndef XPPAUT_ARRAYPLOT_H
#define XPPAUT_ARRAYPLOT_H

#include <stdio.h>
#include <X11/Xlib.h>

/* --- Data --- */
extern int aplot_range;

/* --- Functions --- */
void close_aplot_files(void);
void do_array_plot_events(XEvent ev);
void draw_one_array_plot(char *);
void dump_aplot(FILE *fp, int f);
void edit_aplot(void);
void expose_aplot(Window w);
void init_my_aplot(void);
void make_my_aplot(char *name);
void optimize_aplot(int *plist);

#endif /* XPPAUT_ARRAYPLOT_H */
