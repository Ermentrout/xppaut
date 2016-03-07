#ifndef XPPAUT_SCRNGIF_H
#define XPPAUT_SCRNGIF_H

#include <stdio.h>
#include <X11/Xlib.h>


/* --- Functions --- */
void set_global_map(int flag);
void end_ani_gif(FILE *fp);
void add_ani_gif(Window win, FILE *fp, int count);
void screen_to_gif(Window win, FILE *fp);
void get_global_colormap(Window win);

#endif /* XPPAUT_SCRNGIF_H */
