#ifndef XPPAUT_ANIPARSE_H
#define XPPAUT_ANIPARSE_H

#include <X11/Xlib.h>

/* --- Data --- */
extern int animation_on_the_fly;

/* --- Functions --- */
void ani_expose(Window w);
void ani_zero(void);
int check_ani_pause(XEvent ev);
void do_ani_events(XEvent ev);
void get_ani_file(char *fname);
int getppmbits(Window window, int *wid, int *hgt, unsigned char *out);
void new_vcr(void);
void on_the_fly(int task);
int writeframe(char *filename, Window window, int wid, int hgt);

#endif /* XPPAUT_ANIPARSE_H */
