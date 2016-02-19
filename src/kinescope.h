#ifndef XPPAUT_KINESCOPE_H
#define XPPAUT_KINESCOPE_H

#include <X11/Xlib.h>

/* --- Types --- */
typedef struct {
    unsigned int h,w;
    Pixmap xi;
} MOVIE;

/* --- Functions --- */
void do_movie_com(int c);
void reset_film(void);
int film_clip(void);
int show_frame(int i, int h, int w);
void play_back(void);
void save_kine(void);
void make_anigif(void);
void save_movie(char *basename, int fmat);
void auto_play(void);
void too_small(void);

#endif /* XPPAUT_KINESCOPE_H */


