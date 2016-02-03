#ifndef _kinescope_h_
#define _kinescope_h_

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

#endif


