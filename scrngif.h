#ifndef _scrngif_h_
#define _scrngif_h_

#include <stdio.h>
#include <X11/Xlib.h>

 typedef struct GifTree {
   char typ;             /* terminating, lookup, or search */
   int code;             /* the code to be output */
   unsigned char ix;     /* the color map index */
   struct GifTree **node, *nxt, *alt;
 } GifTree;
 
 
typedef struct {
  unsigned char r,g,b;
} GIFCOL;


void set_global_map(int flag);
int ppmtopix(unsigned char r, unsigned char g, unsigned char b, int *n);
void end_ani_gif(FILE *fp);
void add_ani_gif(Window win, FILE *fp, int count);
void screen_to_gif(Window win, FILE *fp);
void get_global_colormap(Window win);
void local_to_global(void);
int use_global_map(unsigned char *pixels, unsigned char *ppm, int h, int w);
int make_local_map(unsigned char *pixels, unsigned char *ppm, int h, int w);
void gif_stuff(Window win, FILE *fp, int task);
void write_global_header(int cols, int rows, FILE *dst);
void GifLoop(FILE *fout, unsigned int repeats);
void write_local_header(int cols, int rows, FILE *fout, int colflag, int delay);
void make_gif(unsigned char *pixels, int cols, int rows, FILE *dst);
int GifEncode(FILE *fout, unsigned char *pixels, int depth, int siz);
void ClearTree(int cc, GifTree *root);
unsigned char *AddCodeToBuffer(int code, short n, unsigned char *buf);

#endif
