#ifndef XPPAUT_MY_SVG_H
#define XPPAUT_MY_SVG_H

#include <stdio.h>

/* --- Data --- */
extern FILE *svgfile;

/* --- Functions --- */
void special_put_text_svg(int x, int y, char *str, int size);
void svg_bead(int x, int y);
void svg_do_color(int color);
void svg_end(void);
void svg_frect(int x, int y, int w, int h);
int svg_init(char *filename, int color);
void svg_line(int xp1, int yp1, int xp2, int yp2);
void svg_linetype(int linetype);
void svg_point(int x, int y);
void svg_text(int x, int y, char *str);

#endif /* XPPAUT_MY_SVG_H */
