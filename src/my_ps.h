#ifndef XPPAUT_MY_PS_H
#define XPPAUT_MY_PS_H

#include <stdio.h>

/* --- Data --- */
extern char PS_FONT[100];
extern int LastPSX;
extern int LastPSY;
extern int LastPtLine;
extern int NoBreakLine;
extern int PltFmtFlag;
extern int PS_Color;
extern int PS_ColorFlag;
extern int PS_FontSize;
extern int PS_Lines;
extern double PS_LineWidth;
extern int PS_Port;


/* --- Functions --- */
void ps_bead(int x, int y);
void ps_do_color(int color);
void ps_end(void);
void ps_frect(int x, int y, int w, int h);
int ps_init(char *filename, int color);
void ps_line(int xp1, int yp1, int xp2, int yp2);
void ps_linetype(int linetype);
void ps_point(int x, int y);
void ps_restore(void);
void ps_text(int x, int y, char *str);
void special_put_text_ps(int x, int y, char *str, int size);

#endif /* XPPAUT_MY_PS_H */
