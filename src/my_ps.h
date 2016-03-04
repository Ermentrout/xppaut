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
extern int PS_Port;
extern int PSColorFlag;
extern int PS_FONTSIZE;
extern int PSLines;

extern double PS_LW;

/* --- Functions --- */
int ps_init(char *filename, int color);
void ps_stroke(void);
void ps_do_color(int color);
void ps_setcolor(int color);
void ps_end(void);
void ps_bead(int x, int y);
void ps_frect(int x, int y, int w, int h);
void ps_last_pt_off(void);
void ps_line(int xp1, int yp1, int xp2, int yp2);
void chk_ps_lines(void);
void ps_linetype(int linetype);
void ps_point(int x, int y);
void ps_write(char *str);
void ps_fnt(int cf, int scale);
void ps_show(char *str, int type);
void ps_abs(int x, int y);
void ps_rel(int x, int y);
void ps_write_pars(FILE *fp);
void special_put_text_ps(int x, int y, char *str, int size);
void fancy_ps_text(int x, int y, char *str, int size, int font);
void ps_text(int x, int y, char *str);

#endif /* XPPAUT_MY_PS_H */
