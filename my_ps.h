#ifndef _my_ps_h_
#define _my_ps_h_


/* my_ps.c */
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
void special_put_text_ps(int x, int y, char *str, int size);
void fancy_ps_text(int x, int y, char *str, int size, int font);
void ps_text(int x, int y, char *str);


#endif
