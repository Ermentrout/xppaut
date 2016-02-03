
#ifndef _my_svg_h_
#define _my_svg_h_


/* my_ps.c */
int svg_init(char *filename, int color);
void svg_stroke(void);
void svg_do_color(int color);
void svg_setcolor(int color);
void svg_end(void);
void svg_bead(int x, int y);
void svg_frect(int x, int y, int w, int h);
void svg_last_pt_off(void);
void svg_line(int xp1, int yp1, int xp2, int yp2);
void chk_svg_lines(void);
void svg_linetype(int linetype);
void svg_point(int x, int y);
void svg_write(char *str);
void svg_fnt(int cf, int scale);
void svg_show(char *str, int type);
void svg_abs(int x, int y);
void svg_rel(int x, int y);
void special_put_text_svg(int x, int y, char *str, int size);
void fancy_svg_text(int x, int y, char *str, int size, int font);
void svg_text(int x, int y, char *str);


#endif
