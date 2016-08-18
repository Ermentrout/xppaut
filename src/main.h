#ifndef XPPAUT_MAIN_H
#define XPPAUT_MAIN_H

#include <X11/Xlib.h>
#include "load_eqn.h"

/* --- Macros --- */
/* FixWindowSize(flag) values */
#define FIX_MAX_SIZE 1
#define FIX_MIN_SIZE 2
#define FIX_SIZE 3

/* --- Data --- */
extern char anifile[XPP_MAX_NAME];
extern char batchout[XPP_MAX_NAME];
extern char big_font_name[MAXVEC];
extern char PlotFormat[MAXVEC];
extern char small_font_name[MAXVEC];
extern char UserBGBitmap[XPP_MAX_NAME];
extern char UserBlack[8];
extern char UserDrawWinColor[8];
extern char UserMainWinColor[8];
extern char UserOUTFILE[XPP_MAX_NAME];
extern char UserWhite[8];

extern float xppvermaj;
extern float xppvermin;

extern int allwinvis;
extern int batch_range;
extern int CURY_OFF;
extern int CURY_OFFb;
extern int CURY_OFFs;
extern int DCURX;
extern int DCURXs;
extern int DCURY;
extern int DCURYb;
extern int DCURYs;
extern int DisplayHeight;
extern int DisplayWidth;
extern int DoTutorial;
extern int PaperWhite;
extern int SCALEX;
extern int SCALEY;
extern int screen;
extern int tfBell;
extern int TipsFlag;
extern int TrueColorFlag;
extern int use_ani_file;
extern int UserGradients;
extern int UserMinHeight;
extern int UserMinWidth;
extern int XPPBatch;
extern int XPPVERBOSE;
extern int Xup;

extern int (*rhs)();

extern unsigned int GrFore;
extern unsigned int GrBack;
extern unsigned int MyBackColor;
extern unsigned int MyForeColor;
extern unsigned int MyDrawWinColor;
extern unsigned int MyMainWinColor;

extern Atom deleteWindowAtom;
extern Display *display;
extern GC gc;
extern GC font_gc;
extern GC gc_graph;
extern GC small_gc;
extern Window command_pop;
extern Window draw_win;
extern Window info_pop;
extern Window main_win;
extern XFontStruct *big_font;
extern XFontStruct *small_font;

/* --- Functions --- */
void bye_bye(void);
void clr_scrn(void);
void commander(int ch);
void FixWindowSize(Window w, int width, int height, int flag);
int  main(int argc, char **argv);
void redraw_all(void);
void top_button_draw(Window w);

#endif /* XPPAUT_MAIN_H */
