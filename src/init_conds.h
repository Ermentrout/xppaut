#ifndef XPPAUT_INIT_CONDS_H
#define XPPAUT_INIT_CONDS_H

#include <X11/Xlib.h>
#include "read_dir.h"

/* --- Macros --- */
/* type values */
enum xppBoxType {
	PARAMETER_BOX = 1,
	INITIAL_CONDITION_BOX,
	DELAY_BOX,
	BOUNDARY_CONDITION_BOX
};

/* --- Types --- */
typedef struct {
	int use,type,xuse;
	int n,n0;
	int nwin,minwid,minhgt;
	Window up,dn;
	Window pgup,pgdn;
	Window base;
	Window cancel,ok,def,go,close;
	Window xvt,pp,arr;
	Window *w;
	Window *we;
	Window *ck;
	char **value,*iname,*wname;
	int *isck;
	int mc,*off,*pos;
} BoxList;

/* --- Data --- */
extern char SLIDERVAR[3][20];
extern double SLIDERLO[3];
extern double SLIDERHI[3];
extern BoxList ParamBox;

/* --- Functions --- */
void box_buttons(Window w);
void box_enter_events(Window w, int yn);
void box_keypress(XEvent ev, int *used);
void clone_ode(void);
void create_par_sliders(Window base, int x0, int h0, int var_index[3]);
void draw_one_box(BoxList b, int index);
void enter_slides(Window w, int val);
void expose_box(Window w);
void expose_slides(Window w);
int file_selector(char *title, char *file, char *wild);
int find_user_name(int type, char *oname);
void initialize_box(void);
void make_new_bc_box(void);
void make_new_ic_box(void);
void make_new_delay_box(void);
void make_new_param_box(void);
void man_ic(void);
void new_parameter(void);
void redo_stuff(void);
void redraw_bcs(void);
void redraw_ics(void);
void redraw_delays(void);
void redraw_params(void);
void reset_sliders(void);
void resize_par_box(Window win);
void resize_par_slides(int h);
void set_edit_params(BoxList *b, int i, char *string);
void slide_button_press(Window w);
void slider_motion(XEvent ev);
void slide_release(Window w);

#endif /* XPPAUT_INIT_CONDS_H */
