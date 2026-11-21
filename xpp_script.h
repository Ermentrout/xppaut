#ifndef _xpp_script_h
#define _xpp_script_h
#include <X11/Xlib.h>

/* command line */
void script_set_file(const char *name);
void script_set_options(int nopause,int exitflag);   /* -1 leaves a value unchanged */
int  script_requested(void);
void script_run_from_menu(void);   /* File -> Do script */

/* lifecycle */
int  script_load(void);
void script_start(void);        /* call just before the main event loop */
int  script_is_active(void);
void script_set_top(int flag);  /* 1 while the main loop waits for an event */

/* driven from xpp_next_event() */
int  script_poll(int *ms,XEvent *ev);    /* 0 idle, 1 ran a command, 2 timed wait, 3 *ev is a key */
int  script_filter_event(XEvent *ev);    /* 1 = event consumed */

/* prompts ask for a scripted answer before they open:
   1 = answered, 0 = open normally, -1 = script error (act as Cancel) */
int script_take_text(const char *label,char *buf,int cap);
int script_take_choice(const char *prompt,const char *key,int *index);
int script_take_form(const char *title,int n,char **names,char *vals,int stride);
int script_message(const char *msg);     /* 1 = suppress the error box */

/* commands that wait for the mouse ask for a scripted click / drag before they wait:
   1 = supplied (window pixels), 0 = wait for the mouse as usual, -1 = script error (act as Esc) */
int script_mouse_point(int *x,int *y);
int script_mouse_drag(int *x1,int *y1,int *x2,int *y2);

int    script_nopause(void);
double script_key_delay(void);
/* ---- recording: File -> Record script (see xpp_script.c) ----------------
   The hooks below are called from the places where XPP takes input.  All of
   them do nothing unless a recording is on.  A "h" is the handle returned by
   a *_begin call (0 = nothing is being recorded; pass it to the matching *_end). */
void script_record_from_menu(void);   /* File -> Record script (toggle) */
void script_record_stop(void);        /* flush and close (also called at exit) */
int  script_recording(void);
void rec_key_event(XEvent *ev);       /* a key was read (get_key_press) */
void rec_key_kind(int kind);          /* the next key read belongs to: 2 = AUTO window, 3 = data browser */
void rec_key_context(int kind);       /* 0 = none: keys read now (e.g. in a popup menu) belong to the AUTO window (2) */
void rec_nokey(int on);               /* 1/0: keys are being typed into a text field */
void rec_menu_key(int ch);            /* a menu entry was chosen with the mouse */
int  rec_text_begin(const char *init);
void rec_text_end(int h,int accepted,const char *label,const char *final);
int  rec_form_begin(int n,char **names,char *vals,int stride);
void rec_form_end(int h,int accepted,const char *title,int n,char **names,char *vals,int stride);
void rec_choice_end(int h,const char *prompt,int ch,const char *key);
void rec_click(int ix,int iy);        /* window pixels of the plot window */
void rec_drag(int x1,int y1,int x2,int y2);
void rec_esc(void);
void rec_set(const char *name,double oldv,double newv);
void rec_go(void);                    /* the Go button of the initial-conditions window */
void rec_auto_key(int ch);            /* a button of the AUTO window */
void rec_auto_traverse_begin(void);   /* Grab mode: raw keys are replaced by "auto goto" */
void rec_auto_traverse_end(void);
void rec_auto_goto(const char *sym,int label,int branch,int nth);
void rec_browse_key(const char *tok);  /* a button of the data browser */
void rec_command(const char *line);    /* one recorded command line (browser open, aplot ...) */
void rec_comment(const char *text);
int  get_key_press_nr(XEvent *ev);    /* get_key_press() that is not recorded */
#endif
