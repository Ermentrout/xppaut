/*  xpp_script.c  --  scripting for XPPAUT

    A script is a plain text file, one command per line:

        # comment                      (also blank lines; a # after a blank ends the
                                       line, so "keys i g   # run it" is fine.
                                       Quote a literal #:  msg "step #1")
        set a=1.3;v=1.2;atol=1e-10     parameters, initial conditions, numerics,
                                       AUTO and plot options (anything an
                                       @-option or a parameter set accepts)
        print a                        write the value of a parameter / IC /
                                       variable to the terminal
        msg "text"                     show text in the info bar (terminal only
                                       with -nopause)
        pause ["text"]                 wait for Enter or Space (Esc aborts the
                                       script); skipped with -nopause
        wait 2.5                       wait for that many seconds; skipped with
                                       -nopause
        speed 0.2                      delay between injected keys (seconds);
                                       ignored with -nopause
        keys i g                       press keys (menu hot keys are lowercase);
                                       see below
        autokeys r                     the same, sent to the AUTO window
        browser open | close           press the Data button / the Close button of the data browser
        browsekeys s e                 the same as keys, sent to the data browser (arrows, s e r w g f l u t p a d)
        aplot edit|fit|range|redraw|print|gif|surface|close
                                       press that button of the array plot window
        ans text                       answer the next prompt (see below)
        auto goto HB 1                 AUTO diagram: select a point as Grab + Enter
        auto goto lab 5                would (see below);  auto list  prints the
                                       labelled points
        shot fig.ps [options]          save the plot (shot auto d.svg: the AUTO
                                       diagram); .ps or .svg; options color=0|1
                                       orient=land|port fontsize=N font=Name
                                       linewidth=X apply to this file only
        for a in 1 2.5 4  /  for a = 0:0.5:2  /  repeat 3  ...  end
                                       loops; $a or ${a} in the lines inside
        try ... else ... end           run the lines between try and else; if any of
                                       them fails, carry on with the else part (or
                                       after end) instead of stopping.  $error in
                                       a msg is the message of the failure.  The
                                       exit status is 3 if anything was caught.
        click x y                      answer a mouse wait (placing text, a freeze key,
        drag x1 y1 x2 y2               a start point) / a mouse drag (zoom rectangle,
                                       arrow): data coordinates of the plot; with
                                       "pixel" first, window pixels; "?" leaves it
                                       to a person
        ok                             press OK on the next prompt, keeping its
                                       default / current values
        cancel                         cancel the next prompt
        quit                           leave XPPAUT

    keys:  a token is a named key (enter tab esc space bksp del up down left
    right home end pgup pgdn), a quoted string typed literally, a single
    letter (an upper case letter is sent as lower case, since all menu hot keys
    are lower case) or any other text, which is typed literally, e.g.
    "keys 0.5 enter".  Keys are delivered one at a time, only when XPP is idle
    and waiting for input, so a key never races a redraw.

    auto goto:  "auto goto TYPE [N] [branch B]" selects the N-th point (default
    1; negative counts from the end) of type BP EP HB LP MX PD TR or UZ in the
    AUTO diagram; "auto goto lab N" the point with label N; "auto goto first" and
    "auto goto last" the ends.  "branch B" limits the search to AUTO branch B (as numbered by "auto list";
    periodic branches have negative numbers).
    The point becomes the grabbed point: its parameter and state values are loaded
    and the next AUTO command (autokeys r ...) continues from it.

    ans:  the next prompt takes its answer from the following ans line instead
    of opening.  Command-line prompts and dialogs: the text.  Multi-field forms:
    fields separated by |, positional ("0|50|-2|2"; an empty field keeps its
    value) or by label ("Xmin=0|Xmax=50").  Yes/No and two-way choices: y or n.
    "ans ?" leaves that one prompt to the person at the keyboard.  A prompt with
    no ans is left open, to be driven by following keys lines.

    Command line:   xppaut model.ode -script file.xps [-nopause] [-scriptexit]

    The script is driven from xpp_next_event(): whenever XPP is idle and
    waiting for input (no X events pending) the next command is run.  Commands
    that change program state (set, print) only run when XPP is at the top
    level, i.e. not inside a menu or dialog.  A failing command stops the script
    and reports "file:line: message"; with -scriptexit XPP then exits with a
    nonzero status.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <errno.h>
#include <strings.h>
#include <sys/time.h>
#include <sys/select.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <X11/Xutil.h>
#include "xpplim.h"
#include "struct.h"
#include "load_eqn.h"
#include "xpp_script.h"

extern Display *display;
extern int Xup;
Window browser_key_window(void);
int browser_owns_key(Window w);
Window browser_button_window(const char *name);
Window aplot_button_window(const char *name);
Window top_button_window(int i);
int file_selector(char *title,char *file,char *wild);
extern int DLeft,DRight,DTop,DBottom;
void get_draw_area(void);
extern Window main_win;
Window auto_key_window(void);
int auto_goto_point(int kind,const char *symb,int n,int label,int branch,char *err,int errsize);
void auto_list_points(void);
#define AG_TYPE 0
#define AG_LABEL 1
#define AG_FIRST 2
#define AG_LAST 3

/* from the rest of XPP */
void do_intern_set(char *name,char *value);
extern int set_option_unmatched;
int get_val(char *name,double *value);
void convert(char *old,char *new);
void get_graph(void);
void reset_graph(void);
void chk_delay(void);
void redraw_params(void);
void redraw_ics(void);
void bottom_msg(int line,char *msg);
void err_msg(char *string);
int yes_reset_auto(void);
int find_user_name(int type,char *name);
extern int AXES,FOREVER,MAXSTOR,IXPLT,IYPLT,IZPLT,MSStyle;
extern GRAPH *MyGraph;
void check_flags(void);
extern double DELTA_T,DELAY,TEND;
void redo_stuff(void);
extern double x_3d[2],y_3d[2],MY_XLO,MY_XHI,MY_YLO,MY_YHI;
extern double last_ic[];
#define XS_IC 2     /* same value as IC in newpars.h */

enum { C_SET=1, C_PRINT, C_MSG, C_PAUSE, C_WAIT, C_SPEED, C_QUIT, C_KEYS, C_AUTOKEYS, C_ANS, C_OK, C_CANCEL, C_AUTOGOTO, C_AUTOLIST, C_SHOT, C_TRY, C_ELSE, C_END, C_CLICK, C_DRAG, C_BROWSEKEYS, C_BROWSER, C_APLOT };

typedef struct {
  KeyCode kc;
  unsigned int state;
} SKEY;

typedef struct {
  int type;
  int line;
  char *arg;       /* text after the command word (quotes removed for msg/pause/ans) */
  double num;      /* wait / speed */
  SKEY *keys;      /* keys / autokeys */
  int nkeys;
  int ag_kind,ag_n,ag_label,ag_branch;   /* auto goto */
  char ag_sym[4];
  char *ctx;       /* loop iteration this command came from, e.g. "a=1.5" */
  /* shot */
  int shot_auto,shot_svg;
  unsigned int o_mask;                   /* which options were given */
  int o_color,o_port,o_fs;
  double o_lw;
  char o_font[24];
  int jump,link;   /* try / else / end: where to go, and which try a command belongs to */
  double mv[4];    /* click / drag: the coordinates */
  int mpix,mq;     /* click / drag: 1 = pixels, not data coordinates;  1 = "?" (a person does it) */
} SCMD;

static SCMD *cmds=NULL;
static const char *cur_ctx=NULL;   /* context of the line being parsed (loader) */
static const char *run_ctx=NULL;   /* context of the command being run  */

/* " [loop: a=1.5]" or "" -- appended to every message that names a line */
static void ctx_text(char *buf,int n,const char *ctx)
{
  if(ctx&&*ctx)snprintf(buf,n," [loop: %s]",ctx);
  else buf[0]=0;
}
static int ncmds=0,pc=0;
static char script_name[512]="";
static int opt_nopause=0,opt_exit=0;
static int active=0,at_top=0,paused=0,waiting=0;
static double wait_until=0.0;
static double key_delay=0.0;
static double key_gap_until=0.0;   /* next key may be sent after this time */
static int key_pos=0;              /* next key inside the current keys command */
static int wait_prompt_done=0;     /* "ans ?": wait until the person finished the prompt */

static double now_sec(void)
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec+1e-6*tv.tv_usec;
}

void script_set_file(const char *name)
{
  strncpy(script_name,name,sizeof(script_name)-1);
  script_name[sizeof(script_name)-1]=0;
}
void script_set_options(int nopause,int exitflag)
{
  if(nopause>=0)opt_nopause=nopause;
  if(exitflag>=0)opt_exit=exitflag;
  if(opt_exit)opt_nopause=1;     /* unattended: never wait for a key */
}
int script_requested(void){return script_name[0]!=0;}
int script_is_active(void){return active;}
void script_set_top(int flag){at_top=flag;}
int script_nopause(void){return opt_nopause;}
double script_key_delay(void){return opt_nopause?0.0:key_delay;}

static void script_finish(int status)
{
  active=0;
  paused=0;
  waiting=0;
  if(opt_exit){
    fflush(stdout);
    yes_reset_auto();
    exit(status);
  }
}

/* ------------------------------------------------------------- try / else
   A failure inside a "try" block does not stop the script: it is logged,
   remembered (for $error), and the script goes on with the else part (or after
   the block).  XPP may be left with a menu or dialog open by the failed
   command, so afterwards Esc is pressed until it is back at the top level
   (see "recovering" in script_poll).                                         */
#define MAX_TRY 16
static int try_stack[MAX_TRY],try_sp=0;
static int caught_total=0;               /* failures caught so far */
static char last_error[400]="";
static int recovering=0,recover_count=0;

static int catch_failure(const char *full,const char *plain)
{
  int t;
  if(try_sp==0||recovering)return 0;
  t=try_stack[--try_sp];
  caught_total++;
  strncpy(last_error,plain,sizeof(last_error)-1);
  last_error[sizeof(last_error)-1]=0;
  fprintf(stderr,"script caught  %s\n",full);
  fflush(stderr);
  pc=cmds[t].jump;
  paused=0;waiting=0;wait_prompt_done=0;key_gap_until=0.0;key_pos=0;
  recovering=1;recover_count=0;
  return 1;
}

/* replace $error / ${error} in the text of msg and pause by the last caught message */
static void expand_error(const char *in,char *out,int n)
{
  int o=0;
  const char *v;
  while(*in&&o<n-1){
    if(strncmp(in,"${error}",8)==0){
      for(v=last_error;*v&&o<n-1;v++)out[o++]=*v;
      in+=8;
    }
    else if(strncmp(in,"$error",6)==0&&!isalnum((unsigned char)in[6])&&in[6]!='_'){
      for(v=last_error;*v&&o<n-1;v++)out[o++]=*v;
      in+=6;
    }
    else out[o++]=*in++;
  }
  out[o]=0;
}

/* stop the script with a message; no dialog is shown */
static void script_fail_quiet(int line,const char *fmt,...)
{
  char msg[400],full[600];
  va_list ap;
  va_start(ap,fmt);
  vsnprintf(msg,sizeof(msg),fmt,ap);
  va_end(ap);
  {char cx[300];ctx_text(cx,sizeof(cx),run_ctx);
   snprintf(full,sizeof(full),"%s:%d: %s%s",script_name,line,msg,cx);}
  if(catch_failure(full,msg))return;
  fprintf(stderr,"script error  %s\n",full);
  fflush(stderr);
  active=0;
  paused=0;
  waiting=0;
  if(opt_exit){
    yes_reset_auto();
    exit(1);
  }
}

/* report an error and stop the script; interactively also tell the person */
static void script_fail(int line,const char *fmt,...)
{
  char msg[400],full[600];
  va_list ap;
  va_start(ap,fmt);
  vsnprintf(msg,sizeof(msg),fmt,ap);
  va_end(ap);
  {char cx[300];ctx_text(cx,sizeof(cx),run_ctx);
   snprintf(full,sizeof(full),"%s:%d: %s%s",script_name,line,msg,cx);}
  if(catch_failure(full,msg))return;
  fprintf(stderr,"script error  %s\n",full);
  fflush(stderr);
  active=0;
  paused=0;
  waiting=0;
  if(opt_exit){
    yes_reset_auto();
    exit(1);
  }
  if(Xup)err_msg(full);
}

/* ---------------------------------------------------------------- parsing */

static char *skip_ws(char *s){while(*s&&isspace((unsigned char)*s))s++;return s;}

static void rtrim(char *s)
{
  int n=strlen(s);
  while(n>0&&isspace((unsigned char)s[n-1]))s[--n]=0;
}

/* Remove a trailing comment: a '#' that starts the line or follows a blank,
   outside quotes.  (Use quotes for a literal #:  msg "step #1".) */
static void strip_comment(char *s)
{
  int i,q=0;
  for(i=0;s[i];i++){
    if(q){ if(s[i]==q)q=0; continue; }
    if(s[i]=='"'||s[i]=='\''){ q=s[i]; continue; }
    if(s[i]=='#'&&(i==0||isspace((unsigned char)s[i-1]))){ s[i]=0; return; }
  }
}

/* remove one pair of surrounding quotes, if present */
static void unquote(char *s)
{
  int n=strlen(s);
  if(n>=2&&((s[0]=='"'&&s[n-1]=='"')||(s[0]=='\''&&s[n-1]=='\''))){
    memmove(s,s+1,n-2);
    s[n-2]=0;
  }
}

static int cmds_cap=0;   /* allocated size of cmds[]: grown geometrically */

static int add_cmd(int type,int line,const char *arg,double num)
{
  if(ncmds>=cmds_cap){
    int ncap=cmds_cap?2*cmds_cap:64;
    SCMD *c=(SCMD *)realloc(cmds,sizeof(SCMD)*(size_t)ncap);
    if(c==NULL)return 0;
    cmds=c;
    cmds_cap=ncap;
  }
  cmds[ncmds].type=type;
  cmds[ncmds].line=line;
  cmds[ncmds].arg=strdup(arg?arg:"");
  cmds[ncmds].num=num;
  cmds[ncmds].keys=NULL;
  cmds[ncmds].nkeys=0;
  cmds[ncmds].ag_kind=0;
  cmds[ncmds].ag_n=1;
  cmds[ncmds].ag_label=0;
  cmds[ncmds].ag_branch=0;
  cmds[ncmds].ag_sym[0]=0;
  cmds[ncmds].ctx=(cur_ctx&&*cur_ctx)?strdup(cur_ctx):NULL;
  cmds[ncmds].shot_auto=0;
  cmds[ncmds].shot_svg=0;
  cmds[ncmds].o_mask=0;
  cmds[ncmds].o_color=0;
  cmds[ncmds].o_port=0;
  cmds[ncmds].o_fs=0;
  cmds[ncmds].o_lw=0.0;
  cmds[ncmds].o_font[0]=0;
  cmds[ncmds].jump=0;
  cmds[ncmds].link=0;
  cmds[ncmds].mv[0]=cmds[ncmds].mv[1]=cmds[ncmds].mv[2]=cmds[ncmds].mv[3]=0.0;
  cmds[ncmds].mpix=0;
  cmds[ncmds].mq=0;
  ncmds++;
  return 1;
}


/* ------------------------------------------------------------------ keys */

static struct { const char *name; KeySym ks; } named_keys[]={
  {"enter",XK_Return},{"return",XK_Return},{"tab",XK_Tab},
  {"esc",XK_Escape},{"escape",XK_Escape},{"space",XK_space},
  {"bksp",XK_BackSpace},{"backspace",XK_BackSpace},
  {"del",XK_Delete},{"delete",XK_Delete},
  {"up",XK_Up},{"down",XK_Down},{"left",XK_Left},{"right",XK_Right},
  {"home",XK_Home},{"end",XK_End},{"pgup",XK_Prior},{"pgdn",XK_Next},
  {NULL,0}
};

/* find a keycode and shift state that XLookupString() turns back into c */
static int resolve_char(char c,SKEY *k)
{
  KeySym ks=(c==' ')?XK_space:(KeySym)(unsigned char)c;
  KeyCode kc=XKeysymToKeycode(display,ks);
  XKeyEvent e;
  char buf[16];
  KeySym sym;
  unsigned int st[2]={0,ShiftMask};
  int i,n;
  if(kc==0)return 0;
  for(i=0;i<2;i++){
    memset(&e,0,sizeof(e));
    e.type=KeyPress;
    e.display=display;
    e.keycode=kc;
    e.state=st[i];
    n=XLookupString(&e,buf,sizeof(buf),&sym,NULL);
    if(n==1&&buf[0]==c){k->kc=kc;k->state=st[i];return 1;}
  }
  return 0;
}

static int resolve_keysym(KeySym ks,SKEY *k)
{
  KeyCode kc=XKeysymToKeycode(display,ks);
  if(kc==0)return 0;
  k->kc=kc;
  k->state=0;
  return 1;
}

static int push_key(SKEY **keys,int *n,SKEY k)
{
  SKEY *p=(SKEY *)realloc(*keys,sizeof(SKEY)*(*n+1));
  if(p==NULL)return 0;
  p[*n]=k;
  *keys=p;
  (*n)++;
  return 1;
}

/* Parse the argument of keys/autokeys.  Returns 0 and fills *keys,*n, or
   nonzero with a message in err. */
static int parse_keys(const char *arg,SKEY **keys,int *n,char *err)
{
  const char *p=arg;
  char tok[256],low[256];
  int i,len,j,found;
  SKEY k;
  *keys=NULL;*n=0;
  while(*p){
    while(*p&&isspace((unsigned char)*p))p++;
    if(*p==0)break;
    if(*p=='"'||*p=='\''){
      char q=*p++;
      len=0;
      while(*p&&*p!=q&&len<255)tok[len++]=*p++;
      if(*p!=q){sprintf(err,"unterminated quote in keys");return 1;}
      p++;
      tok[len]=0;
      /* literal text, typed exactly */
      for(i=0;i<len;i++){
        if(tok[i]<' '||tok[i]>'~'){sprintf(err,"cannot type a non-printable character");return 1;}
        if(!resolve_char(tok[i],&k)){sprintf(err,"cannot type '%c' on this keyboard",tok[i]);return 1;}
        push_key(keys,n,k);
      }
      continue;
    }
    len=0;
    while(*p&&!isspace((unsigned char)*p)&&len<255)tok[len++]=*p++;
    tok[len]=0;
    for(i=0;i<=len;i++)low[i]=(char)tolower((unsigned char)tok[i]);
    found=0;
    for(j=0;named_keys[j].name!=NULL;j++){
      if(strcmp(low,named_keys[j].name)==0){
        if(!resolve_keysym(named_keys[j].ks,&k)){sprintf(err,"no key for \"%s\" on this keyboard",tok);return 1;}
        push_key(keys,n,k);
        found=1;
        break;
      }
    }
    if(found)continue;
    if(len==1){
      /* a single letter: menu hot keys are lower case */
      char c=(char)tolower((unsigned char)tok[0]);
      if(!resolve_char(c,&k)){sprintf(err,"cannot type '%c' on this keyboard",c);return 1;}
      push_key(keys,n,k);
      continue;
    }
    /* any other token is typed exactly as written */
    for(i=0;i<len;i++){
      if(tok[i]<' '||tok[i]>'~'){sprintf(err,"cannot type a non-printable character");return 1;}
      if(!resolve_char(tok[i],&k)){sprintf(err,"cannot type '%c' on this keyboard",tok[i]);return 1;}
      push_key(keys,n,k);
    }
  }
  if(*n==0){sprintf(err,"keys needs at least one key");return 1;}
  return 0;
}

/* Errors found while loading a script: printed, and the first one is kept so
   that File -> Do script can show it in a box. */
static char load_err_text[400]="";
static void load_error(const char *fmt,...)
{
  char msg[400];
  va_list ap;
  va_start(ap,fmt);
  vsnprintf(msg,sizeof(msg),fmt,ap);
  va_end(ap);
  {char cx[300];ctx_text(cx,sizeof(cx),cur_ctx);strncat(msg,cx,sizeof(msg)-strlen(msg)-1);}
  fprintf(stderr,"script error  %s\n",msg);
  fflush(stderr);
  if(load_err_text[0]==0)strcpy(load_err_text,msg);
}

static void free_cmds(void)
{
  int i;
  for(i=0;i<ncmds;i++){free(cmds[i].arg);free(cmds[i].keys);free(cmds[i].ctx);}
  free(cmds);
  cmds=NULL;
  ncmds=0;
  cmds_cap=0;
}


/* =========================================================================
   Loops.  Before any command is parsed, "for" and "repeat" blocks are expanded
   into plain lines with the loop variables substituted ("$a" or "${a}"; "$$"
   is a literal $).  Every iteration is therefore syntax-checked when the script
   is loaded, and each command remembers which iteration it belongs to so that
   messages can say so.

        for a in 1.0 1.5 2.0          a list (words, numbers, quoted text)
        for a = 0:0.5:2               start:step:stop, or start:stop (step 1)
        repeat 3                      no variable
          ...
        end
   ========================================================================= */
typedef struct { char *text; int line; } RLINE;
typedef struct { char *text; int line; char *ctx; } ELINE;
typedef struct { char name[40]; char value[160]; } LVAR;

static RLINE *rl=NULL;  static int nrl=0;      /* the lines of the file          */
static ELINE *el=NULL;  static int nel=0;      /* the lines after expansion      */
static int *rmatch=NULL;                       /* for/repeat line -> its "end"   */
static int rl_cap=0,el_cap=0;                 /* allocated sizes (grown geometrically) */
#define MAX_EXPANDED 200000
#define MAX_LOOP_VALUES 100000
#define MAX_LOOP_DEPTH 8

static void free_expanded(void)
{
  int i;
  for(i=0;i<nel;i++){free(el[i].text);free(el[i].ctx);}
  free(el);el=NULL;nel=0;el_cap=0;
  for(i=0;i<nrl;i++)free(rl[i].text);
  free(rl);rl=NULL;nrl=0;rl_cap=0;
  free(rmatch);rmatch=NULL;
}

static int rl_add(const char *text,int line)
{
  if(nrl>=rl_cap){
    int ncap=rl_cap?2*rl_cap:64;
    RLINE *p=(RLINE *)realloc(rl,sizeof(RLINE)*(size_t)ncap);
    if(p==NULL)return 0;
    rl=p;
    rl_cap=ncap;
  }
  rl[nrl].text=strdup(text);
  rl[nrl].line=line;
  nrl++;
  return 1;
}

static int el_add(const char *text,int line,const char *ctx)
{
  ELINE *p;
  if(nel>=MAX_EXPANDED){
    load_error("%s:%d: the loops expand to more than %d lines",script_name,line,MAX_EXPANDED);
    return 0;
  }
  if(nel>=el_cap){
    int ncap=el_cap?2*el_cap:256;
    p=(ELINE *)realloc(el,sizeof(ELINE)*(size_t)ncap);
    if(p==NULL)return 0;
    el=p;
    el_cap=ncap;
  }
  el[nel].text=strdup(text);
  el[nel].line=line;
  el[nel].ctx=(ctx&&*ctx)?strdup(ctx):NULL;
  nel++;
  return 1;
}

static void first_word(const char *s,char *w,int n)
{
  int i=0;
  while(*s&&isspace((unsigned char)*s))s++;
  while(s[i]&&!isspace((unsigned char)s[i])&&i<n-1){w[i]=(char)tolower((unsigned char)s[i]);i++;}
  w[i]=0;
}

/* next blank-separated token; a quoted token may contain blanks (quotes removed) */
static int tok_next(const char **pp,char *tok,int size)
{
  const char *p=*pp;
  int n=0;
  char q=0;
  while(*p&&isspace((unsigned char)*p))p++;
  if(*p==0){*pp=p;return 0;}
  if(*p=='"'||*p=='\''){
    q=*p++;
    while(*p&&*p!=q&&n<size-1)tok[n++]=*p++;
    if(*p==q)p++;
  }
  else while(*p&&!isspace((unsigned char)*p)&&n<size-1)tok[n++]=*p++;
  tok[n]=0;
  *pp=p;
  return 1;
}

/* replace $name, ${name} by the loop variable's value; "$$" is a literal $ */
static int substitute(const char *in,char *out,int outsize,const LVAR *vars,int nv,int line)
{
  int i=0,o=0,k;
  while(in[i]){
    if(o>=outsize-2){
      load_error("%s:%d: line too long after substituting the loop variables",script_name,line);
      return 0;
    }
    if(in[i]=='$'){
      if(in[i+1]=='$'){out[o++]='$';i+=2;continue;}
      {
        int braces=(in[i+1]=='{');
        const char *q=in+i+1+braces;
        if(isalpha((unsigned char)*q)||*q=='_'){
          char name[40];
          int nl=0;
          const char *v;
          while((isalnum((unsigned char)q[nl])||q[nl]=='_')&&nl<39){name[nl]=q[nl];nl++;}
          name[nl]=0;
          if(braces&&q[nl]!='}'){
            load_error("%s:%d: unterminated ${%s",script_name,line,name);
            return 0;
          }
          for(k=nv-1;k>=0;k--)if(strcmp(vars[k].name,name)==0)break;
          if(k<0){
            if(strcmp(name,"error")==0){
              /* $error is the message of the failure caught by try; it is filled in
                 when the else part runs, so it is kept as written */
              const char *tt=braces?"${error}":"$error";
              while(*tt&&o<outsize-2)out[o++]=*tt++;
              i+=1+braces+nl+braces;
              continue;
            }
            load_error("%s:%d: undefined loop variable $%s (write $$ for a literal $)",script_name,line,name);
            return 0;
          }
          v=vars[k].value;
          while(*v&&o<outsize-2)out[o++]=*v++;
          i+=1+braces+nl+braces;
          continue;
        }
      }
    }
    out[o++]=in[i++];
  }
  out[o]=0;
  return 1;
}

static int is_loop_word(const char *w){return strcmp(w,"for")==0||strcmp(w,"repeat")==0;}

/* pair every for/repeat with its end, and check try / else / end; 1 if the blocks are well formed */
static int match_blocks(void)
{
  int i,sp=0,*stack,*kind,*had_else;
  char w[16];
  rmatch=(int *)calloc(nrl>0?nrl:1,sizeof(int));
  stack=(int *)calloc(nrl>0?nrl:1,sizeof(int));
  kind=(int *)calloc(nrl>0?nrl:1,sizeof(int));        /* 0 = loop, 1 = try */
  had_else=(int *)calloc(nrl>0?nrl:1,sizeof(int));
  if(rmatch==NULL||stack==NULL||kind==NULL||had_else==NULL){free(stack);free(kind);free(had_else);return 0;}
  for(i=0;i<nrl;i++){
    first_word(rl[i].text,w,16);
    if(is_loop_word(w)){stack[sp]=i;kind[sp]=0;had_else[sp]=0;sp++;}
    else if(strcmp(w,"try")==0){stack[sp]=i;kind[sp]=1;had_else[sp]=0;sp++;}
    else if(strcmp(w,"else")==0){
      if(sp==0||kind[sp-1]!=1){
        load_error("%s:%d: \"else\" without a matching try",script_name,rl[i].line);
        free(stack);free(kind);free(had_else);return 0;
      }
      if(had_else[sp-1]){
        load_error("%s:%d: this try already has an else",script_name,rl[i].line);
        free(stack);free(kind);free(had_else);return 0;
      }
      had_else[sp-1]=1;
    }
    else if(strcmp(w,"end")==0){
      if(sp==0){
        load_error("%s:%d: \"end\" without a matching for, repeat or try",script_name,rl[i].line);
        free(stack);free(kind);free(had_else);return 0;
      }
      sp--;
      if(kind[sp]==0)rmatch[stack[sp]]=i;
    }
  }
  if(sp>0){
    load_error("%s:%d: this %s has no matching \"end\"",script_name,rl[stack[sp-1]].line,kind[sp-1]?"try":"loop");
    free(stack);free(kind);free(had_else);return 0;
  }
  free(stack);free(kind);free(had_else);
  return 1;
}

static void format_value(char *out,int n,double v,double scale)
{
  if(fabs(v)<1e-12*(scale>1.0?scale:1.0))v=0.0;
  snprintf(out,n,"%.12g",v);
}

/* The values of one loop.  Returns how many (>0) with *vals malloc'ed, or 0
   after reporting an error.  name receives the variable ("" for repeat). */
static int loop_values(const char *hdr,int line,char *name,char ***valsout)
{
  const char *p=hdr;
  char t[200],w[16];
  char **vals=NULL;
  int n=0;
  name[0]=0;
  first_word(p,w,16);
  tok_next(&p,t,sizeof(t));              /* the word for / repeat */
  if(strcmp(w,"repeat")==0){
    char *e;
    long cnt;
    if(!tok_next(&p,t,sizeof(t))){load_error("%s:%d: repeat needs a count",script_name,line);return 0;}
    cnt=strtol(t,&e,10);
    if(*e||cnt<1||cnt>MAX_LOOP_VALUES){
      load_error("%s:%d: \"%s\" is not a valid repeat count (1 to %d)",script_name,line,t,MAX_LOOP_VALUES);
      return 0;
    }
    if(tok_next(&p,t,sizeof(t))){load_error("%s:%d: unexpected \"%s\" after repeat count",script_name,line,t);return 0;}
    vals=(char **)calloc((size_t)cnt,sizeof(char *));
    if(vals==NULL)return 0;
    for(n=0;n<(int)cnt;n++)vals[n]=strdup("");
    *valsout=vals;
    return (int)cnt;
  }
  /* for NAME in LIST   |   for NAME = a:b[:c] */
  {
    int nl=0;
    while(isspace((unsigned char)*p))p++;
    while((isalnum((unsigned char)*p)||*p=='_')&&nl<39){name[nl++]=*p++;}
    name[nl]=0;
    if(nl==0||isdigit((unsigned char)name[0])){
      load_error("%s:%d: for needs a variable name: for a in 1 2 3   or   for a = 0:0.5:2",script_name,line);
      return 0;
    }
    while(isspace((unsigned char)*p))p++;
    if(*p=='='){
      char spec[300],*part[3],*c1,*c2,*e;
      double a,b,step,cnt;
      int np=1,k;
      p++;
      while(isspace((unsigned char)*p))p++;
      strncpy(spec,p,sizeof(spec)-1);spec[sizeof(spec)-1]=0;
      rtrim(spec);
      part[0]=spec;
      for(c1=spec;*c1;c1++)if(*c1==':'&&np<3){*c1=0;part[np++]=c1+1;}
      if(np<2){load_error("%s:%d: a range is start:stop or start:step:stop, e.g. for a = 0:0.5:2",script_name,line);return 0;}
      a=strtod(part[0],&e);if(*e||*part[0]==0){load_error("%s:%d: \"%s\" is not a number",script_name,line,part[0]);return 0;}
      if(np==2){step=1.0;b=strtod(part[1],&e);if(*e||*part[1]==0){load_error("%s:%d: \"%s\" is not a number",script_name,line,part[1]);return 0;}}
      else {
        step=strtod(part[1],&e);if(*e||*part[1]==0){load_error("%s:%d: \"%s\" is not a number",script_name,line,part[1]);return 0;}
        b=strtod(part[2],&e);if(*e||*part[2]==0){load_error("%s:%d: \"%s\" is not a number",script_name,line,part[2]);return 0;}
      }
      (void)c2;
      if(step==0.0){load_error("%s:%d: the step of a range cannot be 0",script_name,line);return 0;}
      cnt=(b-a)/step;
      if(cnt<-1e-9){load_error("%s:%d: empty range (start %g, step %g, stop %g): the step has the wrong sign",script_name,line,a,step,b);return 0;}
      cnt=floor(cnt+1e-9)+1.0;
      if(cnt>MAX_LOOP_VALUES){load_error("%s:%d: the range has more than %d values",script_name,line,MAX_LOOP_VALUES);return 0;}
      n=(int)cnt;
      vals=(char **)calloc((size_t)n,sizeof(char *));
      if(vals==NULL)return 0;
      for(k=0;k<n;k++){
        char buf[64];
        format_value(buf,sizeof(buf),a+k*step,fabs(b)>fabs(a)?fabs(b):fabs(a));
        vals[k]=strdup(buf);
      }
      *valsout=vals;
      return n;
    }
    if(!tok_next(&p,t,sizeof(t))||strcmp(t,"in")!=0){
      load_error("%s:%d: expected \"in\" or \"=\" after the variable name",script_name,line);
      return 0;
    }
    vals=NULL;
    while(tok_next(&p,t,sizeof(t))){
      char **nv=(char **)realloc(vals,sizeof(char *)*(n+1));
      if(nv==NULL||n>=MAX_LOOP_VALUES){
        load_error("%s:%d: too many values in the list",script_name,line);
        return 0;
      }
      vals=nv;
      vals[n++]=strdup(t);
    }
    if(n==0){load_error("%s:%d: the list after \"in\" is empty",script_name,line);return 0;}
    *valsout=vals;
    return n;
  }
}

/* expand rl[i0..i1) into el[] */
static int expand(int i0,int i1,const LVAR *vars,int nv,const char *ctx)
{
  int i;
  char word[16];
  static char buf[8192];
  for(i=i0;i<i1;i++){
    first_word(rl[i].text,word,16);
    if(is_loop_word(word)){
      char hdr[2048],name[40],c2[512],piece[260];
      char **vals=NULL;
      int j=rmatch[i],k,nvals,nn,m;
      LVAR nvars[MAX_LOOP_DEPTH];
      if(nv>=MAX_LOOP_DEPTH){
        load_error("%s:%d: loops are nested more than %d deep",script_name,rl[i].line,MAX_LOOP_DEPTH);
        return 0;
      }
      if(!substitute(rl[i].text,hdr,sizeof(hdr),vars,nv,rl[i].line))return 0;
      nvals=loop_values(hdr,rl[i].line,name,&vals);
      if(nvals<=0)return 0;
      if(nv>0)memcpy(nvars,vars,sizeof(LVAR)*nv);
      for(k=0;k<nvals;k++){
        nn=nv;
        if(name[0]){
          strncpy(nvars[nn].name,name,39);nvars[nn].name[39]=0;
          strncpy(nvars[nn].value,vals[k],159);nvars[nn].value[159]=0;
          nn++;
          snprintf(piece,sizeof(piece),"%s=%s",name,vals[k]);
        }
        else snprintf(piece,sizeof(piece),"repeat %d of %d",k+1,nvals);
        if(ctx&&*ctx)snprintf(c2,sizeof(c2),"%s, %s",ctx,piece);
        else snprintf(c2,sizeof(c2),"%s",piece);
        if(!expand(i+1,j,nvars,nn,c2)){
          for(m=0;m<nvals;m++)free(vals[m]);
          free(vals);
          return 0;
        }
      }
      for(m=0;m<nvals;m++)free(vals[m]);
      free(vals);
      i=j;
      continue;
    }
    if(!substitute(rl[i].text,buf,sizeof(buf),vars,nv,rl[i].line))return 0;
    if(!el_add(buf,rl[i].line,ctx))return 0;
  }
  return 1;
}

/* Read and check the script.  Returns 0 on success, otherwise prints the
   problem and returns nonzero.  */
int script_load(void)
{
  FILE *fp;
  char buf[2048],*p,*arg,word[64];
  int line=0,i,n,bad=0,kk;
  int pstack[MAX_TRY],pelse[MAX_TRY],psp=0;
  char *end;
  double v;

  if(script_name[0]==0)return 0;
  free_cmds();
  load_err_text[0]=0;
  free_expanded();
  fp=fopen(script_name,"r");
  if(fp==NULL){
    load_error("cannot open %s: %s",script_name,strerror(errno));
    return 1;
  }
  /* 1. read the lines (comments and blank lines dropped) */
  while(fgets(buf,sizeof(buf),fp)!=NULL){
    line++;
    p=skip_ws(buf);
    strip_comment(p);
    rtrim(p);
    if(*p==0)continue;
    rl_add(p,line);
  }
  fclose(fp);
  /* 2. expand the loops */
  cur_ctx=NULL;
  if(!match_blocks()||!expand(0,nrl,NULL,0,"")){
    free_expanded();
    return 1;
  }
  /* 3. parse every (expanded) line */
  for(kk=0;kk<nel;kk++){
    line=el[kk].line;
    cur_ctx=el[kk].ctx;
    strncpy(buf,el[kk].text,sizeof(buf)-1);
    buf[sizeof(buf)-1]=0;
    p=skip_ws(buf);
    strip_comment(p);
    rtrim(p);
    if(*p==0)continue;
    n=0;
    while(p[n]&&!isspace((unsigned char)p[n])&&n<63){word[n]=(char)tolower((unsigned char)p[n]);n++;}
    word[n]=0;
    arg=skip_ws(p+n);
    if(strcmp(word,"set")==0){
      if(*arg==0){load_error("%s:%d: set needs name=value",script_name,line);bad=1;continue;}
      unquote(arg);
      add_cmd(C_SET,line,arg,0.0);
    }
    else if(strcmp(word,"print")==0){
      if(*arg==0){load_error("%s:%d: print needs a name",script_name,line);bad=1;continue;}
      add_cmd(C_PRINT,line,arg,0.0);
    }
    else if(strcmp(word,"msg")==0){
      unquote(arg);
      add_cmd(C_MSG,line,arg,0.0);
    }
    else if(strcmp(word,"pause")==0){
      unquote(arg);
      add_cmd(C_PAUSE,line,arg,0.0);
    }
    else if(strcmp(word,"wait")==0||strcmp(word,"speed")==0){
      v=strtod(arg,&end);
      if(end==arg||v<0.0){
        load_error("%s:%d: %s needs a number of seconds",script_name,line,word);
        bad=1;continue;
      }
      add_cmd(word[0]=='w'?C_WAIT:C_SPEED,line,arg,v);
    }
    else if(strcmp(word,"click")==0||strcmp(word,"drag")==0){
      /* click [pixel] x y     drag [pixel] x1 y1 x2 y2     click ?   drag ?  */
      int isdrag=(word[0]=='d'),need=isdrag?4:2,got=0,pix=0,ask=0,ok=1,first=1;
      double mv[4]={0.0,0.0,0.0,0.0};
      const char *q=arg;
      char tok[80];
      while(tok_next(&q,tok,sizeof(tok))){
        if(strcmp(tok,"?")==0&&first){ask=1;}
        else if((strcasecmp(tok,"pixel")==0||strcasecmp(tok,"pixels")==0||strcasecmp(tok,"px")==0)&&first){pix=1;}
        else {
          v=strtod(tok,&end);
          if(end==tok||*end!=0||got>=need||ask){ok=0;break;}
          mv[got++]=v;
        }
        first=0;
      }
      if(!ok||(!ask&&got!=need)){
        load_error(isdrag?"%s:%d: drag needs x1 y1 x2 y2 (data coordinates; \"drag pixel ...\" for window pixels; \"drag ?\" leaves it to a person)"
                        :"%s:%d: click needs x y (data coordinates; \"click pixel x y\" for window pixels; \"click ?\" leaves it to a person)",
                   script_name,line);
        bad=1;continue;
      }
      add_cmd(isdrag?C_DRAG:C_CLICK,line,arg,0.0);
      {SCMD *c=&cmds[ncmds-1];c->mv[0]=mv[0];c->mv[1]=mv[1];c->mv[2]=mv[2];c->mv[3]=mv[3];c->mpix=pix;c->mq=ask;}
    }
    else if(strcmp(word,"quit")==0){
      add_cmd(C_QUIT,line,"",0.0);
    }
    else if(strcmp(word,"keys")==0||strcmp(word,"autokeys")==0||strcmp(word,"browsekeys")==0){
      SKEY *kk;int nk;char err[200];
      if(!Xup){load_error("%s:%d: %s needs the GUI",script_name,line,word);bad=1;continue;}
      if(parse_keys(arg,&kk,&nk,err)){
        load_error("%s:%d: %s",script_name,line,err);
        bad=1;continue;
      }
      if(!add_cmd(word[0]=='k'?C_KEYS:(word[0]=='a'?C_AUTOKEYS:C_BROWSEKEYS),line,arg,0.0)){bad=1;continue;}
      cmds[ncmds-1].keys=kk;
      cmds[ncmds-1].nkeys=nk;
    }
    else if(strcmp(word,"browser")==0||strcmp(word,"aplot")==0){
      static const char *bnames[]={"open","close",NULL};
      static const char *anames[]={"edit","fit","range","redraw","print","gif","surface","close",NULL};
      const char **names=(word[0]=='b')?bnames:anames;
      char sub[64];
      int k,found=0;
      if(!Xup){load_error("%s:%d: %s needs the GUI",script_name,line,word);bad=1;continue;}
      strncpy(sub,arg,sizeof(sub)-1);sub[sizeof(sub)-1]=0;
      rtrim(sub);
      for(k=0;sub[k];k++)sub[k]=(char)tolower((unsigned char)sub[k]);
      for(k=0;names[k]!=NULL;k++)if(strcmp(sub,names[k])==0)found=1;
      if(!found){
        load_error(word[0]=='b'?"%s:%d: browser needs open or close":
                   "%s:%d: aplot needs one of edit fit range redraw print gif surface close",script_name,line);
        bad=1;continue;
      }
      add_cmd(word[0]=='b'?C_BROWSER:C_APLOT,line,sub,0.0);
    }
    else if(strcmp(word,"auto")==0){
      char tmp[512],*tk[8],*t;
      int nt=0,idx,k,kind=AG_TYPE,nn=1,lab=0,br=0,okparse=1;
      char sym[4]="";
      static const char *types[]={"BP","EP","HB","LP","MX","PD","TR","UZ",NULL};
      strncpy(tmp,arg,sizeof(tmp)-1);tmp[sizeof(tmp)-1]=0;
      for(t=strtok(tmp," \t");t!=NULL&&nt<8;t=strtok(NULL," \t"))tk[nt++]=t;
      for(k=0;k<nt;k++){char *q;for(q=tk[k];*q;q++)*q=(char)tolower((unsigned char)*q);}
      if(nt==0){load_error("%s:%d: auto needs a subcommand: goto or list",script_name,line);bad=1;continue;}
      if(strcmp(tk[0],"list")==0){
        if(nt>1){load_error("%s:%d: auto list takes no arguments",script_name,line);bad=1;continue;}
        add_cmd(C_AUTOLIST,line,arg,0.0);
        continue;
      }
      if(strcmp(tk[0],"goto")!=0){
        load_error("%s:%d: unknown auto subcommand \"%s\" (use goto or list)",script_name,line,tk[0]);
        bad=1;continue;
      }
      if(nt<2){load_error("%s:%d: auto goto needs a point, e.g. HB 1, lab 3, first or last",script_name,line);bad=1;continue;}
      idx=2;
      if(strcmp(tk[1],"first")==0)kind=AG_FIRST;
      else if(strcmp(tk[1],"last")==0)kind=AG_LAST;
      else if(strcmp(tk[1],"lab")==0||strcmp(tk[1],"label")==0){
        char *e;
        kind=AG_LABEL;
        if(nt<3){load_error("%s:%d: auto goto lab needs a label number",script_name,line);bad=1;continue;}
        lab=(int)strtol(tk[2],&e,10);
        if(*e||lab<1){load_error("%s:%d: \"%s\" is not a valid label number",script_name,line,tk[2]);bad=1;continue;}
        idx=3;
      }
      else {
        for(k=0;types[k]!=NULL;k++){
          char up[4];int m;
          for(m=0;m<2;m++)up[m]=(char)toupper((unsigned char)tk[1][m]);
          up[2]=0;
          if(strlen(tk[1])==2&&strcmp(up,types[k])==0){strcpy(sym,types[k]);break;}
        }
        if(sym[0]==0){
          load_error("%s:%d: unknown point type \"%s\" (use BP EP HB LP MX PD TR UZ, lab N, first or last)",script_name,line,tk[1]);
          bad=1;continue;
        }
        kind=AG_TYPE;
        if(idx<nt&&strcmp(tk[idx],"branch")!=0&&strcmp(tk[idx],"br")!=0){
          char *e;
          nn=(int)strtol(tk[idx],&e,10);
          if(*e||nn==0){load_error("%s:%d: \"%s\" is not a valid point number (1, 2, ... or -1 for the last)",script_name,line,tk[idx]);okparse=0;}
          idx++;
        }
      }
      if(!okparse){bad=1;continue;}
      if(idx<nt){
        char *e;
        if((strcmp(tk[idx],"branch")==0||strcmp(tk[idx],"br")==0)&&idx+1<nt){
          br=(int)strtol(tk[idx+1],&e,10);
          if(*e||br==0){load_error("%s:%d: \"%s\" is not a valid branch number (see auto list; periodic branches are negative)",script_name,line,tk[idx+1]);bad=1;continue;}
          idx+=2;
        }
      }
      if(idx<nt){load_error("%s:%d: unexpected \"%s\" in auto goto",script_name,line,tk[idx]);bad=1;continue;}
      if(!add_cmd(C_AUTOGOTO,line,arg,0.0)){bad=1;continue;}
      cmds[ncmds-1].ag_kind=kind;
      cmds[ncmds-1].ag_n=nn;
      cmds[ncmds-1].ag_label=lab;
      cmds[ncmds-1].ag_branch=br;
      strcpy(cmds[ncmds-1].ag_sym,sym);
    }
    else if(strcmp(word,"shot")==0){
      const char *ap=arg;
      char t[300],file[300];
      int sauto=0,svg=0,okp=1,ocolor=0,oport=0,ofs=0;
      unsigned int mask=0;
      double olw=0.0;
      char ofont[24]="";
      const char *dot;
      if(!tok_next(&ap,t,sizeof(t))){
        load_error("%s:%d: shot needs a file name, e.g. shot fig1.ps   or   shot auto diagram.svg",script_name,line);
        bad=1;continue;
      }
      if(strcasecmp(t,"auto")==0){
        sauto=1;
        if(!tok_next(&ap,t,sizeof(t))){
          load_error("%s:%d: shot auto needs a file name",script_name,line);
          bad=1;continue;
        }
      }
      strncpy(file,t,sizeof(file)-1);file[sizeof(file)-1]=0;
      dot=strrchr(file,'.');
      if(dot==NULL||(strcasecmp(dot,".ps")!=0&&strcasecmp(dot,".svg")!=0)){
        load_error("%s:%d: shot: the file name must end in .ps (Postscript) or .svg",script_name,line);
        bad=1;continue;
      }
      svg=(strcasecmp(dot,".svg")==0);
      while(tok_next(&ap,t,sizeof(t))){
        char *eq=strchr(t,'='),*key,*val,*e;
        int k;
        if(eq==NULL||eq==t||eq[1]==0){
          load_error("%s:%d: shot: \"%s\" is not an option (name=value: color orient fontsize font linewidth)",script_name,line,t);
          okp=0;break;
        }
        *eq=0;key=t;val=eq+1;
        for(k=0;key[k];k++)key[k]=(char)tolower((unsigned char)key[k]);
        if(strcmp(key,"color")==0){
          if(strcmp(val,"0")!=0&&strcmp(val,"1")!=0){load_error("%s:%d: shot: color must be 0 (black and white) or 1",script_name,line);okp=0;break;}
          ocolor=atoi(val);mask|=1;
        }
        else if(strcmp(key,"orient")==0){
          if(strcasecmp(val,"land")==0||strcasecmp(val,"landscape")==0)oport=0;
          else if(strcasecmp(val,"port")==0||strcasecmp(val,"portrait")==0)oport=1;
          else {load_error("%s:%d: shot: orient must be landscape or portrait",script_name,line);okp=0;break;}
          mask|=2;
        }
        else if(strcmp(key,"fontsize")==0){
          ofs=(int)strtol(val,&e,10);
          if(*e||ofs<1){load_error("%s:%d: shot: fontsize must be a positive integer",script_name,line);okp=0;break;}
          mask|=4;
        }
        else if(strcmp(key,"font")==0){
          if(strlen(val)>=sizeof(ofont)){load_error("%s:%d: shot: font name too long",script_name,line);okp=0;break;}
          strcpy(ofont,val);mask|=8;
        }
        else if(strcmp(key,"linewidth")==0){
          olw=strtod(val,&e);
          if(*e||olw<=0.0){load_error("%s:%d: shot: linewidth must be a positive number",script_name,line);okp=0;break;}
          mask|=16;
        }
        else {
          load_error("%s:%d: shot: unknown option \"%s\" (color orient fontsize font linewidth)",script_name,line,key);
          okp=0;break;
        }
      }
      if(!okp){bad=1;continue;}
      if(!add_cmd(C_SHOT,line,file,0.0)){bad=1;continue;}
      cmds[ncmds-1].shot_auto=sauto;
      cmds[ncmds-1].shot_svg=svg;
      cmds[ncmds-1].o_mask=mask;
      cmds[ncmds-1].o_color=ocolor;
      cmds[ncmds-1].o_port=oport;
      cmds[ncmds-1].o_fs=ofs;
      cmds[ncmds-1].o_lw=olw;
      strcpy(cmds[ncmds-1].o_font,ofont);
    }
    else if(strcmp(word,"try")==0||strcmp(word,"else")==0||strcmp(word,"end")==0){
      if(*arg){
        load_error("%s:%d: \"%s\" takes no arguments",script_name,line,word);
        bad=1;continue;
      }
      if(word[0]=='t'){                                    /* try */
        if(psp>=MAX_TRY){
          load_error("%s:%d: try blocks are nested more than %d deep",script_name,line,MAX_TRY);
          bad=1;continue;
        }
        if(!add_cmd(C_TRY,line,"",0.0)){bad=1;continue;}
        pstack[psp]=ncmds-1;
        pelse[psp]=-1;
        psp++;
      }
      else if(word[1]=='l'){                               /* else */
        if(psp==0||pelse[psp-1]>=0){
          load_error("%s:%d: \"else\" without a matching try",script_name,line);
          bad=1;continue;
        }
        if(!add_cmd(C_ELSE,line,"",0.0)){bad=1;continue;}
        pelse[psp-1]=ncmds-1;
        cmds[ncmds-1].link=pstack[psp-1];
      }
      else {                                               /* end */
        int t,e;
        if(psp==0){
          load_error("%s:%d: \"end\" without a matching try",script_name,line);
          bad=1;continue;
        }
        psp--;
        t=pstack[psp];e=pelse[psp];
        if(!add_cmd(C_END,line,"",0.0)){bad=1;continue;}
        cmds[ncmds-1].link=t;
        /* on a failure the try jumps to the start of the else part, or to its end */
        cmds[t].jump=(e>=0)?e+1:ncmds-1;
        /* when the try part finishes normally the else part is skipped */
        if(e>=0)cmds[e].jump=ncmds;
      }
    }
    else if(strcmp(word,"ok")==0||strcmp(word,"cancel")==0){
      if(*arg!=0){load_error("%s:%d: %s takes no argument",script_name,line,word);bad=1;continue;}
      add_cmd(word[0]=='o'?C_OK:C_CANCEL,line,"",0.0);
    }
    else if(strcmp(word,"ans")==0){
      if(*arg==0){load_error("%s:%d: ans needs an answer (use \"\" for an empty one)",script_name,line);bad=1;continue;}
      unquote(arg);
      add_cmd(C_ANS,line,arg,0.0);
    }
    else {
      load_error("%s:%d: unknown command \"%s\"",script_name,line,word);
      bad=1;
    }
  }
  if(psp>0){
    load_error("%s: a try block has no matching \"end\"",script_name);
    bad=1;
  }
  cur_ctx=NULL;
  free_expanded();
  if(bad){
    free_cmds();
    return 1;
  }
  return 0;
}

/* Called once the window exists and just before the main loop starts. */
void script_record_from_env(void);
void script_start(void)
{
  script_record_from_env();
  if(script_name[0]==0)return;
  if(script_load()!=0){
    exit(2);
  }
  pc=0;
  try_sp=0;caught_total=0;recovering=0;last_error[0]=0;
  active=1;
  if(opt_nopause)
    printf("script %s: %d commands (-nopause)\n",script_name,ncmds);
  else
    printf("script %s: %d commands\n",script_name,ncmds);
  fflush(stdout);
}


/* File -> Do script: pick a script file and run it inside the running program.
   Pauses and waits are honoured (a person is at the keyboard) and errors are
   shown in a box.  Scripts cannot start other scripts. */
extern char this_file[];
int file_selector(char *title,char *file,char *wild);

/* default script name: the model name without .ode, plus .xps (model.ode -> model.xps) */
static void default_script_name(char *buf,int n)
{
  char tmp[XPP_MAX_NAME];
  int len;
  strncpy(tmp,this_file,sizeof(tmp)-1);tmp[sizeof(tmp)-1]=0;
  len=(int)strlen(tmp);
  if(len>4&&strcasecmp(tmp+len-4,".ode")==0)tmp[len-4]=0;
  snprintf(buf,n,"%s.xps",tmp);
}

void script_run_from_menu(void)
{
  char filename[XPP_MAX_NAME];
  if(active){
    err_msg("A script is already running -- a script cannot start another script");
    return;
  }
  default_script_name(filename,sizeof(filename));
  if(!file_selector("Run script",filename,"*.xps"))return;
  script_set_file(filename);
  opt_nopause=0;
  opt_exit=0;
  if(script_load()!=0){
    char box[460];
    snprintf(box,sizeof(box),"%s",load_err_text[0]?load_err_text:"The script could not be read");
    script_name[0]=0;
    err_msg(box);
    return;
  }
  pc=0;
  try_sp=0;caught_total=0;recovering=0;last_error[0]=0;
  key_pos=0;
  paused=0;
  waiting=0;
  wait_prompt_done=0;
  key_gap_until=0.0;
  key_delay=0.0;
  active=1;
  printf("script %s: %d commands\n",script_name,ncmds);
  fflush(stdout);
}

/* ------------------------------------------------------------- execution */

/* After a "set": do what XPP does when the same things are changed by hand,
   so its internal state agrees with what was set.
     - redo_stuff(): derived parameters, convolution kernels and tabulated
       functions are recomputed from the new parameter values (the Parameters
       dialog does exactly this);
     - the delay buffer is rebuilt only if the maximal delay or the step size
       changed (rebuilding it also discards the last solution, so doing it
       on every set would break "Initialconds -> Last" and "Continue");
     - for 2D plots XPP keeps xmin,xmax,ymin,ymax equal to xlo,xhi,ylo,yhi (the
       Window dialog, zoom and the loader all do this, and the 2D View dialog
       reads xmin..ymax), so the window is copied there as well, or XPP
       would show -- and later restore -- the old window;
     - the parameter and IC boxes and sliders are redrawn, the graph's
       time/3D flags are refreshed (check_flags) and the graph reset. */
static void graph_commit(int delay_changed)
{
  if(AXES<5){
    x_3d[0]=MY_XLO;
    x_3d[1]=MY_XHI;
    y_3d[0]=MY_YLO;
    y_3d[1]=MY_YHI;
  }
  if(delay_changed)chk_delay();
  redo_stuff();
  redraw_params();
  redraw_ics();
  /* reset_graph() copies the plotted variables into MyGraph and redraws, but
     it does not refresh the flags derived from them: TimeFlag (plotting
     against time, which switches off nullclines, direction fields, ...) and
     ThreeDFlag.  XPP's own dialogs call check_flags() for that.  Do it first,
     so the redraw inside reset_graph() already uses the right flags.       */
  MyGraph->xv[0]=IXPLT;
  MyGraph->yv[0]=IYPLT;
  MyGraph->zv[0]=IZPLT;
  MyGraph->grtype=AXES;
  check_flags();
  reset_graph();
}

static void do_set(SCMD *c)
{
  char *tok,tmp[2048];
  char name[256],value[256],*eq;
  int n,old_maxstor=MAXSTOR,dchg;
  double old_delay=DELAY,old_dt=DELTA_T;
  strncpy(tmp,c->arg,sizeof(tmp)-1);tmp[sizeof(tmp)-1]=0;
  get_graph();
  /* same tokenizer as extract_action(): separators are blanks , ; */
  for(tok=strtok(tmp," ,;\n");tok!=NULL;tok=strtok(NULL," ,;\n")){
    eq=strchr(tok,'=');
    if(eq==NULL||eq==tok||eq[1]==0){
      script_fail(c->line,"bad assignment \"%s\" (expected name=value)",tok);
      graph_commit(DELAY!=old_delay||DELTA_T!=old_dt);
      return;
    }
    n=eq-tok;
    if(n>255)n=255;
    strncpy(name,tok,n);name[n]=0;
    strncpy(value,eq+1,255);value[255]=0;
    set_option_unmatched=0;
    do_intern_set(name,value);
    if(set_option_unmatched){
      script_fail(c->line,"unknown parameter, variable or option \"%s\"",name);
      graph_commit(DELAY!=old_delay||DELTA_T!=old_dt);
      return;
    }
    if(MAXSTOR!=old_maxstor){
      /* the storage arrays are allocated once, at start-up */
      MAXSTOR=old_maxstor;
      script_fail(c->line,"maxstor cannot be changed from a script (it is fixed when XPP starts; use @ maxstor=... in the .ode file)");
      graph_commit(DELAY!=old_delay||DELTA_T!=old_dt);
      return;
    }
    if(strcasecmp(name,"total")==0){
      /* as in the Numerics menu: a negative total means run until stopped */
      if(TEND<0.0){FOREVER=1;TEND=-TEND;}
      else FOREVER=0;
    }
  }
  dchg=(DELAY!=old_delay||DELTA_T!=old_dt);
  graph_commit(dchg);
}

static void do_print(SCMD *c)
{
  char name[80],up[80];
  double v;
  int i;
  strncpy(name,c->arg,79);name[79]=0;
  rtrim(name);
  convert(name,up);
  i=find_user_name(XS_IC,up);
  if(i>-1){
    /* a state variable: report its initial condition (what "set x=..." sets) */
    printf("%s(0) = %.12g\n",name,last_ic[i]);
    fflush(stdout);
    return;
  }
  if(get_val(up,&v)){
    printf("%s = %.12g\n",name,v);
    fflush(stdout);
  }
  else
    script_fail(c->line,"print: \"%s\" is not a parameter, variable or initial condition",name);
}

/* Build the KeyPress event for one key. */
static void fill_key_event(XEvent *ev,Window win,SKEY *k)
{
  XKeyEvent *e;
  memset(ev,0,sizeof(XEvent));
  e=&ev->xkey;
  e->type=KeyPress;
  e->send_event=True;
  e->display=display;
  e->window=win;
  e->root=DefaultRootWindow(display);
  e->subwindow=None;
  e->time=CurrentTime;
  e->x=e->y=e->x_root=e->y_root=1;
  e->state=k->state;
  e->keycode=k->kc;
  e->same_screen=True;
}


static void fill_button_event(XEvent *ev,Window win)
{
  XButtonEvent *e;
  memset(ev,0,sizeof(XEvent));
  e=&ev->xbutton;
  e->type=ButtonPress;
  e->send_event=True;
  e->display=display;
  e->window=win;
  e->root=DefaultRootWindow(display);
  e->subwindow=None;
  e->time=CurrentTime;
  e->x=e->y=e->x_root=e->y_root=1;
  e->button=Button1;
  e->same_screen=True;
}

/* ------------------------------------------------------------------ shot
   Save the picture straight to a file, without going through the dialogs
   (so a script does not depend on which options forms a menu happens to
   have).  The Postscript options given on the line are used for this one file
   only; the settings are put back afterwards.                              */
extern int PS_Color,PS_Port,PS_FONTSIZE;
extern char PS_FONT[];
extern double PS_LW;
int ps_init(char *filename,int color);
int svg_init(char *filename,int color);
void ps_restore(void);
void svg_restore(void);
int auto_diagram_to_file(const char *filename,int svg);

static void do_shot(SCMD *c)
{
  int oc=PS_Color,op=PS_Port,of=PS_FONTSIZE,r=1;
  double ol=PS_LW;
  char ofnt[100];
  strncpy(ofnt,PS_FONT,sizeof(ofnt)-1);ofnt[sizeof(ofnt)-1]=0;
  if(c->o_mask&1)PS_Color=c->o_color;
  if(c->o_mask&2)PS_Port=c->o_port;
  if(c->o_mask&4)PS_FONTSIZE=c->o_fs;
  if(c->o_mask&8){strncpy(PS_FONT,c->o_font,99);PS_FONT[99]=0;}
  if(c->o_mask&16)PS_LW=c->o_lw;
  if(c->shot_auto)
    r=auto_diagram_to_file(c->arg,c->shot_svg);
  else if(c->shot_svg){
    if(svg_init(c->arg,PS_Color))svg_restore();
    else r=-1;
  }
  else {
    if(ps_init(c->arg,PS_Color))ps_restore();
    else r=-1;
  }
  PS_Color=oc;PS_Port=op;PS_FONTSIZE=of;PS_LW=ol;
  strncpy(PS_FONT,ofnt,99);PS_FONT[99]=0;
  if(r==0)script_fail(c->line,"shot auto: there is no AUTO diagram to save");
  else if(r<0)script_fail(c->line,"shot: cannot write %s",c->arg);
  else {
    printf("script: saved %s\n",c->arg);
    fflush(stdout);
  }
}

/* Run the next command if it can run now.  Returns:
     0  nothing to do (let XPP block for events),
     1  a command ran; look again immediately,
     2  timed wait; *ms is the timeout in milliseconds,
     3  *ev holds a synthetic key press; hand it to the caller.           */
int script_poll(int *ms,XEvent *ev)
{
  SCMD *c;
  double t;
  Window win;
  if(!active)return 0;
  if(recovering){
    /* a failure was caught: close whatever menu or dialog the failed command left
       open, by pressing Esc, until XPP is waiting at the top level again */
    if(at_top)recovering=0;
    else if(recover_count>=9){
      recovering=0;
      try_sp=0;
      script_fail(cmds[pc>0&&pc<=ncmds?pc-1:0].line,"could not return to the top level after a failure (a menu or dialog is still open)");
      return 1;
    }
    else {
      /* Esc closes menus, prompts, the file selector and message boxes, but a
         multi-field form has no keyboard cancel: there Esc only moves to the next
         field, and Tab (Enter with MSStyle) accepts the values as they are.  So:
         Esc three times, then the form's accept key twice, then Esc again for
         whatever that opened (e.g. the file selector after the Postscript form). */
      SKEY k;
      KeySym want=XK_Escape;
      if(recover_count==3||recover_count==4)want=MSStyle?XK_Return:XK_Tab;
      if(!resolve_keysym(want,&k)){recovering=0;try_sp=0;return 1;}
      recover_count++;
      fill_key_event(ev,main_win,&k);
      return 3;
    }
  }
  if(paused)return 0;
  if(waiting){
    t=wait_until-now_sec();
    if(t>0.0){*ms=(int)ceil(t*1000.0);return 2;}
    waiting=0;
  }
  if(key_gap_until>0.0){
    t=key_gap_until-now_sec();
    if(t>0.0){*ms=(int)ceil(t*1000.0);return 2;}
    key_gap_until=0.0;
  }
  if(wait_prompt_done){        /* a person is answering a prompt ("ans ?") */
    if(!at_top)return 0;
    wait_prompt_done=0;
  }
  if(pc>=ncmds){
    if(!at_top){
      /* a menu or dialog is still open and nothing is left to close it */
      if(opt_nopause){
        script_fail(cmds[ncmds-1].line,"the script ended while a menu or dialog is still open (missing keys or ans?)");
        return 1;
      }
      return 0;
    }
    if(caught_total>0){
      printf("script %s: finished, but %d failure%s caught by try (exit status 3 with -scriptexit)\n",
             script_name,caught_total,caught_total==1?" was":"s were");
      fflush(stdout);
      script_finish(3);
    }
    else {
      printf("script %s: finished\n",script_name);
      fflush(stdout);
      script_finish(0);
    }
    return 1;
  }
  c=&cmds[pc];
  run_ctx=c->ctx;
  switch(c->type){
  case C_SET:
  case C_PRINT:
  case C_AUTOGOTO:
  case C_AUTOLIST:
  case C_SHOT:
    if(!at_top){
      if(opt_nopause){
        script_fail(c->line,"\"%s\" cannot run while a menu or dialog is open (missing keys or ans?)",
                    c->type==C_SET?"set":(c->type==C_PRINT?"print":(c->type==C_SHOT?"shot":"auto")));
        return 1;
      }
      return 0;
    }
    pc++;
    if(c->type==C_SET)do_set(c);
    else if(c->type==C_PRINT)do_print(c);
    else if(c->type==C_SHOT)do_shot(c);
    else if(c->type==C_AUTOLIST)auto_list_points();
    else {
      char err[300];
      if(auto_goto_point(c->ag_kind,c->ag_sym,c->ag_n,c->ag_label,c->ag_branch,err,sizeof(err)))
        script_fail(c->line,"auto goto: %s",err);
    }
    return 1;
  case C_BROWSER:
  case C_APLOT:
    {
      /* a button of the data browser / array plot window is pressed */
      Window bw=0;
      pc++;
      if(c->type==C_BROWSER&&strcmp(c->arg,"open")==0)bw=top_button_window(5);
      else if(c->type==C_BROWSER)bw=browser_button_window(c->arg);
      else bw=aplot_button_window(c->arg);
      if(bw==0){
        script_fail(c->line,"%s %s: the %s window is not open",c->type==C_BROWSER?"browser":"aplot",c->arg,
                    c->type==C_BROWSER?"data browser":"array plot");
        return 1;
      }
      fill_button_event(ev,bw);
      if(!opt_nopause&&key_delay>0.0)key_gap_until=now_sec()+key_delay;
      return 3;
    }
  case C_KEYS:
  case C_AUTOKEYS:
  case C_BROWSEKEYS:
    if(c->type==C_BROWSEKEYS){
      win=browser_key_window();
      if(win==0){
        script_fail(c->line,"browsekeys: the data browser is not open (open it with: browser open)");
        return 1;
      }
    }
    else if(c->type==C_AUTOKEYS){
      win=auto_key_window();
      if(win==0){
        script_fail(c->line,"autokeys: the AUTO window is not open (open it with keys f a first)");
        return 1;
      }
    }
    else win=main_win;
    fill_key_event(ev,win,&c->keys[key_pos]);
    key_pos++;
    if(key_pos>=c->nkeys){key_pos=0;pc++;}
    if(!opt_nopause&&key_delay>0.0)key_gap_until=now_sec()+key_delay;
    return 3;
  case C_ANS:
  case C_OK:
  case C_CANCEL:
    pc++;
    script_fail(c->line,"%s has no prompt to answer (the line before it must open one)",
                c->type==C_ANS?"ans":(c->type==C_OK?"ok":"cancel"));
    return 1;
  case C_CLICK:
  case C_DRAG:
    pc++;
    script_fail(c->line,"%s has nothing to answer: the line before it must start something that waits for the mouse (for example keys t a, then ans for the text)",
                c->type==C_CLICK?"click":"drag");
    return 1;
  case C_MSG:
    {
      char text[1200];
      pc++;
      expand_error(c->arg,text,sizeof(text));
      printf("%s\n",text);
      fflush(stdout);
      if(!opt_nopause&&Xup)bottom_msg(0,text);
    }
    return 1;
  case C_PAUSE:
    {
      char text[1200];
      pc++;
      if(opt_nopause)return 1;
      expand_error(c->arg,text,sizeof(text));
      if(Xup)bottom_msg(0,text[0]?text:"Paused -- press Enter to continue (Esc stops the script)");
      paused=1;
    }
    return 1;
  case C_TRY:
    pc++;
    if(try_sp>=MAX_TRY){
      script_fail(c->line,"try blocks are nested more than %d deep",MAX_TRY);
      return 1;
    }
    try_stack[try_sp++]=(int)(c-cmds);
    return 1;
  case C_ELSE:
    /* the try part finished without a failure: leave it and skip the else part */
    if(try_sp>0&&try_stack[try_sp-1]==c->link)try_sp--;
    pc=c->jump;
    return 1;
  case C_END:
    pc++;
    if(try_sp>0&&try_stack[try_sp-1]==c->link)try_sp--;
    return 1;
  case C_WAIT:
    pc++;
    if(opt_nopause)return 1;
    wait_until=now_sec()+c->num;
    waiting=1;
    return 1;
  case C_SPEED:
    pc++;
    key_delay=c->num;
    return 1;
  case C_QUIT:
    pc++;
    printf("script %s: quit\n",script_name);
    if(caught_total>0)printf("script %s: %d failure%s caught by try\n",script_name,caught_total,caught_total==1?" was":"s were");
    fflush(stdout);
    yes_reset_auto();
    exit(caught_total>0?3:0);
  }
  return 1;
}

/* ----------------------------------------------------------- prompt answers
   A prompt (see ggets.c, pop_list.c, dialog_box.c) asks for its answer
   before it opens.  If the next script command is "ans" it is consumed here.
   Return values:  1 answer supplied,  0 not scripted -- open the prompt as
   usual (it can then be driven by following keys lines),  -1 script error,
   the prompt should behave as if cancelled.                                 */

static SCMD *next_ans(void)
{
  if(!active||pc>=ncmds)return NULL;
  if(cmds[pc].type!=C_ANS&&cmds[pc].type!=C_OK&&cmds[pc].type!=C_CANCEL)return NULL;
  run_ctx=cmds[pc].ctx;
  return &cmds[pc++];
}

/* "ans ?": hand this prompt to the person */
static int ans_is_ask(SCMD *c)
{
  if(c->type!=C_ANS||strcmp(c->arg,"?")!=0)return 0;
  if(opt_nopause){
    script_fail_quiet(c->line,"\"ans ?\" needs a person at the keyboard (not allowed with -nopause or -scriptexit)");
    return -1;
  }
  wait_prompt_done=1;
  return 1;
}

int script_take_text(const char *label,char *buf,int cap)
{
  SCMD *c=next_ans();
  int r;
  (void)label;
  if(c==NULL)return 0;
  if(c->type==C_CANCEL)return -1;
  if(c->type==C_OK)return 1;          /* keep the default */
  r=ans_is_ask(c);
  if(r==-1)return -1;
  if(r==1)return 0;
  if((int)strlen(c->arg)>=cap){
    script_fail_quiet(c->line,"the answer is too long for this prompt (at most %d characters)",cap-1);
    return -1;
  }
  strcpy(buf,c->arg);
  return 1;
}

/* two-way choice: key holds the two letters, e.g. "yn" */
int script_take_choice(const char *prompt,const char *key,int *index)
{
  SCMD *c=next_ans();
  int r;
  char ch;
  if(c==NULL)return 0;
  if(c->type==C_CANCEL)return -1;
  if(c->type==C_OK){
    script_fail_quiet(c->line,"\"ok\" does not apply to \"%s\"; answer it with ans %c or ans %c",prompt,key[0],key[1]);
    return -1;
  }
  r=ans_is_ask(c);
  if(r==-1)return -1;
  if(r==1)return 0;
  ch=(char)tolower((unsigned char)c->arg[0]);
  if(ch!=0&&ch==tolower((unsigned char)key[0])){*index=0;return 1;}
  if(ch!=0&&ch==tolower((unsigned char)key[1])){*index=1;return 1;}
  script_fail_quiet(c->line,"the answer to \"%s\" must be %c or %c",prompt,key[0],key[1]);
  return -1;
}

static void trim_both(char *s)
{
  int n;
  char *p=s;
  while(*p&&isspace((unsigned char)*p))p++;
  if(p!=s)memmove(s,p,strlen(p)+1);
  n=strlen(s);
  while(n>0&&isspace((unsigned char)s[n-1]))s[--n]=0;
}

/* compare a field label ("*1Xmin:" -> "xmin") with a name from the ans line */
static int label_matches(const char *label,const char *name)
{
  char a[80],b[80];
  int i,n;
  const char *p=label;
  if(*p=='*'&&isdigit((unsigned char)p[1])){p+=2;}   /* "*N" hint prefix */
  strncpy(a,p,79);a[79]=0;
  strncpy(b,name,79);b[79]=0;
  trim_both(a);trim_both(b);
  n=strlen(a);if(n>0&&a[n-1]==':')a[n-1]=0;
  n=strlen(b);if(n>0&&b[n-1]==':')b[n-1]=0;
  for(i=0;a[i];i++)a[i]=(char)tolower((unsigned char)a[i]);
  for(i=0;b[i];i++)b[i]=(char)tolower((unsigned char)b[i]);
  return strcmp(a,b)==0;
}

/* multi-field form.  vals is an array of n strings of stride bytes. */
int script_take_form(const char *title,int n,char **names,char *vals,int stride)
{
  SCMD *c=next_ans();
  int r,j,m=0,target,lit;
  char field[256],*p,*q,*eq,tmp[1024];
  char name[80],*val;
  if(c==NULL)return 0;
  if(c->type==C_CANCEL)return -1;
  if(c->type==C_OK)return 1;          /* press OK with the values as shown */
  r=ans_is_ask(c);
  if(r==-1)return -1;
  if(r==1)return 0;
  strncpy(tmp,c->arg,sizeof(tmp)-1);tmp[sizeof(tmp)-1]=0;
  p=tmp;
  for(;;){
    int len;
    q=strchr(p,'|');
    if(q)*q=0;
    strncpy(field,p,255);field[255]=0;
    trim_both(field);
    /* a quoted field is a literal value (it may contain = or be empty) */
    lit=0;
    len=strlen(field);
    if(len>=2&&(field[0]=='"'||field[0]=='\'')&&field[len-1]==field[0]){
      memmove(field,field+1,len-2);
      field[len-2]=0;
      lit=1;
    }
    target=m;
    val=field;
    eq=lit?NULL:strchr(field,'=');
    if(eq!=NULL){
      int nl=(int)(eq-field);
      if(nl>79)nl=79;
      strncpy(name,field,nl);
      name[nl]=0;
      for(j=0;j<n;j++)if(label_matches(names[j],name))break;
      if(j>=n){
        /* "Label=value" but there is no such field: say so (and list the fields)
           instead of silently storing the text in the wrong field. */
        char list[400]="";
        int k,used=0;
        for(k=0;k<n;k++){
          const char *lb=names[k];
          if(lb[0]=='*'&&isdigit((unsigned char)lb[1]))lb+=2;
          used+=snprintf(list+used,sizeof(list)-used,"%s%s",k?", ":"",lb);
          if(used>=(int)sizeof(list)-40){strcat(list,", ...");break;}
        }
        script_fail_quiet(c->line,"\"%s\" has no field called \"%s\" (its fields are: %s; quote the value to use a literal =)",title,name,list);
        return -1;
      }
      target=j;
      val=eq+1;
      trim_both(val);
    }
    if(target>=n){
      script_fail_quiet(c->line,"\"%s\" has only %d fields",title,n);
      return -1;
    }
    if(!lit&&val==field&&field[0]==0){
      /* empty positional field: keep the current value */
    }
    else {
      if((int)strlen(val)>=stride){
        script_fail_quiet(c->line,"field %d of \"%s\" is too long (at most %d characters)",target+1,title,stride-1);
        return -1;
      }
      strcpy(vals+target*stride,val);
    }
    m++;
    if(q==NULL)break;
    p=q+1;
  }
  return 1;
}

/* An error box is about to open.  While a script runs this means something
   went wrong: stop the script (and exit when unattended).  Returns 1 if the
   box should be suppressed.                                                */

/* ------------------------------------------------------------ mouse waits
   A few commands wait for the mouse: placing text, a label or a freeze key, a start point with
   the mouse, the zoom rectangle, an arrow.  They wait in a loop of their own, which the script
   engine never sees, so they ask for a scripted answer first, like the prompts above.

     get_mouse_xy()  -> script_mouse_point()      the line after the command is   click x y
     rubber()        -> script_mouse_drag()       the line after the command is   drag x1 y1 x2 y2

   x, y are data coordinates of the current plot (click pixel x y: window pixels).  "click ?" and
   "drag ?" leave that one click to the person at the keyboard.  With nothing to answer the wait,
   an unattended script (-nopause, -scriptexit) stops with an error that names the line (instead of
   waiting for ever); an interactive one waits for the person.
   A following "keys esc" ends the wait, as Esc does for a person.  A wait that repeats
   (Initialconds -> mIce) takes one click line per click and then "keys esc".
   Return values:  1 supplied,  0 not scripted -- wait for the mouse as usual,  -1 script error or
   Esc, act as if the person had pressed Esc.                                               */
extern Window draw_win;
void get_draw_area(void);
extern int DLeft,DRight,DTop,DBottom;

static int mouse_to_pixels(SCMD *c,int k,int *ix,int *iy)
{
  if(c->mpix){
    *ix=(int)floor(c->mv[k]+0.5);
    *iy=(int)floor(c->mv[k+1]+0.5);
  }
  else {
    /* the inverse of scale_to_real(), which turns a real click into data coordinates */
    double dx,dy;
    get_draw_area();
    dx=MyGraph->xhi-MyGraph->xlo;
    dy=MyGraph->yhi-MyGraph->ylo;
    if(dx==0.0||dy==0.0){
      script_fail(c->line,"%s: the plot window has no extent (xlo=xhi or ylo=yhi)",
                  c->type==C_CLICK?"click":"drag");
      return 0;
    }
    *ix=(int)floor((c->mv[k]-MyGraph->xlo)*(double)(DRight-DLeft)/dx+DLeft+0.5);
    *iy=(int)floor((c->mv[k+1]-MyGraph->ylo)*(double)(DTop-DBottom)/dy+DBottom+0.5);
  }
  {
    /* the point must be inside the plot window, as a real click has to be */
    Window root;
    int wx,wy;
    unsigned int ww,wh,bw,dp;
    if(Xup&&XGetGeometry(display,draw_win,&root,&wx,&wy,&ww,&wh,&bw,&dp)){
      if(*ix<0||*iy<0||*ix>=(int)ww||*iy>=(int)wh){
        script_fail(c->line,"%s %g %g is outside the plot window (window pixels %d,%d; the window is %u by %u)",
                    c->type==C_CLICK?"click":"drag",c->mv[k],c->mv[k+1],*ix,*iy,ww,wh);
        return 0;
      }
    }
  }
  return 1;
}

static int mouse_request(int want,int *x,int *y,int *x2,int *y2)
{
  SCMD *c;
  int line;
  const char *what=want==C_CLICK?"a mouse click":"a mouse drag (a rectangle or a line)";
  if(!active)return 0;
  line=(pc>0&&pc<=ncmds)?cmds[pc-1].line:(ncmds>0?cmds[ncmds-1].line:0);
  if(pc<ncmds&&(cmds[pc].type==C_CLICK||cmds[pc].type==C_DRAG)){
    c=&cmds[pc];
    run_ctx=c->ctx;
    if(c->type!=want){
      pc++;
      script_fail(c->line,"XPP is waiting for %s here, but this line is a %s: use %s",what,
                  c->type==C_CLICK?"click":"drag",want==C_CLICK?"click x y":"drag x1 y1 x2 y2");
      return -1;
    }
    pc++;
    if(c->mq){                         /* "click ?": the person does it */
      if(opt_nopause){
        script_fail(c->line,"\"%s ?\" leaves the mouse to a person, which is not possible with -nopause or -scriptexit",
                    c->type==C_CLICK?"click":"drag");
        return -1;
      }
      return 0;
    }
    if(!mouse_to_pixels(c,0,x,y))return -1;
    if(want==C_DRAG&&!mouse_to_pixels(c,2,x2,y2))return -1;
    return 1;
  }
  if(pc<ncmds&&cmds[pc].type==C_KEYS&&key_pos<cmds[pc].nkeys){
    /* the next key is Esc: it ends the wait, as it does for a person (a wait that repeats, such as
       Initialconds -> mIce, is ended this way after its last click) */
    SKEY esc;
    if(resolve_keysym(XK_Escape,&esc)&&cmds[pc].keys[key_pos].kc==esc.kc){
      key_pos++;
      if(key_pos>=cmds[pc].nkeys){key_pos=0;pc++;}
      return -1;
    }
  }
  if(opt_nopause){
    if(pc>0&&pc<=ncmds)run_ctx=cmds[pc-1].ctx;
    script_fail(line,"XPP is waiting for %s and the script has nothing to give it: add a \"%s\" line after this one, or \"%s ?\" to leave it to a person",
                what,want==C_CLICK?"click x y":"drag x1 y1 x2 y2",want==C_CLICK?"click":"drag");
    return -1;
  }
  return 0;                            /* interactive: the person does it */
}

int script_mouse_point(int *x,int *y)
{
  return mouse_request(C_CLICK,x,y,NULL,NULL);
}

int script_mouse_drag(int *x1,int *y1,int *x2,int *y2)
{
  return mouse_request(C_DRAG,x1,y1,x2,y2);
}

int script_message(const char *msg)
{
  int line;
  char cx[300],full[700];
  if(!active)return 0;
  line=(pc>0&&pc<=ncmds)?cmds[pc-1].line:0;
  ctx_text(cx,sizeof(cx),(pc>0&&pc<=ncmds)?cmds[pc-1].ctx:NULL);
  snprintf(full,sizeof(full),"%s:%d: XPP reported: %s%s",script_name,line,msg,cx);
  if(catch_failure(full,msg))return 1;       /* inside try: the error box is not shown */
  fprintf(stderr,"script %s: stopped at line %d%s because XPP reported: %s\n",script_name,line,cx,msg);
  fflush(stderr);
  active=0;paused=0;waiting=0;
  if(opt_exit){yes_reset_auto();exit(1);}
  return 0;   /* interactive: let the person see the box */
}

/* Look at every event that xpp_next_event() is about to return.  While the
   script is paused, Enter/Space ends the pause (the key is swallowed) and Esc
   stops the script.  Returns 1 if the event was consumed.                  */
int script_filter_event(XEvent *ev)
{
  KeySym ks;
  if(!active||!paused||!at_top)return 0;
  if(ev->type!=KeyPress)return 0;
  ks=XLookupKeysym(&ev->xkey,0);
  if(ks==XK_Return||ks==XK_KP_Enter||ks==XK_space||ks==XK_Linefeed){
    paused=0;
    if(Xup)bottom_msg(0,"");
    return 1;
  }
  if(ks==XK_Escape){
    {char cx[300];ctx_text(cx,sizeof(cx),pc>0?cmds[pc-1].ctx:NULL);
     printf("script %s: stopped by user at line %d%s\n",script_name,pc>0?cmds[pc-1].line:0,cx);}
    fflush(stdout);
    active=0;paused=0;waiting=0;
    if(Xup)bottom_msg(0,"Script stopped");
    if(opt_exit){yes_reset_auto();exit(1);}
    return 1;
  }
  return 0;
}


/* =========================================================================
   RECORDING  --  File -> Record script

   The recorder is the mirror image of the script engine above.  It writes
   the commands that, replayed, repeat what the person did:

     menu keys, typed or chosen with the mouse      keys i g
     keys typed in the AUTO window                   autokeys r
     the answer to a prompt (text, form, yes/no)     ans ... / ok / cancel
     a click or a drag in the plot window            click x y / drag x1 y1 x2 y2
     a parameter or initial condition that was edited  set a=1.3
     a point picked in the AUTO diagram               auto goto lab 5

   Keys typed into a prompt are NOT recorded as keys: the prompt is recorded
   as one "ans" line when it closes.  Keys typed into the parameter boxes are
   not recorded either; the value that was committed is (as a "set").
   Consecutive keys are collected into one line.  Nothing is recorded while a
   script runs (XPP_REC_TEST=1 in the environment lifts this; used for testing)
   and XPP_RECORD=file starts a recording at start-up.  XPP_REC_ECHO=1 echoes
   every recorded line to stderr.
   ========================================================================= */
#include <time.h>
#include <unistd.h>
extern char this_file[];
extern char upar_names[MAXPAR][11],uvar_names[MAXODE][12];
extern int NUPAR,NODE;
int TwoChoice(char *c1,char *c2,char *q,char *key);
int auto_owns_key(Window w);

#define REC_MAXDEPTH 8
static FILE *rec_fp=NULL;
static char rec_name[512]="";
static int rec_nokeys=0;        /* >0: keys are text being typed into a field */
static int rec_depth=0;         /* prompts that are open */
static int rec_test=0,rec_echo=0;
static int rec_kind=0;          /* pending key line: 0 none, 1 keys, 2 autokeys */
static char rec_buf[200];
static char *rec_init[REC_MAXDEPTH];    /* what each open prompt showed when it opened */
static int rec_traverse=0;

/* the last key event that was written (a key can pass through two readers) */
static unsigned long rec_ls=0;
static Time rec_lt=0;
static unsigned int rec_lk=0;
static Window rec_lw=0;

int script_recording(void){return rec_fp!=NULL;}
static int rec_on(void){return rec_fp!=NULL&&(!active||rec_test);}

static void rec_write(const char *line)
{
  if(rec_fp==NULL)return;
  fseek(rec_fp,0,SEEK_END);
  fputs(line,rec_fp);
  fputc('\n',rec_fp);
  fflush(rec_fp);
  if(rec_echo){fprintf(stderr,"rec: %s\n",line);fflush(stderr);}
}

static long rec_pos=0;          /* where the pending key line starts in the file */

/* The pending key line is kept on disk at all times (rewritten as it grows),
   so a crash or a kill loses nothing. */
static void rec_rewrite(void)
{
  if(rec_fp==NULL)return;
  fseek(rec_fp,rec_pos,SEEK_SET);
  if(rec_kind&&rec_buf[0]){fputs(rec_buf,rec_fp);fputc('\n',rec_fp);}
  fflush(rec_fp);
  if(ftruncate(fileno(rec_fp),(off_t)ftell(rec_fp))!=0){}
}

static void rec_flush(void)
{
  if(rec_kind&&rec_buf[0]&&rec_echo){fprintf(stderr,"rec: %s\n",rec_buf);fflush(stderr);}
  rec_kind=0;
  rec_buf[0]=0;
}

static void rec_token(int kind,const char *tok)
{
  const char *pre=(kind==1)?"keys":(kind==2?"autokeys":"browsekeys");
  if(rec_kind!=kind||strlen(rec_buf)+strlen(tok)+2>=76){   /* new line (they are kept short) */
    rec_flush();
    rec_kind=kind;
    fseek(rec_fp,0,SEEK_END);
    rec_pos=ftell(rec_fp);
    strcpy(rec_buf,pre);
  }
  strcat(rec_buf," ");
  strcat(rec_buf,tok);
  rec_rewrite();
}

/* the token that "keys" turns back into the character c */
static void char_token(int c,char *out)
{
  if(c==' ')strcpy(out,"space");
  else if(islower(c)||isdigit(c))sprintf(out,"%c",c);
  else if(c=='"')strcpy(out,"'\"'");
  else sprintf(out,"\"%c\"",c);     /* upper case and punctuation are typed literally */
}

void rec_nokey(int on)
{
  if(on)rec_nokeys++;
  else if(rec_nokeys>0)rec_nokeys--;
}

int get_key_press(XEvent *ev);
int get_key_press_nr(XEvent *ev)
{
  int r;
  rec_nokeys++;
  r=get_key_press(ev);
  rec_nokeys--;
  return r;
}

static int rec_force_kind=0;
static int rec_context=0;

/* The AUTO and data-browser key handlers say so when they take a key: whether a
   key reaches them depends on the window under the pointer and on the input
   focus, which only they know. */
void rec_key_kind(int kind){rec_force_kind=kind;if(kind==2)rec_context=2;}

/* While the AUTO window is handling a key or button, the keys its popup menus
   read (r s: Run, Steady state) belong to the AUTO window too. */
void rec_key_context(int kind){rec_context=kind;}

void rec_key_event(XEvent *ev)
{
  XKeyEvent *k;
  char buf[16],tok[16];
  KeySym ks;
  int n,kind;
  kind=rec_force_kind?rec_force_kind:(rec_context?rec_context:1);
  rec_force_kind=0;
  if(rec_fp==NULL||!rec_on()||rec_nokeys>0||rec_depth>0)return;
  if(ev->type!=KeyPress)return;
  k=&ev->xkey;
  if(k->serial!=0&&k->serial==rec_ls&&k->time==rec_lt&&k->keycode==rec_lk&&k->window==rec_lw)return;
  n=XLookupString(k,buf,sizeof(buf)-1,&ks,NULL);
  switch(ks){
  case XK_Escape: strcpy(tok,"esc");break;
  case XK_Return: case XK_KP_Enter: case XK_Linefeed: strcpy(tok,"enter");break;
  case XK_BackSpace: strcpy(tok,"bksp");break;
  case XK_Delete: strcpy(tok,"del");break;
  case XK_Tab: strcpy(tok,"tab");break;
  case XK_Up: strcpy(tok,"up");break;
  case XK_Down: strcpy(tok,"down");break;
  case XK_Left: strcpy(tok,"left");break;
  case XK_Right: strcpy(tok,"right");break;
  case XK_Home: strcpy(tok,"home");break;
  case XK_End: strcpy(tok,"end");break;
  case XK_Prior: strcpy(tok,"pgup");break;
  case XK_Next: strcpy(tok,"pgdn");break;
  default:
    if(n==1&&buf[0]>=' '&&buf[0]<='~')char_token((unsigned char)buf[0],tok);
    else return;                      /* Shift, Control, function keys, ... */
  }
  rec_ls=k->serial;rec_lt=k->time;rec_lk=k->keycode;rec_lw=k->window;
  rec_token(kind,tok);
}

/* a menu entry chosen with the mouse stands for its hot key */
void rec_menu_key(int ch)
{
  char tok[16];
  if(!rec_on()||rec_depth>0||rec_nokeys>0)return;
  if(ch==27)strcpy(tok,"esc");
  else if(ch==13)strcpy(tok,"enter");
  else if(ch>' '&&ch<='~')char_token(ch,tok);
  else return;
  rec_token(rec_context?rec_context:1,tok);
}

void rec_browse_key(const char *tok)
{
  if(!rec_on()||rec_depth>0)return;
  rec_token(3,tok);
}

void rec_command(const char *line)
{
  if(!rec_on()||rec_depth>0)return;
  rec_flush();
  rec_write(line);
}

void rec_comment(const char *text)
{
  char c[300];
  if(!rec_on()||rec_depth>0)return;
  rec_flush();
  snprintf(c,sizeof(c),"# %s",text);
  rec_write(c);
}

void rec_esc(void)
{
  if(!rec_on()||rec_depth>0)return;
  rec_token(1,"esc");
}

void rec_go(void)
{
  if(!rec_on()||rec_depth>0)return;
  rec_token(1,"i");
  rec_token(1,"g");
}

void rec_auto_key(int ch)
{
  char tok[16];
  if(!rec_on()||rec_depth>0)return;
  char_token(ch,tok);
  rec_token(2,tok);
  rec_context=2;
}

/* Grab mode in the AUTO window: the keys that move the cursor along the
   diagram are replaced by one "auto goto" line when a point is chosen.  The
   key that started Grab mode (g) is taken back. */
void rec_auto_traverse_begin(void)
{
  if(!rec_on())return;
  if(rec_kind==2){
    int n=(int)strlen(rec_buf);
    if(n>=2&&strcmp(rec_buf+n-2," g")==0){
      rec_buf[n-2]=0;
      if(strcmp(rec_buf,"autokeys")==0){rec_buf[0]=0;rec_kind=0;}
      rec_rewrite();
    }
  }
  rec_traverse++;
  rec_nokeys++;
}

void rec_auto_traverse_end(void)
{
  if(rec_traverse>0){rec_traverse--;if(rec_nokeys>0)rec_nokeys--;}
}

void rec_auto_goto(const char *sym,int label,int branch,int nth)
{
  char line[160];
  if(!rec_on())return;
  rec_flush();
  if(label>0)snprintf(line,sizeof(line),"auto goto lab %d",label);
  else if(sym!=NULL&&nth>0)snprintf(line,sizeof(line),"auto goto %s %d branch %d",sym,nth,branch);
  else {
    rec_write("# a point without a label or a type was grabbed here; it cannot be written as \"auto goto\"");
    return;
  }
  rec_write(line);
}

/* ---- prompts ---------------------------------------------------------- */
static int rec_push(const char *init)
{
  if(!rec_on()||rec_depth>=REC_MAXDEPTH)return 0;
  rec_flush();
  rec_init[rec_depth]=init?strdup(init):NULL;
  rec_depth++;
  return 1;
}

static char *rec_pop(void)
{
  char *s;
  if(rec_depth<=0)return NULL;
  rec_depth--;
  s=rec_init[rec_depth];
  rec_init[rec_depth]=NULL;
  return s;
}

/* a short comment that says what the following answer is for */
static void rec_label(const char *label)
{
  char c[200];
  int n,i;
  if(label==NULL)return;
  while(*label=='*'||isdigit((unsigned char)*label)){
    if(*label=='*'&&isdigit((unsigned char)label[1])){label+=2;break;}
    break;
  }
  snprintf(c,sizeof(c),"# %s",label);
  for(i=0;c[i];i++)if(c[i]=='\n'||c[i]=='\r')c[i]=' ';
  n=(int)strlen(c);
  while(n>2&&(isspace((unsigned char)c[n-1])||c[n-1]==':'))c[--n]=0;
  rec_write(c);
}

static void rec_ans_text(const char *text)
{
  char line[1200];
  int n=(int)strlen(text);
  int q=(n==0||isspace((unsigned char)text[0])||isspace((unsigned char)text[n-1])
         ||strchr(text,'#')!=NULL||text[0]=='?'||text[0]=='"'||text[0]=='\'');
  if(n>1000){rec_write("# the answer was too long to be recorded");return;}
  if(q){
    char qc=(strchr(text,'"')==NULL)?'"':'\'';
    if(qc=='\''&&strchr(text,'\'')!=NULL)q=0;      /* both kinds: leave it bare */
    else {snprintf(line,sizeof(line),"ans %c%s%c",qc,text,qc);rec_write(line);return;}
  }
  snprintf(line,sizeof(line),"ans %s",text);
  rec_write(line);
}

int rec_text_begin(const char *init){return rec_push(init);}

void rec_text_end(int h,int accepted,const char *label,const char *final)
{
  char *init;
  if(!h)return;
  init=rec_pop();
  rec_label(label);
  if(!accepted)rec_write("cancel");
  else if(init!=NULL&&strcmp(init,final)==0)rec_write("ok");
  else rec_ans_text(final);
  if(init)free(init);
}

/* a form keeps its first values in one block of n*stride bytes */
static char *rec_formvals[REC_MAXDEPTH];

int rec_form_begin(int n,char **names,char *vals,int stride)
{
  char *copy;
  (void)names;
  if(!rec_on()||rec_depth>=REC_MAXDEPTH)return 0;
  copy=(char *)malloc((size_t)n*(size_t)stride);
  if(copy==NULL)return 0;
  memcpy(copy,vals,(size_t)n*(size_t)stride);
  if(!rec_push(NULL)){free(copy);return 0;}
  rec_formvals[rec_depth-1]=copy;
  return 1;
}

void rec_form_end(int h,int accepted,const char *title,int n,char **names,char *vals,int stride)
{
  char line[1200],*first;
  int i,used=0,changed=0;
  if(!h)return;
  first=rec_formvals[rec_depth-1];
  rec_formvals[rec_depth-1]=NULL;
  { char *s=rec_pop(); if(s)free(s); }
  rec_label(title);
  if(!accepted)rec_write("cancel");
  else {
    used=snprintf(line,sizeof(line),"ans ");
    for(i=0;i<n;i++){
      const char *lb=names[i],*v=vals+(size_t)i*stride,*v0=first?first+(size_t)i*stride:"";
      char name[80];
      int k,m;
      if(strcmp(v,v0)==0)continue;
      if(strchr(v,'|')!=NULL||strchr(v,'#')!=NULL){
        rec_write("# a field contained | or # and could not be recorded");
        continue;
      }
      if(lb[0]=='*'&&isdigit((unsigned char)lb[1]))lb+=2;
      strncpy(name,lb,sizeof(name)-1);name[sizeof(name)-1]=0;
      m=(int)strlen(name);
      while(m>0&&(isspace((unsigned char)name[m-1])||name[m-1]==':'))name[--m]=0;
      for(k=0;name[k];k++)if(name[k]=='|'||name[k]=='=')name[k]='_';
      if(used+strlen(name)+strlen(v)+4>=sizeof(line))continue;
      used+=snprintf(line+used,sizeof(line)-used,"%s%s=%s",changed?"|":"",name,v);
      changed++;
    }
    if(changed)rec_write(line);
    else rec_write("ok");
  }
  if(first)free(first);
}

void rec_choice_end(int h,const char *prompt,int ch,const char *key)
{
  char *s;
  char line[40];
  if(!h)return;
  s=rec_pop();
  if(s)free(s);
  rec_label(prompt);
  ch=tolower((unsigned char)ch);
  if(ch!=0&&((key[0]&&ch==tolower((unsigned char)key[0]))||(key[0]&&key[1]&&ch==tolower((unsigned char)key[1])))){
    snprintf(line,sizeof(line),"ans %c",ch);
    rec_write(line);
  }
  else rec_write("cancel");
}

/* ---- mouse ------------------------------------------------------------ */
/* the same arithmetic as the click/drag lines use when a script is run */
static void pixel_to_data(int ix,int iy,double *x,double *y)
{
  double dx,dy;
  get_draw_area();
  dx=MyGraph->xhi-MyGraph->xlo;
  dy=MyGraph->yhi-MyGraph->ylo;
  *x=MyGraph->xlo+(double)(ix-DLeft)*dx/(double)(DRight-DLeft);
  *y=MyGraph->ylo+(double)(iy-DBottom)*dy/(double)(DTop-DBottom);
}

void rec_click(int ix,int iy)
{
  double x,y;
  char line[120];
  if(!rec_on()||rec_depth>0)return;
  rec_flush();
  if(DRight==DLeft||DTop==DBottom){snprintf(line,sizeof(line),"click pixel %d %d",ix,iy);rec_write(line);return;}
  pixel_to_data(ix,iy,&x,&y);
  snprintf(line,sizeof(line),"click %.10g %.10g",x,y);
  rec_write(line);
}

void rec_drag(int x1,int y1,int x2,int y2)
{
  double a,b,c,d;
  char line[200];
  if(!rec_on()||rec_depth>0)return;
  rec_flush();
  if(DRight==DLeft||DTop==DBottom){snprintf(line,sizeof(line),"drag pixel %d %d %d %d",x1,y1,x2,y2);rec_write(line);return;}
  pixel_to_data(x1,y1,&a,&b);
  pixel_to_data(x2,y2,&c,&d);
  snprintf(line,sizeof(line),"drag %.10g %.10g %.10g %.10g",a,b,c,d);
  rec_write(line);
}

/* a parameter or initial condition was changed by hand */
void rec_set(const char *name,double oldv,double newv)
{
  char line[160];
  if(!rec_on()||rec_depth>0||oldv==newv)return;
  rec_flush();
  snprintf(line,sizeof(line),"set %s=%.12g",name,newv);
  rec_write(line);
}

/* ---- start and stop ---------------------------------------------------- */
static void rec_snapshot(void)
{
  char line[300];
  int i,used;
  double v;
  rec_write("# state when the recording started");
  used=snprintf(line,sizeof(line),"set ");
  for(i=0;i<NUPAR;i++){
    char item[160];
    if(get_val(upar_names[i],&v)==0)continue;
    snprintf(item,sizeof(item),"%s%s=%.12g",used>4?";":"",upar_names[i],v);
    if(used+(int)strlen(item)>=72){rec_write(line);used=snprintf(line,sizeof(line),"set ");snprintf(item,sizeof(item),"%s=%.12g",upar_names[i],v);}
    used+=snprintf(line+used,sizeof(line)-used,"%s",item);
  }
  if(used>4)rec_write(line);
  used=snprintf(line,sizeof(line),"set ");
  for(i=0;i<NODE;i++){
    char item[160];
    snprintf(item,sizeof(item),"%s%s=%.12g",used>4?";":"",uvar_names[i],last_ic[i]);
    if(used+(int)strlen(item)>=72){rec_write(line);used=snprintf(line,sizeof(line),"set ");snprintf(item,sizeof(item),"%s=%.12g",uvar_names[i],last_ic[i]);}
    used+=snprintf(line+used,sizeof(line)-used,"%s",item);
  }
  if(used>4)rec_write(line);
}

static void rec_atexit(void)
{
  /* leave a complete file behind however XPP ends */
  if(rec_fp!=NULL){
    rec_flush();
    rec_write("# end of recording");
    fclose(rec_fp);
    rec_fp=NULL;
  }
}

static int rec_open(const char *name,int snapshot)
{
  static int registered=0;
  time_t t=time(NULL);
  char line[600],*nl;
  rec_fp=fopen(name,"w");
  if(rec_fp==NULL)return 0;
  if(!registered){atexit(rec_atexit);registered=1;}
  strncpy(rec_name,name,sizeof(rec_name)-1);rec_name[sizeof(rec_name)-1]=0;
  rec_kind=0;rec_buf[0]=0;rec_depth=0;rec_nokeys=0;rec_traverse=0;
  rec_ls=0;rec_lt=0;rec_lk=0;rec_lw=0;
  rec_echo=(getenv("XPP_REC_ECHO")!=NULL);
  rec_test=(getenv("XPP_REC_TEST")!=NULL);
  {char *ts=ctime(&t);nl=strchr(ts,'\n');if(nl)*nl=0;
   snprintf(line,sizeof(line),"# XPPAUT script recorded %s",ts);}
  rec_write(line);
  snprintf(line,sizeof(line),"# model: %s",this_file);
  rec_write(line);
  rec_write("# Edit freely.  Run it with File -> Do script, or:  xppaut model.ode -script this_file.xps");
  if(snapshot)rec_snapshot();
  return 1;
}

void script_record_stop(void)
{
  if(rec_fp==NULL)return;
  /* the keys that chose "Stop recording" are not part of the recording */
  if(rec_kind==1){
    int n=(int)strlen(rec_buf);
    if(n>=4&&strcmp(rec_buf+n-4," f o")==0){
      rec_buf[n-4]=0;
      if(strcmp(rec_buf,"keys")==0){rec_buf[0]=0;rec_kind=0;}
      rec_rewrite();
    }
  }
  rec_flush();
  rec_write("# end of recording");
  fclose(rec_fp);
  rec_fp=NULL;
  while(rec_depth>0){char *s=rec_pop();if(s)free(s);}
  if(Xup){
    char m[600];
    snprintf(m,sizeof(m),"Recording saved to %s",rec_name);
    bottom_msg(0,m);
  }
}

void script_record_from_menu(void)
{
  char filename[XPP_MAX_NAME];
  int snap;
  if(rec_fp!=NULL){script_record_stop();return;}
  if(active){
    err_msg("A script is running -- it cannot be recorded while it runs");
    return;
  }
  default_script_name(filename,sizeof(filename));
  if(!file_selector("Record script",filename,"*.xps"))return;
  {
    FILE *t=fopen(filename,"r");
    if(t!=NULL){
      fclose(t);
      if(TwoChoice("YES","NO","File exists. Overwrite?","yn")!='y')return;
    }
  }
  snap=(TwoChoice("YES","NO","Start with the current parameters?","yn")=='y');
  if(!rec_open(filename,snap)){
    err_msg("Cannot write the script file");
    return;
  }
  if(Xup){
    char m[600];
    snprintf(m,sizeof(m),"Recording to %s -- File -> recOrd to stop",filename);
    bottom_msg(0,m);
  }
}

/* start a recording from the environment (used for testing and demos) */
void script_record_from_env(void)
{
  const char *f=getenv("XPP_RECORD");
  if(f!=NULL&&*f&&rec_fp==NULL){
    if(!rec_open(f,0))fprintf(stderr,"XPP_RECORD: cannot write %s\n",f);
  }
}
