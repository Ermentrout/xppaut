/* This makes a big box with windows that have the names of the variables and
 * their current initial data, parameters, BCs etc
 *
 * This also has the slider boxes
 * This also has the file selector gadget
 * This also has the clone gadget
 */
#include "init_conds.h"

#ifndef HAVE_WCTYPE_H
# include <ctype.h>
#else
# include <wctype.h>
#endif
#include <dirent.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "arrayplot.h"
#include "calc.h"
#include "delay_handle.h"
#include "derived.h"
#include "eig_list.h"
#include "form_ode.h"
#include "ggets.h"
#include "graf_par.h"
#include "integrate.h"
#include "load_eqn.h"
#include "lunch-new.h"
#include "main.h"
#include "many_pops.h"
#include "markov.h"
#include "menudrive.h"
#include "mykeydef.h"
#include "nullcline.h"
#include "parserslow.h"
#include "pop_list.h"
#include "read_dir.h"
#include "util/strutil.h"
#include "tabular.h"
#include "util/timeutil.h"
#include "solver/volterra2.h"
#include "xpplim.h"

#include "bitmap/bc.bitmap"
#include "bitmap/delay.bitmap"
#include "bitmap/filebrowse.bitmap"
#include "bitmap/home.bitmap"
#include "bitmap/ic.bitmap"
#include "bitmap/lineup.bitmap"
#include "bitmap/linedn.bitmap"
#include "bitmap/pageup.bitmap"
#include "bitmap/pagedn.bitmap"
#include "bitmap/param.bitmap"
#include "bitmap/start.bitmap"


/* --- Macros --- */
#define FILESELNWIN 10

#define HOTWILD 2
#define HOTFILE 1

#define BOXEVENT   (ButtonPressMask 	|\
	KeyPressMask		|\
	ExposureMask		|\
	StructureNotifyMask	|\
	LeaveWindowMask		|\
	EnterWindowMask)

#define EDIT_WAIT 0
#define EDIT_NEXT 1
#define EDIT_ESC 2
#define EDIT_DONE 3

#define SB_DIM 5
#define SB_SPC 2


/* --- Types --- */
typedef struct {
	int n,n0,here;
	Window base,cancel,ok,up,dn,pgup,pgdn,file,wild,w[FILESELNWIN],dir,home,start;
	Window fw,ww;
	char wildtxt[DEFAULT_STRING_LENGTH],filetxt[DEFAULT_STRING_LENGTH];
	int nwin,minwid,minhgt;
	int off,pos,hot;
	char title[DEFAULT_STRING_LENGTH];
} FILESEL;

typedef struct {
	int use,pos,l;
	char parname[20];
	double lo,hi,val;
	int hgt;
	int type,index;
	Window left,right,top,main,slide,go;
} PAR_SLIDER;


/* --- Forward declarations --- */
static void add_edit_float(BoxList *b, int i, double z);
static void add_editval(BoxList *b, int i, char *string);
static void box_enter(BoxList b, Window w, int val);
static void box_list_scroll(BoxList *b, int i);
static int button_selector(Window w);
static void check_box_cursor(void);
static void create_file_selector(char *title, char *file, char *wild);
static void crossing_selector(Window w, int c);
static void destroy_box(BoxList *b);
static void destroy_selector(void);
static void display_box(BoxList b, Window w);
static void display_file_sel(FILESEL f, Window w);
static void do_box_button(BoxList *b, Window w);
static void do_box_key(BoxList *b, XEvent ev, int *used);
static int do_file_select_events(void);
static void do_slide_button(int w, PAR_SLIDER *p);
static void do_slide_motion(Window w, int x, PAR_SLIDER *p,int state);
static void do_slide_release(int w, PAR_SLIDER *p);
static void draw_editable(Window win, char *string, int off, int cursor, int mc);
static void draw_slider(Window w, int x, int hgt, int l);
static int edit_bitem(BoxList *b, int i, int ch);
static int edit_fitem(int ch, char *string, Window w, int *off1, int *pos1, int mc);
static void enter_slider(Window w, PAR_SLIDER *p, int val);
static void expose_selector(Window w);
static void expose_slider(Window w, PAR_SLIDER *p);
static void fs_scroll(int i);
static void get_nrow_from_hgt(int h, int *n, int *w);
static void justify_string(Window w1, char *s1);
static void load_entire_box(BoxList *b);
static void make_box_list(BoxList *b, char *wname, char *iname, int n, int type, int use);
static void make_box_list_window(BoxList *b, int type);
static void make_par_slider(Window base, int x, int y, int width, int index, int var_index);
static void new_wild(void);
static void put_edit_cursor(Window w, int pos);
static void redraw_directory(void);
static void redraw_entire_box(BoxList *b);
static void redraw_file_list(void);
static void redraw_fs_text(char *string, Window w, int flag);
static void redraw_slide(PAR_SLIDER *p);
static int selector_key(XEvent ev);
static void set_default_ics(void);
static void set_default_params(void);
static void set_slide_pos(PAR_SLIDER *p);
static void set_up_xvt(void);
static void set_up_pp(void);
static void set_up_arry(void);
static void set_value_from_box(BoxList *b, int i);
static int to_float(char *s, double *z);


/* --- Data --- */
static int HotBoxItem=-1;
static FILESEL filesel;
static PAR_SLIDER my_par_slide[3];
static BoxList *HotBox;
static BoxList ICBox;
static BoxList DelayBox;
static BoxList BCBox;

char SLIDERVAR[3][20];
double SLIDERLO[3] = {0, 0, 0};
double SLIDERHI[3] = {1, 1, 1};

BoxList ParamBox;


/* --- Functions --- */
void box_buttons(Window w) {

	if(ICBox.xuse) {
		do_box_button(&ICBox,w);
	}
	if(BCBox.xuse) {
		do_box_button(&BCBox,w);
	}
	if(DelayBox.xuse) {
		do_box_button(&DelayBox,w);
	}
	if(ParamBox.xuse) {
		do_box_button(&ParamBox,w);
	}
}


void box_enter_events(Window w, int yn) {
	int i;
	int val;
	if(yn==1) {
		val=2;
	} else {
		val=1;
	}
	if(ICBox.xuse) {
		box_enter(ICBox,w,val);
	}
	if(BCBox.xuse) {
		box_enter(BCBox,w,val);
	}
	if(ParamBox.xuse) {
		box_enter(ParamBox,w,val);
	}
	if(DelayBox.xuse) {
		box_enter(DelayBox,w,val);
	}
	if(ICBox.xuse && (w==ICBox.xvt || w==ICBox.pp || w==ICBox.arr)) {
		XSetWindowBorderWidth(display,w,val);
	}
	if(ICBox.xuse==0) {
		return;
	}
	for(i=0;i<ICBox.nwin;i++) {
		if(w==ICBox.ck[i]) {
			XSetWindowBorderWidth(display,w,val);
		}
	}
}


void box_keypress(XEvent ev, int *used) {
	if(ICBox.xuse) {
		do_box_key(&ICBox,ev,used);
		if(*used) {
			return;
		}
	}
	if(BCBox.xuse) {
		do_box_key(&BCBox,ev,used);
		if(*used) {
			return;
		}
	}
	if(DelayBox.xuse) {
		do_box_key(&DelayBox,ev,used);
		if(*used) {
			return;
		}
	}
	if(ParamBox.xuse) {
		do_box_key(&ParamBox,ev,used);
		if(*used) {
			return;
		}
	}
}


/* CLONE */
void clone_ode(void) {
	int i,j,x,y;
	FILE *fp;
	char clone[XPP_MAX_NAME];

	char *s;
	time_t ttt;
	double z;
	clone[0]=0;
	if(!file_selector("Clone ODE file",clone,"*.ode")) {
		return;
	}
	if((fp=fopen(clone,"w"))==NULL) {
		err_msg(" Cant open clone file");
		return;
	}
	ttt=time(0);
	fprintf(fp,"# clone of %s on %s",this_file,ctime(&ttt));
	for(i=0;i<NLINES;i++) {
		s=save_eqn[i];
		if(s[0]=='p' || s[0]=='P' || s[0]=='b' || s[0]=='B') {
			x=find_char(s,"'",0,&j);
			y=find_char(s,"=",0,&j);
			if(x!=0 || y!=0) {
				fprintf(fp,"# original\n# %s\n",s);
				continue;
			}
		}
		if(strncasecmp("done",s,4)==0) {
			continue;
		}
		fprintf(fp,"%s\n",s);
	}
	fprintf(fp,"# Cloned parameters etc here\n");
	/* now we do parameters boundary conds and ICs */
	j=0;
	fprintf(fp,"init ");
	for(i=0;i<(NODE+NMarkov);i++) {
		if(j==8) {
			fprintf(fp,"\ninit ");
			j=0;
		}
		fprintf(fp,"%s=%g ",uvar_names[i],last_ic[i]);
		j++;
	}
	fprintf(fp,"\n");

	/* BDRY conds */
	if(my_bc[0].string[0]!='0') {
		for(i=0;i<NODE;i++) {
			fprintf(fp,"bdry %s\n",my_bc[i].string);
		}
	}
	j=0;
	if(NUPAR>0) {
		fprintf(fp,"par ");
		for(i=0;i<NUPAR;i++) {
			if(j==8) {
				fprintf(fp,"\npar ");
				j=0;
			}
			get_val(upar_names[i],&z);
			fprintf(fp,"%s=%g ",upar_names[i],z);
			j++;
		}
	}
	fprintf(fp,"\n");
	fprintf(fp,"done \n");
	fclose(fp);
}


void create_par_sliders(Window base, int x0, int h0, int var_ind[3]) {
	int i;
	for(i=0;i<3;i++) {
		make_par_slider(base,x0+i*36*DCURXs,h0,100,i, var_ind[i]);
	}
}


/* new code is a bit tricky here - we dont want
   to draw it if it is not visible
	there are nwin windows covering indexes
	n0,n0+1,...n0+nwin-1
	if the index is beyond this dont draw it
*/
void draw_one_box(BoxList b, int index) {
	Window w,we;
	int n0=b.n0;
	int n1=n0+b.nwin-1;
	int i;
	if(b.xuse==0) {
		return;
	}
	if(index<n0 || index>n1) {
		return; /* don't draw the ones out of range*/
	}
	i=index-n0;
	w=b.w[i];
	we=b.we[i];
	switch(b.type) {
	case PARAMETER_BOX:
		draw_editable(we,b.value[index],
					  b.off[index],b.pos[index],b.mc);
		justify_string(w,upar_names[index]);
		break;
	case BOUNDARY_CONDITION_BOX:
		justify_string(w,my_bc[index].name);
		draw_editable(we,b.value[index],b.off[index],
					  b.pos[index],b.mc);
		break;
	case INITIAL_CONDITION_BOX:
		draw_editable(we,b.value[index],b.off[index],
					  b.pos[index],b.mc);
		justify_string(w,uvar_names[index]);
		break;
	case DELAY_BOX:
		justify_string(w,uvar_names[index]);
		draw_editable(we,b.value[index],b.off[index],
					  b.pos[index],b.mc);
		break;
	}
}


void enter_slides(Window w, int val) {
	int i;
	for(i=0;i<3;i++) {
		enter_slider(w,&my_par_slide[i],val);
	}
}


/* this is added to take care of making sure
	 exposure of the boxes is easily taken care of
  */
void expose_box(Window w) {
	if(ICBox.xuse) {
		display_box(ICBox,w);
	}
	if(BCBox.xuse) {
		display_box(BCBox,w);
	}
	if(ParamBox.xuse) {
		display_box(ParamBox,w);
	}
	if(DelayBox.xuse) {
		display_box(DelayBox,w);
	}
}


void expose_slides(Window w) {
	int i;
	for(i=0;i<3;i++) {
		expose_slider(w,&my_par_slide[i]);
	}
}


int file_selector(char *title, char *file, char *wild) {
	int i;
	if(!get_directory(cur_dir)) {
		return 0;
	}
	if(!get_fileinfo(wild,cur_dir,&my_ff)) {
		return 0;
	}
	create_file_selector(title,file,wild);
	i=do_file_select_events();
	destroy_selector();
	XFlush(display);/*Need to do this otherwise the file dialog hangs around*/
	if(i==0) {
		return 0;
	}
	strcpy(file,filesel.filetxt);
	return 1; /* got a file name */
}


int find_user_name(int type, char *oname) {
	char name[25];
	int j=0,k=0,i=-1;
	for(j=0;j<strlen(oname);j++) {
		if(!isspace(oname[j])) {
			name[k]=oname[j];
			k++;
		}
	}
	name[k]=0;
	for(i=0;i<NUPAR;i++) {
		if((type==PARAMETER_BOX) && (strcasecmp(upar_names[i],name)==0)) {
			break;
		}
	}
	if(i<NUPAR) {
		return(i);
	}
	for(i=0;i<NEQ;i++) {
		if((type==INITIAL_CONDITION_BOX) && (strcasecmp(uvar_names[i],name)==0)) {
			break;
		}
	}
	if(i<NEQ) {
		return(i);
	}
	return(-1);
}


void make_new_bc_box(void) {
	if(BCBox.xuse) {
		XRaiseWindow(display,BCBox.base);
		return;
	}
	make_box_list_window(&BCBox,BOUNDARY_CONDITION_BOX);
	make_icon((char*)bc_bits,bc_width,bc_height,BCBox.base);
}


void initialize_box(void) {
	make_box_list(&ICBox,"Initial Data","ICs",NODE+NMarkov,INITIAL_CONDITION_BOX,1);
	if(NUPAR>0) {
		make_box_list(&ParamBox,"Parameters","Par",NUPAR,PARAMETER_BOX,1);
	} else {
		ParamBox.use=0;
	}
	if(NDELAYS>0) {
		make_box_list(&DelayBox,"Delay ICs","Delay", NODE,DELAY_BOX,1);
	} else {
		DelayBox.use=0;
	}
	make_box_list(&BCBox,"Boundary Conds","BCs",NODE,BOUNDARY_CONDITION_BOX,1);
}


void make_new_delay_box(void) {
	if(DelayBox.use==0) {
		return;
	}
	if( DelayBox.xuse==1) {
		XRaiseWindow(display,DelayBox.base);
		return;
	}
	make_box_list_window(&DelayBox,DELAY_BOX);
	make_icon((char*)delay_bits,delay_width,delay_height,DelayBox.base);
}


void make_new_ic_box(void) {
	if(ICBox.xuse) {
		XRaiseWindow(display,ICBox.base);
		return;
	}
	make_box_list_window(&ICBox,INITIAL_CONDITION_BOX);
	make_icon((char*)ic_bits,ic_width,ic_height,ICBox.base);
}


void make_new_param_box(void) {
	if(ParamBox.use==0) {
		return;
	}
	if(ParamBox.xuse==1) {
		XRaiseWindow(display,ParamBox.base);
		return;
	}
	make_box_list_window(&ParamBox,PARAMETER_BOX);
	make_icon((char*)param_bits,param_width,param_height,ParamBox.base);
}


void man_ic(void) {
	int done,index=0;
	double z;
	char name[XPP_MAX_NAME],junk[DEFAULT_STRING_LENGTH];
	while(1) {
		sprintf(name,"%s :",uvar_names[index]);
		z=last_ic[index];
		done=new_float(name,&z);
		if(done==0) {
			last_ic[index]=z;
			sprintf(junk,"%.16g",z);
			set_edit_params(&ICBox,index,junk);
			draw_one_box(ICBox,index);
			index++;
			if(index>=NODE+NMarkov) {
				return;
			}
		}
		if(done==-1) {
			return;
		}
	}
}


void new_parameter(void) {
	int done,index;
	double z;
	char name[DEFAULT_STRING_LENGTH],value[DEFAULT_STRING_LENGTH],junk[DEFAULT_STRING_LENGTH];
	while(1) {
		name[0]=0;
		done=new_string("PARAMETER_BOX:",name);
		if(strlen(name)==0 || done==0) {
			redo_stuff();
			return;
		}
		if(strncasecmp(name,"DEFAULT",7  )==0) {
			set_default_params();
			continue;
		}

		if(strncasecmp(name,"!LOAD", 5 )==0) {
			io_parameter_file(name,READEM);
			continue;
		}
		if(strncasecmp(name,"!SAVE", 5 )==0) {
			io_parameter_file(name,WRITEM);
			continue;
		} else {
			index=find_user_name(PARAMETER_BOX,name);
			if(index>=0) {
				get_val(upar_names[index],&z);
				sprintf(value,"%s :",name);
				done=new_float(value,&z);
				if(done==0) {
					set_val(upar_names[index],z);
					sprintf(junk,"%.16g",z);
					set_edit_params(&ParamBox,index,junk);
					draw_one_box(ParamBox,index);
					reset_sliders();
				}
				if(done==-1) {
					redo_stuff();
					return;
				}
			}
		}
	}
}


void redo_stuff(void) {
	evaluate_derived();
	re_evaluate_kernels();
	redo_all_fun_tables();
	evaluate_derived();
}


void redraw_bcs(void) {
	int i;
	for(i=0;i<NODE;i++) {
		draw_one_box(BCBox,i);
	}
}


void redraw_delays(void) {
	int i;
	if(DelayBox.use) {
		for(i=0;i<NODE;i++) {
			draw_one_box(DelayBox,i);
		}
	}
}


void redraw_ics(void) {
	int i,in;
	for(i=0;i<NODE+NMarkov;i++) {
		add_edit_float(&ICBox,i,last_ic[i]);
		draw_one_box(ICBox,i);
	}
	reset_sliders();
	if(ICBox.xuse==0) {
		return;
	}
	for(i=0;i<ICBox.nwin;i++) {
		in=i+ICBox.n0;
		if(ICBox.isck[in]) {
			XDrawString(display,ICBox.ck[i],small_gc,0,CURY_OFFs,"*",1);
		} else {
			XClearWindow(display,ICBox.ck[i]);
		}
	}
}


void redraw_params(void) {
	int i;
	double z;
	evaluate_derived();
	if(ParamBox.use) {
		for(i=0;i<NUPAR;i++) {
			get_val(upar_names[i],&z);
			add_edit_float(&ParamBox,i,z);
			draw_one_box(ParamBox,i);
		}
	}
	reset_sliders();
}


void reset_sliders(void) {
	int i;
	double val;
	PAR_SLIDER *p;
	for(i=0;i<3;i++) {
		p=&my_par_slide[i];
		if(p->use) {
			if(p->type==INITIAL_CONDITION_BOX) {
				val=last_ic[p->index];
			} else {
				get_val(p->parname,&val);
			}
			p->val=val;
			set_slide_pos(p);
			expose_slider(p->slide,p);
			expose_slider(p->top,p);
		}
	}
}


void resize_par_box(Window win) {
	unsigned int h,w;
	int nwin = 0;
	int ok=0;
	BoxList *b;
	if(ICBox.xuse==1 && win==ICBox.base) {
		ok=1;
		b=&ICBox;
		get_new_size(win,&w,&h);
		get_nrow_from_hgt(h,&nwin,(int*)&w);
	}
	if(ParamBox.xuse==1 && win==ParamBox.base) {
		ok=2;
		b=&ParamBox;
		waitasec(ClickTime);
		get_new_size(win,&w,&h);
		get_nrow_from_hgt(h,&nwin,(int*)&w);
	}
	if(BCBox.xuse==1 && win==BCBox.base) {
		ok=3;
		b=&BCBox;
		get_new_size(win,&w,&h);
		get_nrow_from_hgt(h,&nwin,(int*)&w);
	}
	if(DelayBox.xuse==1 && win==DelayBox.base) {
		ok=4;
		b=&DelayBox;
		get_new_size(win,&w,&h);
		get_nrow_from_hgt(h,&nwin,(int*)&w);
	}
	if(ok==0) {
		return;
	}
	if(nwin>b->n) {
		nwin=b->n;
	}
	if(nwin==b->nwin) {
		return;
	}
	b->nwin=nwin;
	b->n0=(b->n-b->nwin);
	/* b->n0=0;This is a work-around for the following bug (still happening):
			 i) User scrolls to bottom of parameter list
		 ii) then resizes window
		 iii) then Segmentation Fault
*/
	switch(ok) {
	case INITIAL_CONDITION_BOX:
		destroy_box(&ICBox);
		make_new_ic_box();
		break;
	case PARAMETER_BOX:
		destroy_box(&ParamBox);
		make_new_param_box();
		break;
	case DELAY_BOX:
		destroy_box(&DelayBox);
		make_new_delay_box();
		break;
	case BOUNDARY_CONDITION_BOX:
		destroy_box(&BCBox);
		make_new_bc_box();
		break;
	}
}


void resize_par_slides(int h) {
	int i;
	for(i=0;i<3;i++) {
		XMoveResizeWindow(display,my_par_slide[i].main,10+36*i*DCURXs,
						  h,32*DCURXs,3*(DCURYs+2));
	}
}


void set_edit_params(BoxList *b, int i, char *string) {
	int l=strlen(string);
	strcpy(b->value[i],string);
	b->off[i]=0;
	if(l>b->mc) {
		b->pos[i]=b->mc;
	} else {
		b->pos[i]=l;
	}
}


void slide_button_press(Window w) {
	int i;
	for(i=0;i<3;i++) {
		do_slide_button(w,&my_par_slide[i]);
	}
}


void slide_release(Window w) {
	int i;
	for(i=0;i<3;i++) {
		do_slide_release(w,&my_par_slide[i]);
	}
}


void slider_motion(XEvent ev) {
	int x,i;
	Window w;
	w=ev.xmotion.window;
	x=ev.xmotion.x;
	/* printf(" state=%d\n",ev.xmotion.state); */
	for(i=0;i<3;i++) {
		do_slide_motion(w,x,&my_par_slide[i],ev.xmotion.state);
	}
}


/* --- static functions --- */
static void add_edit_float(BoxList *b, int i, double z) {
	char junk[DEFAULT_STRING_LENGTH];
	sprintf(junk,"%.16g",z);
	add_editval(b,i,junk);
}


static void add_editval(BoxList *b, int i, char *string) {
	int n0=b->n0,n1=b->n0+b->nwin-1;
	int iw;
	set_edit_params(b,i,string);
	if(i<n0 || i>n1) {
		return;
	}
	iw=i-n0;
	if(b->xuse) {
		draw_editable(b->we[iw],string,b->off[i],b->pos[i],b->mc);
	}
}


static void box_enter(BoxList b, Window w, int val) {
	if(w==b.ok || w==b.cancel || w==b.def || w==b.go || w==b.close ||
	   w==b.dn || w==b.up || w==b.pgdn || w==b.pgup) {
		XSetWindowBorderWidth(display,w,val);
	}
}


static void box_list_scroll(BoxList *b, int i) {
	int n0=b->n0;
	int new_int;
	int nw=b->nwin,n=b->n;
	int nend;
	if(n<=nw) {
		return;  /* do nothing - there is nothing to do */
	}
	new_int =n0-i;
	nend=new_int+nw;
	if(new_int<0) {
		new_int=0;
	}
	if(nend>n) {
		new_int=n-nw;
	}
	b->n0=new_int;
	switch(b->type) {
	case PARAMETER_BOX:
		load_entire_box(b);
		redraw_params();
		reset_sliders();
		break;
	case BOUNDARY_CONDITION_BOX:
		load_entire_box(b);
		redraw_bcs();
		break;
	case INITIAL_CONDITION_BOX:
		load_entire_box(b);
		redraw_ics();
		reset_sliders();
		break;
	case DELAY_BOX:
		load_entire_box(b);
		redraw_delays();
		break;
	}
}


static int button_selector(Window w) {
	int i,i0;
	int k,n=filesel.n;
	if(w==filesel.ok) {
		return 1;
	}
	if(w==filesel.cancel) {
		return 2;
	}
	if(w==filesel.up) {
		fs_scroll(1);
	}
	if(w==filesel.dn) {
		fs_scroll(-1);
	}
	if(w==filesel.pgup) {
		fs_scroll(filesel.nwin);
	}
	if(w==filesel.pgdn) {
		fs_scroll(-filesel.nwin);
	}
	if(w==filesel.home) {
		char* HOMEDIR = getenv("HOME");
		if((HOMEDIR == NULL)  ||  (strlen(HOMEDIR)==0)) {
			plintf("User's HOME environment variable not set.\n");
			return 0;
		}
		change_directory(HOMEDIR);

		get_directory(cur_dir);
		redraw_directory();
		free_finfo(&my_ff); /* delete the old file info */
		filesel.n0=0; /* back to the top of the list */
		get_fileinfo(filesel.wildtxt,cur_dir,&my_ff);
		filesel.n=my_ff.ndirs+my_ff.nfiles;
		strcpy(filesel.filetxt,cur_dir);

		int m = strlen(filesel.filetxt);
		if(filesel.filetxt[m-1] != '/') {
			strcat(filesel.filetxt,"/");
		}
		redraw_fs_text(filesel.filetxt,filesel.file,1);
		redraw_file_list();
		XFlush(display);
		return 0;
	}
	if(w==filesel.start) {
		char* START = getenv("XPPSTART");
		if((START == NULL)  ||  (strlen(START)==0)) {
			plintf("User's XPPSTART environment variable not set.\n");
			return 0;
		}

		change_directory(START);
		get_directory(cur_dir);
		redraw_directory();
		free_finfo(&my_ff); /* delete the old file info */
		filesel.n0=0; /* back to the top of the list */
		get_fileinfo(filesel.wildtxt,cur_dir,&my_ff);
		filesel.n=my_ff.ndirs+my_ff.nfiles;
		strcpy(filesel.filetxt,cur_dir);

		int m = strlen(filesel.filetxt);
		if(filesel.filetxt[m-1] != '/') {
			strcat(filesel.filetxt,"/");
		}

		redraw_fs_text(filesel.filetxt,filesel.file,1);
		redraw_file_list();
		XFlush(display);

		return 0;
	}
	if(w==filesel.file) { /* selected the file text */
		if(filesel.hot!=HOTFILE) {
			filesel.pos=strlen(filesel.filetxt);
		}
		filesel.hot=HOTFILE;
		redraw_fs_text(filesel.filetxt,filesel.file,1);
		redraw_fs_text(filesel.wildtxt,filesel.wild,0);
		/* set up text stuff */
		return 0;
	}
	if(w==filesel.wild) {
		if(filesel.hot!=HOTWILD) {
			filesel.pos=strlen(filesel.wildtxt);
		}
		filesel.hot=HOTWILD;
		redraw_fs_text(filesel.filetxt,filesel.file,0);
		redraw_fs_text(filesel.wildtxt,filesel.wild,1);
		return 0;
	}
	i0=-1;
	for(i=0;i<filesel.nwin;i++) {
		if(w==filesel.w[i]) {
			i0=i;
		}
	}
	if(i0>-1) { /* clicked on a file or directory */
		k=i0+filesel.n0;
		if(k<my_ff.ndirs) { /* it is a directory so we should reset */
			change_directory(my_ff.dirnames[k]);
			get_directory(cur_dir);
			redraw_directory();
			free_finfo(&my_ff); /* delete the old file info */
			filesel.n0=0; /* back to the top of the list */
			get_fileinfo(filesel.wildtxt,cur_dir,&my_ff);
			filesel.n=my_ff.ndirs+my_ff.nfiles;
			strcpy(filesel.filetxt,cur_dir);

			int m = strlen(filesel.filetxt);
			if(filesel.filetxt[m-1] != '/') {
				strcat(filesel.filetxt,"/");
			}

			redraw_fs_text(filesel.filetxt,filesel.file,1);
			redraw_file_list();
			XFlush(display);
			return 0;
		}
		if(k<n) {
			k=k-my_ff.ndirs;
			strcpy(filesel.filetxt,my_ff.filenames[k]);
			return 1; /* got a file */
		}
	}
	return 0;
}


static void check_box_cursor(void) {
	int n0;
	if(HotBoxItem<0 || HotBox->xuse==0) {
		return;
	}
	n0=HotBox->n0;
	draw_editable(HotBox->we[HotBoxItem],HotBox->value[HotBoxItem+n0],
			HotBox->off[HotBoxItem],HotBox->pos[HotBoxItem],
			HotBox->mc);
	HotBoxItem=-1;
}


static void create_file_selector(char *title, char *file, char *wild) {
	int n=my_ff.ndirs+my_ff.nfiles;
	int nwin=FILESELNWIN;
	/*int wid,hgt,i;
 */
	int hgt,i;
	int width,height;

	Window base;
	XTextProperty winname;
	XSizeHints size_hints;
	filesel.n=n;
	filesel.n0=0;
	filesel.nwin=nwin;
	strcpy(filesel.title,title);
	strcpy(filesel.wildtxt,wild);
	strcpy(filesel.filetxt,file);
	width=80*DCURXs;
	/*wid=30*DCURXs;*/
	hgt=DCURYs+4;
	height=(5+nwin)*hgt;
	filesel.minwid=width;
	filesel.minhgt=height;
	/* plintf("Title=%s\n",title); */
	/* printf("Here now 23!\n");*/
	base=make_plain_window(RootWindow(display,screen),0,0,width,height,4);
	/* printf("Here now 23!\n"); */
	filesel.base=base;
	XStringListToTextProperty(&title,1,&winname);
	size_hints.flags=PPosition|PSize|PMinSize;
	size_hints.x=0;
	size_hints.y=0;
	size_hints.width=width;
	size_hints.height=height;
	size_hints.min_width=width;
	size_hints.min_height=height;
	size_hints.max_width=width;
	size_hints.max_height=height;

	make_icon((char*)filebrowse_bits,filebrowse_width,filebrowse_height,base);

	XClassHint class_hints;
	class_hints.res_name="";
	class_hints.res_class="";

	XSetWMProperties(display,base,&winname,NULL,NULL,0,&size_hints,NULL,&class_hints);

	filesel.up=make_icon_window(base,DCURXs,2+3*hgt+72+15,32,24,1,0,0,lineup_bits);
	filesel.dn=make_icon_window(base,DCURXs,2+3*hgt+96+18,32,24,1,0,0,linedn_bits);
	filesel.pgup=make_icon_window(base,DCURXs,2+3*hgt+48+12,32,24,1,0,0,pageup_bits);
	filesel.pgdn=make_icon_window(base,DCURXs,2+3*hgt+120+21,32,24,1,0,0,pagedn_bits);
	filesel.home=make_icon_window(base,DCURXs,2+3*hgt,32,24,1,0,0,home_bits);
	filesel.start=make_icon_window(base,DCURXs,2+3*hgt+24+3,32,24,1,0,0,start_bits);

	filesel.dir=make_plain_window(base,7*DCURXs,2,width-7*DCURXs-5,DCURYs,0);
	filesel.wild=make_plain_window(base,7*DCURXs,2+hgt,width-7*DCURXs-5,DCURYs,1);
	filesel.ww=make_window(base,2,2+hgt,6*DCURXs+2,DCURYs,0);
	filesel.file=make_plain_window(base,7*DCURXs,2+2*hgt,width-7*DCURXs-5,DCURYs,1);
	filesel.fw=make_window(base,2,2+2*hgt,6*DCURXs+2,DCURYs,0);
	for(i=0;i<nwin;i++) {
		filesel.w[i]=make_plain_window(base,6*DCURXs+5,2+(3+i)*hgt,width-6*DCURXs-10,DCURYs,0);
	}
	filesel.ok=make_window(base,width/2-7*DCURXs-3,height-hgt,7*DCURXs,DCURYs,1);
	filesel.cancel=make_window(base,width/2+3,height-hgt,7*DCURXs,DCURYs,1);
	filesel.here=1;
	filesel.hot=HOTFILE;
	filesel.pos=strlen(filesel.filetxt);
	filesel.off=0;
}


static void crossing_selector(Window w, int c) {
	int t1=1,t2=2,i;
	if(c==1) {
		t1=0;
		t2=1;
	}
	for(i=0;i<filesel.nwin;i++) {
		if(w==filesel.w[i]) {
			XSetWindowBorderWidth(display,w,t1);
			return;
		}
	}
	if(w==filesel.ok || w==filesel.cancel || w==filesel.pgup || w==filesel.pgdn ||
	   w==filesel.up || w==filesel.dn || w==filesel.file || w==filesel.wild ||
	   w==filesel.home || w==filesel.start) {
		XSetWindowBorderWidth(display,w,t2);
	}
}


static void destroy_box(BoxList *b) {
	if(b->xuse==0) {
		return;
	}
	b->xuse=0;
	XFlush(display);
	XSetInputFocus(display,main_win,RevertToParent,CurrentTime);
	if(b->use==0) {
		return;
	}
	waitasec(ClickTime);

	XDestroySubwindows(display,b->base);
	XDestroyWindow(display,b->base);

	free(b->w);
	free(b->we);
	if(b->type==INITIAL_CONDITION_BOX) {
		free(b->ck);
		free(b->isck);
	}
	waitasec(200);
	XFlush(display);
}


static void destroy_selector(void) {
	filesel.here=0;
	waitasec(ClickTime);
	XDestroySubwindows(display,filesel.base);
	XDestroyWindow(display,filesel.base);
	free_finfo(&my_ff);
}


static void display_box(BoxList b, Window w) {
	int i;
	int n0=b.n0;
	int n1=n0+b.nwin;
	int index;
	if(b.xuse==0) {
		return;
	}
	if(b.close==w) {
		XDrawString(display,w,small_gc,5,CURY_OFFs, "Close",5);
	}
	if(b.go==w) {
		XDrawString(display,w,small_gc,5,CURY_OFFs, "Go",2);
	}
	if(b.ok==w) {
		XDrawString(display,w,small_gc,5,CURY_OFFs, "Ok",2);
	}
	if(b.cancel==w) {
		XDrawString(display,w,small_gc,5,CURY_OFFs, "Cancel",6);
	}
	if(b.def==w) {
		XDrawString(display,w,small_gc,5,CURY_OFFs, "Default",7);
	}
	if(b.up==w) {
		/*XDrawString(display,w,small_gc,5+DCURX,CURY_OFFs,
	"^",1);
	*/
	}
	if(b.dn==w) {
		/*XDrawString(display,w,small_gc,5+DCURX,CURY_OFFs,
	"v",1);*/
	}
	if(b.pgup==w) {
		/*XDrawString(display,w,small_gc,5,CURY_OFFs,
		"^^",2);*/
	}
	if(b.pgdn==w) {
		/*XDrawString(display,w,small_gc,5,CURY_OFFs,
	"vv",2);*/
	}
	if(b.type==INITIAL_CONDITION_BOX) {
		if(b.xvt==w) {
			XDrawString(display,w,small_gc,3,CURY_OFFs,"xvst",4);
		}
		if(b.pp==w) {
			XDrawString(display,w,small_gc,3,CURY_OFFs,"xvsy",4);
		}
		if(b.arr==w) {
			XDrawString(display,w,small_gc,3,CURY_OFFs,"arry",4);
		}
	}
	for(i=0;i<b.nwin;i++) {
		if(b.w[i]==w || b.we[i]==w) {
			draw_one_box(b,i+b.n0);
			return;
		}
	}
	if(b.type==INITIAL_CONDITION_BOX) {
		for(i=0;i<b.nwin;i++) {
			index=i+b.n0;
			if(index>=n0 && index<n1) {
				if(b.ck[i]==w && b.isck[index]==1) {
					XDrawString(display,w,small_gc,5,CURY_OFFs,"*",1);
				}
			}
		}
	}
}


static void display_file_sel(FILESEL f, Window w) {
	int i,i0;
	Window root;
	int xloc;
	int yloc;

	unsigned int cwid;
	unsigned int chgt;
	unsigned int cbwid;
	unsigned int cdepth;

	XGetGeometry(display,f.base,&root,&xloc,&yloc,&cwid,&chgt,&cbwid,&cdepth);
	XResizeWindow(display,f.wild,cwid-7*DCURXs-5,DCURYs);
	XResizeWindow(display,f.file,cwid-7*DCURXs-5,DCURYs);
	for(i=0;i<f.nwin;i++) {
		XResizeWindow(display,f.w[i],cwid-6*DCURXs-10,DCURYs);
	}
	int hgt=DCURYs+4;
	XMoveResizeWindow(display,f.ok,cwid/2-7*DCURXs-3,chgt-hgt,7*DCURXs,DCURYs);
	XMoveResizeWindow(display,f.cancel,cwid/2+3,chgt-hgt,7*DCURXs,DCURYs);

	char t[DEFAULT_STRING_LENGTH];
	if(f.here!=1) {
		return;
	}

	if(f.ok==w) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Ok",2);
	}
	if(f.cancel==w) {
		XDrawString(display,w,small_gc,5,CURY_OFFs,"Cancel",6);
	}
	if(f.up==w) {
		/*XDrawString(display,w,small_gc,5+DCURX/2,CURY_OFFs,"^",1);
	*/
	}
	if(f.dn==w) {
		/*XDrawString(display,w,small_gc,5+DCURX/2,CURY_OFFs,"vv",1);
	*/
	}
	if(f.pgup==w) {
		/*XDrawString(display,w,small_gc,5,CURY_OFFs,"^^",2);
	*/
	}
	if(f.pgdn==w) {
		/* XDrawString(display,w,small_gc,5,CURY_OFFs,"vv",2);
	*/
	}
	if(f.file==w) {
		XClearWindow(display,w);
		XDrawString(display,w,small_gc,2,CURY_OFFs,
					f.filetxt,strlen(f.filetxt));
	}
	if(f.wild==w) {
		XClearWindow(display,w);
		XDrawString(display,w,small_gc,2, CURY_OFFs,f.wildtxt,strlen(f.wildtxt));
	}
	if(f.fw==w) {
		XDrawString(display,w,small_gc,5, CURY_OFFs,"File: ",6);
	}
	if(f.ww==w) {
		XDrawString(display,w,small_gc,5, CURY_OFFs,"Wild: ",6);
	}
	if(f.dir==w) {
		sprintf(t," %s",f.title);
		XDrawString(display,w,small_gc,0,CURY_OFFs,t,strlen(t));
		XTextProperty windowName;
		sprintf(t,"%s - %s",f.wildtxt,cur_dir);
		char *nameit[]={t};
		XStringListToTextProperty(nameit,1,&windowName);
		XSetWMName(display,f.base,&windowName);
	}
	for(i=0;i<f.nwin;i++) {
		if(w==f.w[i]) {
			i0=i+f.n0;
			if(i0>=f.n) {
				XDrawString(display,w,small_gc,5,CURY_OFFs," ",1);
			} else {
				if(i0<my_ff.ndirs) {
					sprintf(t,"<>%s",my_ff.dirnames[i0]);
				} else {
					sprintf(t,"%s",my_ff.filenames[i0-my_ff.ndirs]);
				}
				XDrawString(display,w,small_gc,5,CURY_OFFs,t,strlen(t));
			}
		}
	}
}


static void do_box_button(BoxList *b, Window w) {
	int i,n=b->nwin;
	if(b->xuse==0) {
		return;
	}
	if(w==b->close) {
		destroy_box(b);
		return;
	}
	if(w==b->ok || w==b->go) {
		load_entire_box(b);
	}
	if(w==b->cancel) {
		redraw_entire_box(b);
	}
	if(w==b->go) {
		run_now();
	}
	if(w==b->def && b->type==PARAMETER_BOX) {
		set_default_params();
	}
	if(w==b->def && b->type==INITIAL_CONDITION_BOX) {
		set_default_ics();
	}
	/* now for the "scrolling" */
	if(w==b->up) {
		box_list_scroll(b,1);
	}
	if(w==b->pgup) {
		box_list_scroll(b,b->nwin);
	}
	if(w==b->dn) {
		box_list_scroll(b,-1);
	}
	if(w==b->pgdn) {
		box_list_scroll(b,-b->nwin);
	}
	for(i=0;i<n;i++) {
		if(w==b->we[i]) {
			XSetInputFocus(display,w,RevertToParent,CurrentTime);
			check_box_cursor();
			HotBoxItem=i;
			HotBox=b;
			draw_editable(w,b->value[i+b->n0],b->off[i+b->n0],b->pos[i+b->n0],b->mc);
		}
	}

	if(b->type==INITIAL_CONDITION_BOX) {
		for(i=0;i<b->nwin;i++) {
			if(w==b->ck[i]) {
				b->isck[i+b->n0]=1-b->isck[i+b->n0];
				if(b->isck[i+b->n0]) {
					XDrawString(display,w,small_gc,0,CURY_OFFs,"*",1);
				} else {
					XClearWindow(display,w);
				}
			}
		}
		if(w==b->xvt) {
			set_up_xvt();
		}
		if(w==b->pp) {
			set_up_pp();
		}
		if(w==b->arr) {
			set_up_arry();
		}
	}
}


static void do_box_key(BoxList *b, XEvent ev, int *used) {
	Window w=ev.xkey.window;
	char ch;
	Window focus;
	int rev,n=b->nwin,i,j,flag;
	*used=0;
	if(b->xuse==0) {
		return;
	}
	for(i=0;i<n;i++) {
		if(b->we[i]==w) {
			XGetInputFocus(display,&focus,&rev);
			if(w==focus) {
				*used=1;
				ch=get_key_press(&ev);
				flag=edit_bitem(b,i,ch);
				if(flag==EDIT_NEXT && n>1) {
					j=i+1;
					if(j==n) {
						j=0;
					}
					XSetInputFocus(display,b->we[j],RevertToParent,CurrentTime);
					set_value_from_box(b,i);

					HotBoxItem=j;
					draw_editable(b->we[i],b->value[i+b->n0],
							b->off[i+b->n0],b->pos[i+b->n0],b->mc);
					draw_editable(b->we[j],b->value[j+b->n0],
							b->off[j+b->n0],b->pos[j+b->n0],b->mc);
					if(b->type==PARAMETER_BOX || b->type==INITIAL_CONDITION_BOX) {
						reset_sliders();
					}
				}
				if(flag==EDIT_DONE) {
					HotBoxItem=-1;
					XSetInputFocus(display,main_win,RevertToParent,CurrentTime);
					load_entire_box(b);
				}
				if(flag==EDIT_ESC) {
					HotBoxItem=-1;
					XSetInputFocus(display,main_win,RevertToParent,CurrentTime);
				}
			}
		}
	}
}


static int do_file_select_events(void) {
	int done;
	XEvent ev;
	/* plintf("Xup=%d\n",Xup); */
	while(1) {
		XNextEvent(display,&ev);
		switch(ev.type) {
		case ConfigureNotify:
		case Expose:
		case MapNotify:
			if(Xup) {
				do_expose(ev);
			}
			expose_selector(ev.xany.window);
			break;
		case ButtonPress:
			done=button_selector(ev.xbutton.window);
			if(done==1) {
				return 1; /* OK made a selection */
			}
			if(done==2) {
				return 0; /* canceled the whole thing */
			}
			break;
		case EnterNotify:
			crossing_selector(ev.xcrossing.window,0);
			break;
		case LeaveNotify:
			crossing_selector(ev.xcrossing.window,1);
			break;
		case KeyPress:
			done=selector_key(ev);
			if(done==2) {
				return 0;
			}
			if(done==1) {
				return 1;
			}
			break;
		}
	}
}


static void do_slide_button(int w, PAR_SLIDER *p) {
	static char *n[]={"*3Par/Var","Value","Low","High"};
	char values[4][MAX_LEN_SBOX];
	int status;
	double lo,hi,val;
	if(w==p->go && p->use==1) {
		run_now();
	}
	if(w!=p->top) {
		return;
	}
	strcpy(values[0],p->parname);
	sprintf(values[1],"%.16g",p->val);
	sprintf(values[2],"%.16g",p->lo);
	sprintf(values[3],"%.16g",p->hi);
	status=do_string_box(4,4,1,"Set Sliders",n,values,35);
	if(status==0) {
		return;
	}
	if(strlen(values[0])==0) { /* empty string cancels */
		p->use=0;
		return;
	}
	status=find_user_name(PARAMETER_BOX,values[0]);
	if(status==-1) {
		status=find_user_name(INITIAL_CONDITION_BOX,values[0]);
		if(status==-1) {
			err_msg("Not a parameter or variable !");
			return;
		}
		p->type=INITIAL_CONDITION_BOX;
		p->index=status;
	} else {
		p->type=PARAMETER_BOX;
		p->index=status;
	}
	lo=atof(values[2]);
	hi=atof(values[3]);
	val=atof(values[1]);
	if(val<lo || val>hi || hi<=lo) {
		err_msg(" low <= value <= high ");
		return;
	}
	p->val=val;
	p->hi=hi;
	p->lo=lo;
	strcpy(p->parname,values[0]);
	set_val(p->parname,val);
	if(p->type==INITIAL_CONDITION_BOX) {
		last_ic[p->index]=val;
	}
	redraw_params();
	redraw_ics();
	p->use=1;
	set_slide_pos(p);
	redraw_slide(p);
}


static void do_slide_motion(Window w, int x, PAR_SLIDER *p,int state) {
	int sp=SuppressBounds;
	if(w==p->slide) {
		p->pos=x;
		if(x<2) {
			p->pos=2;
		}
		if(x>(p->l-2)) {
			p->pos=p->l-2;
		}
		expose_slider(p->slide,p);
		if(p->use) {
			p->val=p->lo+ (p->hi-p->lo)*(double)(p->pos-2)/(double)(p->l-4);
			expose_slider(p->top,p);
			set_val(p->parname,p->val);
			if(p->type==INITIAL_CONDITION_BOX) {
				last_ic[p->index]=p->val;
			}
			if(state<300) {
				clr_all_scrns();
				redraw_dfield();
				create_new_cline();
				draw_label(draw_win);
				SuppressBounds=1;
				run_now();
				SuppressBounds=sp;}
		}
	}
}


static void do_slide_release(int w, PAR_SLIDER *p) {
	if(p->use==0) {
		return;
	}
	if(p->slide==w) {
		set_val(p->parname,p->val);
		if(p->type==INITIAL_CONDITION_BOX) {
			last_ic[p->index]=p->val;
		}
		redraw_ics();
		redraw_params();
	}
}


/* cursor position in letters to the left */
/* first character of string is off */
static void draw_editable(Window win, char *string, int off, int cursor, int mc) {
	int l=strlen(string)-off,rev,cp;
	Window focus;
	if(l>mc) {
		l=mc;
	}
	XClearWindow(display,win);
	XDrawString(display,win,small_gc,0,CURY_OFF,string+off,l);
	XGetInputFocus(display,&focus,&rev);
	if(focus==win) {
		cp=DCURXs*(cursor-off); /* must be fixed */
		put_edit_cursor(win,cp);
	}
}


static void draw_slider(Window w, int x, int hgt, int l) {
	int x0=x-2,i;
	if(x0<0) {
		x0=0;
	}
	if(x0>(l-4)) {
		x0=l-4;
	}
	XClearWindow(display,w);
	for(i=0;i<4;i++) {
		XDrawLine(display,w,small_gc,x0+i,0,x0+i,hgt);
	}
}


static int edit_bitem(BoxList *b, int i, int ch) {
	Window win=b->we[i];
	int i0=i+b->n0;
	char *string=b->value[i0];
	int off=b->off[i0];
	int pos=b->pos[i0];
	int mc=b->mc;
	int l=strlen(string),wpos=pos-off;
	switch(ch) {
	case LEFT:
		if(pos>0) {
			pos--;
			wpos--;
			if(wpos<0) {
				off=off-4;
				if(off<0) {
					off=0;
				}
				wpos=pos-off;
			}
		} else {
			ping();
		}
		break;
	case RIGHT:
		if(pos<l) {
			pos++;
			wpos++;
			if(wpos>mc) {
				off=off+4;
				if(off+mc>l) {
					off=l-mc;
				}
				wpos=pos-off;
			}
		} else {
			ping();
		}
		break;
	case HOME:
		pos=0;
		wpos=0;
		break;
	case END:
		pos=l;
		wpos=mc;
		break;
	case BADKEY:
		return 0;
	case DOWN: box_list_scroll(b,-1);
		return 0;
	case UP:
		box_list_scroll(b,1);
		return 0;
	case PGUP:
		box_list_scroll(b,b->nwin);
		return 0;
	case PGDOWN:
		box_list_scroll(b,-b->nwin);
		return 0;    /* junk key  */
	case ESC:
		return EDIT_ESC;
	case FINE:
		return EDIT_NEXT;
	case BACKSPACE:
	case DELETE:
		if(pos>0) {
			memmov(&string[pos-1],&string[pos],l-pos+1);
			pos--;
			wpos--;
			if(wpos<0) {
				off=off-4;
				if(off<0) {
					off=0;
				}
				wpos=pos-off;
			}
			l--;
		} else {
			ping();
		}
		break;
	case TAB:
		return EDIT_DONE;
	default:
		if( (ch>=' ') && (ch <= '~')) {
			if(strlen(string)>=256) {
				ping();
			} else {
				movmem(&string[pos+1],&string[pos],l-pos+1);
				string[pos]=ch;
				pos=pos+1;
				wpos++;
				l++;
				if(wpos>mc) {
					off=off+4;
					if(off+mc>l) {
						off=l-mc;
					}
					wpos=pos-off;
				}
			}
		}
		break;
	}
	/* all done lets save everything */
	off=pos-wpos;

	b->off[i0]=off;
	b->pos[i0]=pos;
	draw_editable(win,string,off,pos,mc);
	return 0;
}


static int edit_fitem(int ch, char *string, Window w, int *off1, int *pos1, int mc) {
	int l=strlen(string),cp;
	int off=*off1,pos=*pos1,wpos=pos-off;
	/*  plintf(" pos=%d off=%d s=%s \n",
	  pos,off,string); */
	switch(ch) {
	case LEFT:
		if(pos>0) {
			pos--;
			wpos--;
			if(wpos<0) {
				off=off-4;
				if(off<0) {
					off=0;
				}
				wpos=pos-off;
			}
		} else {
			ping();
		}
		break;
	case RIGHT:
		if(pos<l) {
			pos++;
			wpos++;
			if(wpos>mc) {
				off=off+4;
				if(off+mc>l) {
					off=l-mc;
				}
				wpos=pos-off;
			}
		} else {
			ping();
		}
		break;
	case HOME:
		pos=0;
		wpos=0;
		break;
	case END:
		pos=l;
		wpos=mc;
		break;
	case BADKEY:
		return EDIT_WAIT;
	case DOWN:
		fs_scroll(-1);
		return EDIT_WAIT;
	case UP:
		fs_scroll(1);
		return EDIT_WAIT;
	case PGUP:
		fs_scroll(filesel.nwin);
		return EDIT_WAIT;
	case PGDOWN:
		fs_scroll(-filesel.nwin);
		return EDIT_WAIT;    /* junk key  */
	case ESC:
		return EDIT_ESC;
	case FINE:
		return EDIT_DONE;
	case BACKSPACE:
	case DELETE:
		if(pos>0) {
			memmov(&string[pos-1],&string[pos],l-pos+1);
			pos--;
			wpos--;
			if(wpos<0) {
				off=off-4;
				if(off<0) {
					off=0;
				}
				wpos=pos-off;
			}
			l--;
		} else {
			ping();
		}
		break;
	case TAB:/*TAB completion of file names */ {
		struct dirent *dp;
		char ft[XPP_MAX_NAME];
		char ftpath[XPP_MAX_NAME];

		/*User may have typed ahead (maybe they remember the path they want)"*/
		/*Try to change to that new directory if it is one.*/
		if((dp=(struct dirent*)opendir(filesel.filetxt))!=NULL) {
			if(strcmp(cur_dir,filesel.filetxt)!=0) {
				change_directory(filesel.filetxt);
				get_directory(cur_dir);
				redraw_directory();
				free_finfo(&my_ff); /* delete the old file info */
				filesel.n0=0; /* back to the top of the list */
				get_fileinfo(filesel.wildtxt,cur_dir,&my_ff);
				filesel.n=my_ff.ndirs+my_ff.nfiles;
				strcpy(filesel.filetxt,cur_dir);
				int m = strlen(filesel.filetxt);
				if(filesel.filetxt[m-1] != '/') {
					strcat(filesel.filetxt,"/");
				}
				redraw_fs_text(filesel.filetxt,filesel.file,1);
				redraw_file_list();
				XFlush(display);
			}
			return(EDIT_WAIT); /*Wait for further instruction...*/
		}

		int m=strlen(filesel.filetxt)+1;

		if(m>1) {
			strcpy(ft,filesel.filetxt);
		} else {/*We are already at the root level of file system!*/
			strcpy(filesel.filetxt,"/");
			strcpy(ft,filesel.filetxt);
		}
		while ((dp=(struct dirent*)opendir(ft))==NULL)	 {
			/*This WHILE is perhaps a bit clunky but since we can't be sure of
	path separator user will type in the box a trial-by-error
	approach may be more robust.*/
			ft[m] = '\0';
			m--;
			if((ft[m] != '/') & (ft[m] != '\\')) {
				ft[m] = '\0';
				m--;
			}
			if(m==0) {
				break;
			}
		}
		int n=0;
		ft[0] = '\0';
		if(m>  strlen(filesel.filetxt)) {
			return(EDIT_WAIT);
		}
		for (n=0;n<strlen(filesel.filetxt)-m;n++) {
			ft[n]=filesel.filetxt[m+n+1];
			ft[n+1] = '\0';
		}
		strcat(ft,"*");
		strcpy(ftpath,filesel.filetxt);
		ftpath[m+1]='\0';

		/*Make sure we are in the correct directory now
	 since user could have moved cursor back up _several_
	 branches in directory tree before hitting tab key.*/
		change_directory(ftpath);
		free_finfo(&my_ff);
		filesel.n0=0;
		get_fileinfo_tab(ft,ftpath,&my_ff,filesel.wildtxt);
		filesel.n=my_ff.ndirs+my_ff.nfiles;
		if((my_ff.ndirs+my_ff.nfiles)==1) {
			if(my_ff.ndirs == 1) {/*Only possible directory -- take it.*/
				change_directory(my_ff.dirnames[0]);
				get_directory(cur_dir);
				redraw_directory();
				free_finfo(&my_ff); /* delete the old file info */
				filesel.n0=0; /* back to the top of the list */
				get_fileinfo(filesel.wildtxt,cur_dir,&my_ff);
				filesel.n=my_ff.ndirs+my_ff.nfiles;
				strcpy(filesel.filetxt,cur_dir);
				int m = strlen(filesel.filetxt);
				if(filesel.filetxt[m-1] != '/') {
					strcat(filesel.filetxt,"/");
				}
				redraw_fs_text(filesel.filetxt,filesel.file,1);
			} else {
				/*Must be that (my_ff.nfiles == 1) -- but best wait for user to confim
		 * they actually want this file by clicking on it.	*/
				/*Copy the only remaining choice into the file box and wait for user to
		make final choice.
		*/
				strcpy(filesel.filetxt,cur_dir);
				strcat(filesel.filetxt,"/");
				strcat(filesel.filetxt,my_ff.filenames[0]);
				redraw_fs_text(filesel.filetxt,filesel.file,1);
			}
		} else if(filesel.n > 1) {/*Expand the file text to most generic currently represented*/
			int j=0;
			char U[DEFAULT_STRING_LENGTH];
			if(my_ff.ndirs > 0) {
				strcpy(U,my_ff.dirnames[0]);
			} else {
				strcpy(U,my_ff.filenames[0]);
			}
			for (j=0;j<filesel.n;j++) {
				if(j<my_ff.ndirs) {
					stringintersect(U,my_ff.dirnames[j]);
				} else {
					stringintersect(U,my_ff.filenames[j-my_ff.ndirs]);
				}
				if(strlen(U) == 0) {/*No common substring*/
					break;
				}
			}
			strcpy(filesel.filetxt,ftpath);
			strcat(filesel.filetxt,U);
			/*Actually don't want to force appending of path separator here
		since we want user to decide between matching directory and files.
		*/
			redraw_fs_text(filesel.filetxt,filesel.file,1);
		}
		redraw_file_list();
		XFlush(display);
		return 0;
	}
	default:
		if( (ch>=' ') && (ch <= '~')) {
			if(strlen(string)>=256) {
				ping();
			} else {
				movmem(&string[pos+1],&string[pos],l-pos+1);
				string[pos]=ch;
				pos=pos+1;
				wpos++;
				l++;
				if(wpos>mc) {
					off=off+4;
					if(off+mc>l) {
						off=l-mc;
					}
					wpos=pos-off;
				}
			}
		}
		break;
	}
	/* all done lets save everything */
	off=pos-wpos;
	*off1=off;
	*pos1=pos;
	XClearWindow(display,w);
	XDrawString(display,w,small_gc,0,CURY_OFF,string+off,strlen(string)-off);
	cp=DCURXs*(pos-off);
	put_edit_cursor(w,cp);
	return 0;
}


static void enter_slider(Window w, PAR_SLIDER *p, int val) {
	if(w==p->top || w==p->go) {
		XSetWindowBorderWidth(display,w,val+1);
	}
}


static void expose_selector(Window w) {
	display_file_sel(filesel,w);
}


static void expose_slider(Window w, PAR_SLIDER *p) {
	int x,len=12*DCURXs;
	char top[DEFAULT_STRING_LENGTH];
	if(w==p->slide) {
		draw_slider(w,p->pos,p->hgt,p->l);
		return;
	}
	if(w==p->go) {
		XDrawString(display,w,small_gc,2,0.75*CURY_OFFs,"go",2);
		return;
	}
	if(p->use) {
		if(w==p->left) {
			sprintf(top,"%.16g",p->lo);
			x=1;
			XClearWindow(display,w);
			XDrawString(display,w,small_gc,x,CURY_OFFs,top,strlen(top));
			return;
		}
		if(w==p->right) {
			sprintf(top,"%.16g",p->hi);
			x=1;
			if(strlen(top)<12) {
				x=len-DCURXs*strlen(top)-1;
			}
			XClearWindow(display,w);
			XDrawString(display,w,small_gc,x,CURY_OFFs,top,strlen(top));
			return;
		}
		if(w==p->top) {
			sprintf(top,"%s=%.16g",p->parname,p->val);
			XClearWindow(display,w);
			XDrawString(display,w,small_gc,2,CURY_OFFs,top,strlen(top));
		}
	} else {
		if(w==p->top) {
			sprintf(top,"Par/Var?");
			x=1;
			XClearWindow(display,w);
			XDrawString(display,w,small_gc,x,CURY_OFFs,top,strlen(top));
		}
	}
}


static void fs_scroll(int i) {
	int n0=filesel.n0;
	int new_int,nend;
	int nw=filesel.nwin,n=filesel.n;
	if(n<=nw) {
		return;
	}
	new_int=n0-i;
	nend=new_int+nw;
	if(new_int<0) {
		new_int=0;
	}
	if(nend>n) {
		new_int=n-nw;
	}
	filesel.n0=new_int;
	redraw_file_list();
}


/* this returns the fixed width, the number of entries allowed */
static void get_nrow_from_hgt(int h, int *n, int *w) {
	int hgt=DCURYs+4;
	*w=28*DCURXs;
	*n= h/(hgt+4)-5;
}


static void justify_string(Window w1, char *s1) {
	int n1=strlen(s1)*DCURXs,nt=10*DCURXs;
	int i=0;
	if(n1<nt) {
		i=nt-n1;
	}
	XClearWindow(display,w1);
	XDrawString(display,w1,small_gc,i,CURY_OFFs,s1,strlen(s1));
}


static void load_entire_box(BoxList *b) {
	int i,n=b->n;

	for(i=0;i<n;i++) {
		set_value_from_box(b,i);
	}
	if(b->type==PARAMETER_BOX) {
		re_evaluate_kernels();
		redo_all_fun_tables();
		reset_sliders();
	}
	if(b->type==DELAY_BOX) {
		do_init_delay(DELAY);
	}
}


static void make_box_list(BoxList *b, char *wname, char *iname, int n, int type, int use) {
	int nrow,i;
	char sss[XPP_MAX_NAME];
	double z;
	nrow=10;
	if(n<10) {
		nrow=n;
	}
	b->xuse=0;
	b->use=use;
	b->mc=21;
	b->type=type;
	b->n=n;
	b->n0=0;
	b->nwin=nrow;
	b->value=(char **)malloc(n*sizeof(char*));
	b->pos=(int *)malloc(n*sizeof(int));
	b->off=(int *)malloc(n*sizeof(int));
	b->iname=(char *)malloc(strlen(iname)+5);
	strcpy(b->iname,iname);
	b->wname=(char *)malloc(strlen(wname)+5);
	strcpy(b->wname,wname);
	for(i=0;i<n;i++) {
		b->value[i]=(char *)malloc(256);
		switch(type) {
		case PARAMETER_BOX:
			get_val(upar_names[i],&z);
			sprintf(sss,"%.16g",z);
			set_edit_params(b,i,sss);
			break;
		case INITIAL_CONDITION_BOX:
			sprintf(sss,"%.16g",last_ic[i]);
			set_edit_params(b,i,sss);
			break;
		case BOUNDARY_CONDITION_BOX:
			set_edit_params(b,i,my_bc[i].string);
			break;
		case DELAY_BOX:
			set_edit_params(b,i,delay_string[i]);
			break;
		}
	}
}


static void make_box_list_window(BoxList *b,int type) {
	int nrow,n;
	int x,y;
	int xb1,xb2,xb3,xb4;
	int i,wid1,wid2;
	int width,height,wid,hgt;

	Window base;
	XTextProperty winname,iconame;
	XSizeHints size_hints;
	n=b->n;
	nrow=b->nwin;

	wid1=16*DCURXs;
	wid2=22*DCURXs;
	wid=wid1+wid2+DCURXs;
	hgt=DCURYs+4;
	height=(nrow+4)*(hgt+4)+2*hgt;
	width=wid+8*DCURXs;
	b->minwid=width;
	b->minhgt=height;
	base=make_plain_window(RootWindow(display,screen),0,0,width,height,4);
	b->base=base;
	XStringListToTextProperty(&b->wname,1,&winname);
	XStringListToTextProperty(&b->iname,1,&iconame);
	size_hints.flags=PPosition|PSize|PMinSize;
	size_hints.x=0;
	size_hints.y=0;
	size_hints.width=width;
	size_hints.height=height;
	size_hints.min_width=width;
	size_hints.min_height=height;
	size_hints.max_width=width;
	size_hints.max_height=height;
	XClassHint class_hints;
	class_hints.res_name="";
	class_hints.res_class="";
	XSetWMProperties(display,base,&winname,&iconame,NULL,0,&size_hints,NULL,&class_hints);
	b->w = (Window *)malloc(nrow*sizeof(Window));
	b->we = (Window *)malloc(nrow*sizeof(Window));
	if(type==INITIAL_CONDITION_BOX) {
		b->ck=(Window *)malloc(nrow*sizeof(Window));
		b->isck=(int *)malloc(n*sizeof(int));
		for(i=0;i<n;i++) {
			b->isck[i]=0;
		}
	}
	xb1 = 2+7*DCURXs+14;
	xb2=xb1+7*DCURXs+14;
	xb3=xb2+7*DCURXs+14;
	xb4=xb3+7*DCURXs+14;

	b->close=make_window(base,2,5,7*DCURXs+10,DCURYs,1);
	b->ok=make_window(base,xb1,5,7*DCURXs+10,DCURYs,1);
	b->def=make_window(base,xb2,5,7*DCURXs+10,DCURYs,1);
	b->cancel=make_window(base,xb3,5,7*DCURXs+10,DCURYs,1);
	b->go=make_window(base,xb4,5,7*DCURXs+10,DCURYs,1);
	xb1=DCURXs+wid1+wid2+12;

	b->up=make_icon_window(base,xb1,1.75*DCURYs+24+3,32,24,1,0,0,lineup_bits);
	b->dn=make_icon_window(base,xb1,1.75*DCURYs+48+6,32,24,1,0,0,linedn_bits);
	b->pgup=make_icon_window(base,xb1,1.75*DCURYs,32,24,1,0,0,pageup_bits);
	b->pgdn=make_icon_window(base,xb1,1.75*DCURYs+72+9,32,24,1,0,0,pagedn_bits);

	for(i=0;i<nrow;i++) {
		x=DCURXs;
		y=DCURYs+(hgt+4)*i+1.5*hgt;
		b->w[i]=make_plain_window(base,x,y,wid1,hgt,0);
		b->we[i]=make_plain_window(base,x+wid1+2,y,wid2,hgt,1);
		XSelectInput(display,b->w[i],BOXEVENT);
		if(type==INITIAL_CONDITION_BOX) {
			b->ck[i]=make_plain_window(base,1,y,6,DCURYs,1);
		}
	}
	y=DCURYs+(hgt+4)*nrow+1.5*hgt;
	x=(width-24)/3;
	if(type==INITIAL_CONDITION_BOX) {
		b->xvt=make_window(base,x,y,5*DCURXs,DCURYs,1);
		b->pp=make_window(base,x+6*DCURXs,y,5*DCURXs,DCURYs,1);
		b->arr=make_window(base,x+12*DCURXs,y,5*DCURXs,DCURYs,1);
	}
	b->xuse=1;
}


static void make_par_slider(Window base, int x, int y, int width, int index, int var_index) {
	int mainhgt=3*(DCURYs+2);
	int mainwid=32*DCURXs;
	int xs;
	Window w;
	if(mainwid<(width+4)) {
		mainwid=width+4;
	}
	w=make_plain_window(base,x,y,mainwid,mainhgt,1);
	my_par_slide[index].main=w;
	xs=(mainwid-width-4)/2;
	my_par_slide[index].slide=make_window(w,xs,DCURYs+5,width+4,DCURYs-4,1);
	my_par_slide[index].go=make_window(w,xs+width+8,DCURYs+5,3*DCURXs,DCURYs-3,1);
	my_par_slide[index].top=make_window(w,2,2,mainwid-6,DCURYs,1);
	my_par_slide[index].left=make_window(w,2,2*DCURYs+3,12*DCURXs,DCURYs,0);
	my_par_slide[index].right=make_window(w,mainwid-12*DCURXs-4,2*DCURYs+3,
										  12*DCURXs,DCURYs,0);
	my_par_slide[index].lo=0.0;
	my_par_slide[index].hi=1.0;
	my_par_slide[index].val=0.5;
	my_par_slide[index].use=0;
	my_par_slide[index].l=width+4;
	my_par_slide[index].pos=(width+4)/2;
	my_par_slide[index].parname[0]=0;
	my_par_slide[index].hgt=DCURYs-4;

	if (var_index >= 0) {
		strcpy(my_par_slide[index].parname, SLIDERVAR[var_index]);
		my_par_slide[index].use = 1;
		my_par_slide[index].lo = SLIDERLO[var_index];
		my_par_slide[index].hi = SLIDERHI[var_index];
		get_val(my_par_slide[index].parname, &my_par_slide[index].val);
	}
}


static void new_wild(void) {
	free_finfo(&my_ff); /* delete the old file info */
	filesel.n0=0; /* back to the top of the list */
	get_fileinfo(filesel.wildtxt,cur_dir,&my_ff);
	filesel.n=my_ff.ndirs+my_ff.nfiles;
	redraw_file_list();
	XFlush(display);
}


static void put_edit_cursor(Window w, int pos) {
	int x1=pos;
	int x2=x1+1;
	XDrawLine(display,w,small_gc,x1,1,x1,DCURYs-1);
	XDrawLine(display,w,small_gc,x2,1,x2,DCURYs-1);
}


static void redraw_directory(void) {
	XClearWindow(display,filesel.dir);
	expose_selector(filesel.dir);
}


static void redraw_entire_box(BoxList *b) {
	if(b->xuse==0) {
		return;
	}
	switch(b->type) {
	case PARAMETER_BOX:
		redraw_params();
		return;
	case BOUNDARY_CONDITION_BOX:
		redraw_bcs();
		return;
	case INITIAL_CONDITION_BOX:
		redraw_ics();
		return;
	case DELAY_BOX:
		redraw_delays();
		return;
	}
}


static void redraw_file_list(void) {
	int i;
	for(i=0;i<filesel.nwin;i++) {
		XClearWindow(display,filesel.w[i]);
		expose_selector(filesel.w[i]);
	}
}


static void redraw_fs_text(char *string, Window w, int flag) {
	XClearWindow(display,w);
	filesel.off=0;
	if(flag) {
		filesel.pos=strlen(string);
	}
	XDrawString(display,w,small_gc,0,CURY_OFF,string,strlen(string));
	if(flag) {
		put_edit_cursor(w,DCURXs*strlen(string));
	}
}


static void redraw_slide(PAR_SLIDER *p) {
	expose_slider(p->slide,p);
	expose_slider(p->top,p);
	expose_slider(p->left,p);
	expose_slider(p->right,p);
}


static int selector_key(XEvent ev) {
	char ch;
	int flag;
	/* plintf(" hot=%d pos=%d \n",filesel.hot,filesel.pos); */
	ch=get_key_press(&ev);
	switch(filesel.hot) {
	case HOTFILE:
		flag=edit_fitem(ch,filesel.filetxt,filesel.file,
						&(filesel.off),&(filesel.pos),29);
		if(flag==EDIT_DONE) {
			return 1;
		}
		if(flag==EDIT_ESC) {
			return 2;
		}
		return(0);
	case HOTWILD:
		flag=edit_fitem(ch,filesel.wildtxt,filesel.wild,
						&(filesel.off),&(filesel.pos),29);
		if(flag==EDIT_DONE) {
			new_wild();
			return 0;
		}
		if(flag==EDIT_ESC) {
			return 2;
		}
		return 0;
	}
	return 0;
}


static void set_default_ics(void) {
	int i;
	for(i=0;i<NODE+NMarkov;i++) {
		last_ic[i]=default_ic[i];
	}
	redraw_ics();
}


static void set_default_params(void) {
	int i;
	char junk[DEFAULT_STRING_LENGTH];
	for(i=0;i<NUPAR;i++) {
		set_val(upar_names[i],default_val[i]);
		sprintf(junk,"%.16g",default_val[i]);
		set_edit_params(&ParamBox,i,junk);
	}
	redraw_params();
	re_evaluate_kernels();
	redo_all_fun_tables();
}


static void set_slide_pos(PAR_SLIDER *p) {
	double pos;
	int ip;
	pos=2. + (p->l-4)*(p->val-p->lo)/(p->hi-p->lo);
	ip=(int)pos;
	if(ip<2) {
		ip=2;
	}
	if(ip>(p->l-2)) {
		ip=p->l-2;
	}
	p->pos=ip;
}


static void set_up_arry(void) {
	int i;
	int plot_list[2],n=0;

	for(i=0;i<ICBox.n;i++) {
		if(ICBox.isck[i]) {
			if(n<2) {
				plot_list[n]=i+1;
				n++;
			}
			ICBox.isck[i]=0;
		}
	}
	for(i=0;i<ICBox.nwin;i++) {
		XClearWindow(display,ICBox.ck[i]);
	}
	if(n==2) {
		optimize_aplot(plot_list);
	}
}


static void set_up_pp(void) {
	int i;
	int plot_list[3],n=0;

	for(i=0;i<ICBox.n;i++) {
		if(ICBox.isck[i]) {
			if(n<3) {
				plot_list[n]=i+1;
				n++;
			}
			ICBox.isck[i]=0;
		}
	}
	for(i=0;i<ICBox.nwin;i++) {
		XClearWindow(display,ICBox.ck[i]);
	}
	if(n>1) {
		graph_all(plot_list,n,1);
	}
}


static void set_up_xvt(void) {
	int i;
	int plot_list[10];
	int n=0;
	for(i=0;i<ICBox.n;i++) {
		if(ICBox.isck[i]) {
			if(n<10) {
				plot_list[n]=i+1;
				n++;
			}
			ICBox.isck[i]=0;
		}
	}
	for(i=0;i<ICBox.nwin;i++) {
		XClearWindow(display,ICBox.ck[i]);
	}
	if(n>0) {
		graph_all(plot_list,n,0);
	}
}


static void set_value_from_box(BoxList *b, int i) {
	char *s;
	double z;
	s=b->value[i];
	switch(b->type) {
	case INITIAL_CONDITION_BOX:
		if(to_float(s,&z)==-1) {
			return;
		}
		last_ic[i]=z;
		add_edit_float(b,i,z);
		break;
	case PARAMETER_BOX:
		if(to_float(s,&z)==-1) {
			return;
		}
		set_val(upar_names[i],z);
		add_edit_float(b,i,z);
		break;
	case BOUNDARY_CONDITION_BOX:
		strcpy(my_bc[i].string,s);
		add_editval(b,i,s);
		break;
	case DELAY_BOX:
		strcpy(delay_string[i],s);
		add_editval(b,i,s);
		break;
	}
}


static int to_float(char *s, double *z) {
	int flag;
	*z=0.0;
	if(s[0]=='%') {
		flag=do_calc(&s[1],z);
		if(flag==-1) {
			return -1;
		}
		return 0;
	}
	*z=atof(s);
	return(0);
}
