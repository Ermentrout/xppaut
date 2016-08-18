#include "browse.h"

#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "aniparse.h"
#include "comline.h"
#include "dialog_box.h"
#include "eig_list.h"
#include "form_ode.h"
#include "ggets.h"
#include "init_conds.h"
#include "integrate.h"
#include "load_eqn.h"
#include "main.h"
#include "many_pops.h"
#include "markov.h"
#include "menus.h"
#include "menudrive.h"
#include "mykeydef.h"
#include "newhome.h"
#include "parserslow.h"
#include "pop_list.h"
#include "storage.h"
#include "util/strutil.h"
#include "util/timeutil.h"
#include "xpplim.h"

#include "bitmap/browse.bitmap"

/* --- Forward Declarations --- */
static int add_stor_col(char *name, char *formula, BROWSER *b);
static Window br_button_data(Window root, int row, int col, char *name, int iflag);
static void	browse_but_on(BROWSER *b, int i, Window w, int yn);
static void browse_button(XEvent ev, BROWSER *b);
static void browse_keypress(XEvent ev, int *used, BROWSER *b);
static int check_for_stor(float **data);
static void chk_seq(char *f, int *seq, double *a1, double *a2);
static void data_add_col(BROWSER *b);
static void data_del_col(BROWSER *b);
static void data_down(BROWSER *b);
static void data_end(BROWSER *b);
static void data_find(BROWSER *b);
static void data_first(BROWSER *b);
static void data_get(BROWSER *b);
static void data_home(BROWSER *b);
static void data_last(BROWSER *b);
static void data_left(BROWSER *b);
static void data_pgdn(BROWSER *b);
static void data_pgup(BROWSER *b);
static void data_read(BROWSER *b);
static void data_replace(BROWSER *b);
static void data_restore(BROWSER *b);
static void data_right(BROWSER *b);
static void data_table(BROWSER *b);
static void data_unreplace(BROWSER *b);
static void data_up(BROWSER *b);
static void data_write(BROWSER *b);
static void del_stor_col(char *var, BROWSER *b);
static void display_browser(Window w, BROWSER b);
static void draw_data(BROWSER b);
static void enter_browser(XEvent ev, BROWSER *b, int yn);
static void expose_browser(XEvent ev, BROWSER b);
static void find_value(int col, double val, int *row, BROWSER b);
static void kill_browser(BROWSER *b);
static void make_browser(BROWSER *b, char *wname, char *iname, int row, int col);
static void make_d_table(double xlo, double xhi, int col, char *filename, BROWSER b);
static void redraw_browser(BROWSER b);
static void replace_column(char *var, char *form, float **dat, int n);
static void resize_browser(Window win, BROWSER *b);
static void unreplace_column(void);
static void write_browser_data(FILE *fp, BROWSER *b);

/* --- Data --- */
BROWSER my_browser;
static float *old_rep;
static int REPLACE=0,R_COL=0;

float **storage;

/* --- Functions --- */
Window br_button(Window root, int row, int col, char *name, int iflag) {
	Window win;
	int dcol=12*DCURXs;
	int drow=(DCURYs+6);
	/*int width=strlen(name)*DCURXs;
  */
	int width=8*DCURXs;
	int x;
	int y;
	if(iflag==1) {
		dcol=14*DCURXs;
	}
	x=dcol*col+4;
	y=drow*row+4;
	win=make_window(root,x,y,width+5,DCURYs+1,1);
	XSelectInput(display,win,MYMASK);
	return(win);
}


void data_get_mybrowser(int row) {
	my_browser.row0=row;
	data_get(&my_browser);
}

void enter_my_browser(XEvent ev, int yn) {
	if(my_browser.xflag==0) {
		return;
	}
	enter_browser(ev,&my_browser,yn);
}


/*   These are the global exporters ...   */
void expose_my_browser(XEvent ev) {
	if(my_browser.xflag==0) {
		return;
	}
	expose_browser(ev,my_browser);
}


void find_variable(char *s, int *col) {
	*col=-1;
	if(strcasecmp("T",s)==0) {
		*col=0;
		return;
	}
	*col=find_user_name(2,s);
	if(*col>-1) {
		*col=*col+1;
	}
}


float **get_browser_data(void) {
	return my_browser.data;
}


float *get_data_col(int c) {
	return my_browser.data[c];
}


void get_data_xyz(float *x, float *y, float *z,
				  int i1, int i2, int i3, int off) {
	int in=my_browser.row0+off;
	*x=my_browser.data[i1][in];
	*y=my_browser.data[i2][in];
	*z=my_browser.data[i3][in];
}


int get_maxrow_browser(void) {
	return my_browser.maxrow;
}


void init_browser(void) {
	my_browser.dataflag=0;
	my_browser.data=storage;
	my_browser.maxcol=NEQ+1;
	my_browser.maxrow=0;
	my_browser.col0=1;
	my_browser.row0=0;
	my_browser.istart=0;
	my_browser.iend=0;
	strcpy(my_browser.hinttxt,"hint");
}


void make_new_browser(void) {
	if(my_browser.xflag==1) {
		XRaiseWindow(display,my_browser.base);
		return;
	}
	make_browser(&my_browser,"Data Viewer","Data",20,5);
	my_browser.xflag=1;
}


void my_browse_button(XEvent ev) {
	if(my_browser.xflag==0) {
		return;
	}
	browse_button(ev,&my_browser);
}


void my_browse_keypress(XEvent ev, int *used) {
	if(my_browser.xflag==0) {
		return;
	}
	browse_keypress(ev,used,&my_browser);
}


void open_write_file(FILE **fp, char *fil, int *ok) {
	char ans;
	*ok=0;
	*fp=fopen(fil,"r");
	if(*fp!=NULL) {
		fclose(*fp);
		ans=(char)TwoChoice("Yes","No",
							"File Exists! Overwrite?","yn");
		if(ans!='y') {
			return;
		}
	}

	*fp=fopen(fil,"w");
	if(*fp==NULL) {
		respond_box("Ok",
					"Cannot open file");
		*ok=0;
	} else {
		*ok=1;
	}
	return;
}


void refresh_browser(int length) {
	my_browser.dataflag=1;
	my_browser.maxrow=length;
	my_browser.iend=length;
	if(Xup && my_browser.xflag==1) {
		draw_data(my_browser);
	}
}


void reset_browser(void) {
	my_browser.maxrow=0;
	my_browser.dataflag=0;
}


void resize_my_browser(Window win) {
	if(my_browser.xflag==0) {
		return;
	}
	resize_browser(win,&my_browser);
}


void set_browser_data(float **data,int col0) {
	my_browser.data=data;
	my_browser.col0=col0;
}


void wipe_rep(void) {
	if(!REPLACE) {
		return;
	}
	free(old_rep);
	REPLACE=0;
}


void write_mybrowser_data(FILE *fp) {
	write_browser_data(fp,&my_browser);
}


/* --- Static functions --- */
static int add_stor_col(char *name, char *formula,BROWSER *b) {
	int com[4000],i,j;

	if(add_expr(formula,com,&i)) {
		err_msg("Bad Formula .... ");
		return(0);
	}
	if((my_ode[NEQ+FIX_VAR]=(int *)malloc((i+2)*sizeof(int)))==NULL) {
		err_msg("Cant allocate formula space");
		return(0);
	}
	if((storage[NEQ+1]=(float *)malloc(MAXSTOR * sizeof(float)))==NULL) {
		err_msg("Cant allocate space ....");
		free(my_ode[NEQ]);
		return(0);
	}
	if((ode_names[NEQ]=(char *)malloc(80))==NULL) {
		err_msg("Cannot allocate space ...");
		free(my_ode[NEQ]);
		free(storage[NEQ+1]);
		return(0);
	}
	strcpy(ode_names[NEQ],formula);
	strupr(ode_names[NEQ]);
	for(j=0;j<=i;j++)
		my_ode[NEQ+FIX_VAR][j]=com[j];

	strcpy(uvar_names[NEQ],name);
	strupr(uvar_names[NEQ]);
	for(i=0;i<b->maxrow;i++)
		storage[NEQ+1][i]=0.0;   /*  zero it all   */

	for(i=0;i<b->maxrow;i++) {
		for(j=0;j<NODE+1;j++) {
			set_ivar(j,(double)storage[j][i]);
		}
		for(j=NODE;j<NEQ;j++) {
			set_val(uvar_names[j],(double)storage[j+1][i]);
		}
		storage[NEQ+1][i]=(float)evaluate(com);
	}
	add_var(uvar_names[NEQ],0.0);  /*  this could be trouble .... */
	NEQ++;
	b->maxcol=NEQ+1;
	redraw_browser(*b);
	return(1);
}


static Window br_button_data(Window root, int row, int col, char *name, int iflag) {
	Window win;
	int dcol=12*DCURXs;
	int drow=(DCURYs+6);
	int width=strlen(name)*DCURXs;

	int x;
	int y;
	if(iflag==1) {
		dcol=14*DCURXs;
	}
	x=dcol*col+4;
	y=drow*row+4;
	win=make_window(root,x,y,width+5,DCURYs+1,1);
	XSelectInput(display,win,MYMASK);
	return(win);
}

/*  if button is pressed in the browser
	then do the following  */
static void browse_button(XEvent ev, BROWSER *b) {
	XEvent zz;
	int done=1;
	Window w=ev.xbutton.window;
	if(my_browser.xflag==0) {
		return;
	}
	if(w==b->up || w==b->down || w==b->pgup || w==b->pgdn || w==b->left || w==b->right) {
		done=1;
		while(done) {
			if(w==b->up) {
				data_up(b);
			}
			if(w==b->down) {
				data_down(b);
			}
			if(w==b->pgup) {
				data_pgup(b);
			}
			if(w==b->pgdn) {
				data_pgdn(b);
			}
			if(w==b->left) {
				data_left(b);
			}
			if(w==b->right) {
				data_right(b);
			}
			waitasec(100);
			if(XPending(display)>0) {
				XNextEvent(display,&zz);
				switch(zz.type) {
				case ButtonRelease:
					done=0;
					break;
				}
			}
		}
		return;
	}

	if(w==b->home) {
		data_home(b);
		return;
	}
	if(w==b->end) {
		data_end(b);
		return;
	}
	if(w==b->first) {
		data_first(b);
		return;
	}
	if(w==b->last) {
		data_last(b);
		return;
	}
	if(w==b->restore) {
		data_restore(b);
		return;
	}
	if(w==b->write) {
		data_write(b);
		return;
	}
	if(w==b->get) {
		data_get(b);
		return;
	}
	if(w==b->find) {
		data_find(b);
		return;
	}
	if(w==b->repl) {
		data_replace(b);
		return;
	}
	if(w==b->load) {
		data_read(b);
		return;
	}
	if(w==b->addcol) {
		data_add_col(b);
		return;
	}
	if(w==b->delcol) {
		data_del_col(b);
		return;
	}
	if(w==b->unrepl) {
		data_unreplace(b);
		return;
	}
	if(w==b->table) {
		data_table(b);
		return;
	}
	if(w==b->close) {
		kill_browser(b);
		return;
	}
}


static void browse_but_on(BROWSER *b, int i, Window w, int yn) {
	int val=1;
	if(yn) {
		val=2;
	}
	XSetWindowBorderWidth(display,w,val);
	if(yn && TipsFlag && i>=0) {
		strcpy(b->hinttxt,browse_hint[i]);
		display_browser(b->hint,*b);
	}
}


static void browse_keypress(XEvent ev, int *used, BROWSER *b) {
	Window w=ev.xkey.window;
	char ks;
	Window w2;
	int rev;
	*used=0;
	if(my_browser.xflag==0) {
		return;
	}
	XGetInputFocus(display,&w2,&rev);

	if(w==b->main || w==b->base || w==b->upper || w2==b->base) {
		*used=1;
		ks=(char)get_key_press(&ev);
		if(ks==UP) {
			data_up(b);
			return;
		}
		if(ks==DOWN) {
			data_down(b);
			return;
		}
		if(ks==PGUP) {
			data_pgup(b);
			return;
		}
		if(ks==PGDOWN) {
			data_pgdn(b);
			return;
		}
		if(ks==LEFT) {
			data_left(b);
			return;
		}
		if(ks==RIGHT) {
			data_right(b);
			return;
		}
		if(ks==HOME) {
			data_home(b);
			return;
		}
		if(ks==END) {
			data_end(b);
			return;
		}
		if(ks=='s' || ks=='S') {
			data_first(b);
			return;
		}
		if(ks=='e' || ks=='E') {
			data_last(b);
			return;
		}
		if(ks=='r' || ks=='R') {
			data_restore(b);
			return;
		}
		if(ks=='W' || ks=='w') {
			data_write(b);
			return;
		}
		if(ks=='g' || ks=='G') {
			data_get(b);
			return;
		}
		if(ks=='f' || ks=='F') {
			data_find(b);
			return;
		}
		if(ks=='l' || ks=='L') {
			data_read(b);
			return;
		}
		if(ks=='u' || ks=='U') {
			data_unreplace(b);
			return;
		}
		if(ks=='t' || ks=='T') {
			data_table(b);
			return;
		}
		if(ks=='p' || ks=='P') {
			data_replace(b);
			return;
		}
		if(ks=='a' || ks=='A') {
			data_add_col(b);
			return;
		}
		if(ks=='d' || ks=='D') {
			data_del_col(b);
			return;
		}
		if(ks==ESC) {
			XSetInputFocus(display,command_pop,
						   RevertToParent,CurrentTime);
			return;
		}
	} /* end of cases */
}


static int check_for_stor(float **data) {
	if(data!=storage) {
		err_msg("Only data can be in browser");
		return(0);
	}
	else return(1);
}


static void chk_seq(char *f,int *seq, double *a1, double *a2) {
	int i,j=-1;
	char n1[DEFAULT_STRING_LENGTH],n2[256];
	int n=strlen(f);
	*seq=0;
	*a1=0.0;
	*a2=0.0;
	for(i=0;i<n;i++)
	{
		if(f[i]==':') {
			*seq=1;
			j=i;
		}

		if(f[i]==';') {
			*seq=2;
			j=i;
		}
	}
	if(j>-1) {
		for(i=0;i<j;i++) {
			n1[i]=f[i];
		}
		n1[j]=0;
		for(i=j+1;i<n;i++) {
			n2[i-j-1]=f[i];
		}
		n2[n-j-1]=0;
		*a1=atof(n1);
		*a2=atof(n2);
	}
	/*      plintf("seq=%d a1=%g a2=%g\n",*seq,*a1,*a2); */
}


static void data_add_col(BROWSER *b) {
	Window w;
	int rev,status;
	char var[20],form[80];
	if(check_for_stor(b->data)==0) {
		return;
	}
	XGetInputFocus(display,&w,&rev);
	strcpy(var,"");
	strcpy(form,"");
	status=get_dialog("Add Column","Name",var,"Ok","Cancel",20);
	if(status!=FORGET_ALL) {
		status=get_dialog("Add Column","Formula:",form,"Add it","Cancel",80);
		if(status!=FORGET_ALL) {
			add_stor_col(var,form,b);
		}
	}
}


/*  this only works with storage  */
static void data_del_col(BROWSER *b) {
	Window w;
	int rev,status;
	char var[20];
	if(check_for_stor(b->data)==0) {
		return;
	}
	XGetInputFocus(display,&w,&rev);
	err_msg("Sorry - not working very well yet...");
	return;
	strcpy(var,"");
	status=get_dialog("Delete","Name",var,"Ok","Cancel",20);
	if(status!=FORGET_ALL) {
		del_stor_col(var,b);
	}
}


static void data_down(BROWSER *b) {
	if(b->row0<(b->maxrow-1)) {
		b->row0++;
		draw_data(*b);
	}
}


static void data_end(BROWSER *b) {
	b->row0=b->maxrow-1;
	draw_data(*b);
}


static void data_find(BROWSER *b) {
	Window w;
	int rev,status;

	static char *name[]={"*0Variable","Value"};
	char value[2][25];
	int col,row=0;
	float val;

	sprintf(value[0],"%s",uvar_names[0]);
	sprintf(value[1],"0.00");
	XGetInputFocus(display,&w,&rev);
	status=do_string_box(2,2,1,"Find Data",name,value,40);

	XSetInputFocus(display,w,rev,CurrentTime);

	if(status==0) {
		return;
	}
	val=(float)atof(value[1]);
	find_variable(value[0],&col);
	if(col>=0) {
		find_value(col,val,&row,*b);
	}
	if(row>=0) {
		b->row0=row;
		draw_data(*b);
	}
}


static void data_first(BROWSER *b) {
	b->istart=b->row0;
}


static void data_get(BROWSER *b) {
	int i,in=b->row0;
	set_ivar(0,(double)storage[0][in]);
	for(i=0;i<NODE;i++) {
		last_ic[i]=(double)storage[i+1][in];
		set_ivar(i+1,last_ic[i]);
	}
	for(i=0;i<NMarkov;i++) {
		last_ic[i+NODE]=(double)storage[i+NODE+1][in];
		set_ivar(i+1+NODE+FIX_VAR,last_ic[i+NODE]);
	}
	for(i=NODE+NMarkov;i<NEQ;i++) {
		set_val(uvar_names[i],storage[i+1][in]);
	}
	redraw_ics();
}


static void data_home(BROWSER *b) {
	b->row0=0;
	b->istart=0;
	b->iend=b->maxrow;
	draw_data(*b);
}


static void data_last(BROWSER *b) {
	b->iend=b->row0+1;
}


static void data_left(BROWSER *b) {
	int i=b->col0;
	if(i>1) {
		b->col0--;
		redraw_browser(*b);
	}
}


static void data_pgdn(BROWSER *b) {
	int i=b->row0+b->nrow;
	if(i<(b->maxrow-1)) {
		b->row0=i;
	} else {
		b->row0=b->maxrow-1;
	}
	draw_data(*b);
}


static void data_pgup(BROWSER *b) {
	int i=b->row0-b->nrow;
	if(i>0) {
		b->row0=i;
	} else {
		b->row0=0;
	}
	draw_data(*b);
}


static void data_read(BROWSER *b) {
	int status;
	char fil[DEFAULT_STRING_LENGTH];
	char ch;
	FILE *fp;
	int k;
	int len,count=0,white=1;
	float z;

	strcpy(fil,"test.dat");
	status=file_selector("Load data",fil,"*.dat");
	if(status==0) {
		return;
	}
	fp=fopen(fil,"r");
	if(fp==NULL) {
		respond_box("Ok", "Cannot open file");
		return;
	}
	/*  Now we establish the width of the file and read it.
	  If there are more columns than available we
	  ignore them.

	 if there are fewer rows we read whats necessary
	 if there are more rows then read until we
	 are done or MAX_STOR_ROW.
	 This data can be plotted etc like anything else
	*/
	do {
		if(fscanf(fp,"%c",&ch)<1) {
			plintf("Error reading file %s", fp);
		}
		if( !isspace(ch) && (white)) {
			white=0;
			++count;
		}
		if(isspace(ch)  && (1-white)) {
			white=1;
		}
	} while(ch != '\n');
	rewind(fp);
	len=0;
	while(!feof(fp)) {
		for(k=0;k<count;k++) {
			if(fscanf(fp,"%f ",&z)<1) {
				plintf("Error reading file %s", fp);
			}
			if(k<b->maxcol) {
				b->data[k][len]=z;
			}
		}
		++len;
		if(len>=MAXSTOR) {
			break;
		}
	}
	fclose(fp);
	refresh_browser(len);
	storind=len;
}


static void data_replace(BROWSER *b) {
	Window w;
	int rev,status;
	char var[20],form[80];
	XGetInputFocus(display,&w,&rev);
	strcpy(var,uvar_names[0]);
	strcpy(form,uvar_names[0]);
	status=get_dialog("Replace","Variable:",var,"Ok","Cancel",20);
	if(status!=FORGET_ALL) {
		status=get_dialog("Replace","Formula:",form,"Replace","Cancel",80);
		if(status!=FORGET_ALL) {
			replace_column(var,form,b->data,b->maxrow);
		}
		draw_data(*b);
	}
	XSetInputFocus(display,w,rev,CurrentTime);
}


static void data_restore(BROWSER *b)  {
	restore(b->istart,b->iend);
}


static void data_right(BROWSER *b) {
	int i=b->col0+b->ncol;
	if(i<=b->maxcol) {
		b->col0++;
		redraw_browser(*b);
	}
}


static void data_table(BROWSER *b) {
	Window w;
	int rev,status;

	static char *name[]={"Variable","Xlo","Xhi","File"};
	char value[4][25];

	double xlo=0,xhi=1;
	int col;
	sprintf(value[0],"%s",uvar_names[0]);
	sprintf(value[1],"0.00");
	sprintf(value[2],"1.00");
	sprintf(value[3],"%s.tab",value[0]);
	XGetInputFocus(display,&w,&rev);
	status=do_string_box(4,4,1,"Tabulate",name,value,40);
	XSetInputFocus(display,w,rev,CurrentTime);
	if(status==0) {
		return;
	}
	xlo=atof(value[1]);
	xhi=atof(value[2]);
	find_variable(value[0],&col);
	if(col>=0) {
		make_d_table(xlo,xhi,col,value[3],*b);
	}
}


static void data_unreplace(BROWSER *b) {
	unreplace_column();
	draw_data(*b);
}


static void data_up( BROWSER *b) {
	if(b->row0>0) {
		b->row0--;
		draw_data(*b);
	}
}


static void data_write(BROWSER *b) {
	int status;
	char fil[DEFAULT_STRING_LENGTH];
	FILE *fp;
	int i,j;
	int ok;

	strcpy(fil,"test.dat");
	status=file_selector("Write data",fil,"*.dat");
	if(status==0) {
		return;
	}
	open_write_file(&fp,fil,&ok);
	if(!ok) {
		return;
	}
	for(i=b->istart;i<b->iend;i++) {
		for(j=0;j<b->maxcol;j++) {
			fprintf(fp,"%.8g ",b->data[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


static void del_stor_col(char *var, BROWSER *b) {
	int nc;
	int i,j;

	find_variable(var,&nc);

	if(nc<0) {
		err_msg("No such column....");
		return;
	}
	if(nc<=NEQ_MIN) { /* NEQ_MIN = NODE+NAUX */
		err_msg("Can't delete that column");
		return;
	}
	if(check_active_plot(nc)==1) {
		err_msg("This variable is still actively plotted! - Cant delete!");
		return;
	}
	/*  plintf(" nc=%d NEQ= %d\n",nc,NEQ); */
	change_plot_vars(nc);
	if(nc<NEQ) {
		for(j=nc;j<NEQ;j++) {
			for(i=0;i<b->maxrow;i++) {
				storage[j][i]=storage[j+1][i];
			}
			for(i=0;i<400;i++) {
				my_ode[j-1+FIX_VAR][i]=my_ode[j+FIX_VAR][i];
			}
			strcpy(uvar_names[j-1],uvar_names[j]);
			strcpy(ode_names[j-1],ode_names[j]);
		}
	}
	free(storage[NEQ+1]);
	free(ode_names[NEQ]);
	free(my_ode[NEQ+FIX_VAR]);
	NEQ--;
	b->maxcol=NEQ+1;
	redraw_browser(*b);
}


static void display_browser(Window w, BROWSER b) {
	int i,i0;
	if(w==b.hint) {
		XClearWindow(display,b.hint);
		XDrawString(display,w,small_gc,8,CURY_OFFs,b.hinttxt,strlen(b.hinttxt));
		return;
	}

	if(w==b.find) {
		xds("Find")
	}
	if(w==b.up) {
		xds("Up")
	}
	if(w==b.down) {
		xds("Down")
	}
	if(w==b.pgup) {
		xds("PgUp")
	}
	if(w==b.pgdn) {
		xds("PgDn")
	}
	if(w==b.left) {
		xds("Left")
	}
	if(w==b.right) {
		xds("Right")
	}
	if(w==b.home) {
		xds("Home")
	}
	if(w==b.end) {
		xds("End")
	}
	if(w==b.first) {
		xds("First")
	}
	if(w==b.last) {
		xds("Last")
	}
	if(w==b.restore) {
		xds("Restore")
	}
	if(w==b.write) {
		xds("Write")
	}
	if(w==b.get) {
		xds("Get")
	}
	if(w==b.repl) {
		xds("Replace")
	}
	if(w==b.unrepl) {
		xds("Unrepl")
	}
	if(w==b.table) {
		xds("Table")
	}
	if(w==b.load) {
		xds("Load")
	}
	if(w==b.time) {
		xds("Time")
	}
	if(w==b.addcol) {
		xds("Add col")
	}
	if(w==b.close) {
		xds("Close")
	}
	if(w==b.delcol) {
		xds("Del col")
	}
	for(i=0;i<BMAXCOL;i++) {
		if(w==b.label[i]) {
			i0=i+b.col0-1;
			if(i0<b.maxcol-1) {
				XDrawString(display,w,small_gc,5,CURY_OFFs,
							uvar_names[i0],strlen(uvar_names[i0]));
			}
		}
	}
	if(w==b.main) {
		draw_data(b);
	}
}



static void draw_data(BROWSER b) {
	int i,i0,j,j0;
	int x0;
	char string[50];
	int dcol=DCURXs*14;
	int drow=(DCURYs+6);
	if(b.dataflag==0) {
		return;  /*   no data  */
	}
	XClearWindow(display,b.main);

	/* Do time data first  */
	for(i=0;i<b.nrow;i++) {
		i0=i+b.row0;
		if(i0<b.maxrow) {
			sprintf(string,"%.8g",b.data[0][i0]);
			XDrawString(display,b.main,small_gc,DCURXs/2+5,i*drow+DCURYs
						,string,strlen(string));
		}
	}

	/* Do data stuff   */
	for(j=0;j<b.ncol;j++) {
		x0=(j+1)*dcol+DCURXs/2;
		j0=j+b.col0;
		if(j0>=b.maxcol) {
			return;  /* if this one is too big, they all are  */
		}
		for(i=0;i<b.nrow;i++) {
			i0=i+b.row0;
			if(i0<b.maxrow) {
				sprintf(string,"%.7g",b.data[j0][i0]);
				XDrawString(display,b.main,small_gc,x0+5,i*drow+DCURYs,
							string,strlen(string));
			}
		}
	}
}


static void enter_browser(XEvent ev, BROWSER *b, int yn) {
	Window w=ev.xexpose.window;
	if(w==b->find) {
		browse_but_on(b,0,w,yn);
	}
	if(w==b->up) {
		browse_but_on(b,1,w,yn);
	}
	if(w==b->down) {
		browse_but_on(b,2,w,yn);
	}
	if(w==b->pgup) {
		browse_but_on(b,3,w,yn);
	}
	if(w==b->pgdn) {
		browse_but_on(b,4,w,yn);
	}
	if(w==b->left) {
		browse_but_on(b,5,w,yn);
	}
	if(w==b->right) {
		browse_but_on(b,6,w,yn);
	}
	if(w==b->home) {
		browse_but_on(b,7,w,yn);
	}
	if(w==b->end) {
		browse_but_on(b,8,w,yn);
	}
	if(w==b->first) {
		browse_but_on(b,9,w,yn);
	}
	if(w==b->last) {
		browse_but_on(b,10,w,yn);
	}
	if(w==b->restore) {
		browse_but_on(b,11,w,yn);
	}
	if(w==b->write) {
		browse_but_on(b,12,w,yn);
	}
	if(w==b->get) {
		browse_but_on(b,13,w,yn);
	}
	if(w==b->repl) {
		browse_but_on(b,14,w,yn);
	}
	if(w==b->unrepl) {
		browse_but_on(b,15,w,yn);
	}
	if(w==b->table) {
		browse_but_on(b,16,w,yn);
	}
	if(w==b->load) {
		browse_but_on(b,17,w,yn);
	}
	if(w==b->time) {
		browse_but_on(b,18,w,yn);
	}
	if(w==b->addcol) {
		browse_but_on(b,19,w,yn);
	}
	if(w==b->delcol) {
		browse_but_on(b,20,w,yn);
	}
	if(w==b->close) {
		browse_but_on(b,-1,w,yn);
	}
}


static void expose_browser(XEvent ev, BROWSER b) {
	if(my_browser.xflag==0) {
		return;
	}
	if(ev.type!=Expose) {
		return;
	}
	display_browser(ev.xexpose.window,b);
}


static void find_value(int col, double val, int *row, BROWSER b) {
	int n=b.maxrow;
	int i;
	int ihot=0;
	float err,errm;
	errm=fabs(b.data[col][0]-val);
	for(i=b.row0;i<n;i++) {
		err=fabs(b.data[col][i]-val);
		if(err<errm) {
			ihot=i;
			errm=err;
		}
	}
	*row=ihot;
}


static void kill_browser(BROWSER *b) {
	b->xflag=0;
	waitasec(ClickTime);
	XDestroySubwindows(display,b->base);
	XDestroyWindow(display,b->base);
}


static void make_browser(BROWSER *b, char *wname, char *iname, int row, int col) {
	int i;
	int ncol=col;
	int width,height;
	Window base;
	/* XWMHints wm_hints;
*/
	XTextProperty winname,iconname;
	XSizeHints size_hints;
	int dcol=DCURXs*17;
	int drow=(DCURYs+6);
	int ystart=8;

	if(ncol<5) {
		ncol=5;
	}
	height=drow*(row+6);
	width=ncol*dcol;
	b->nrow=row;
	b->ncol=ncol;
	base=make_plain_window(RootWindow(display,screen),0,0,width,height,4);
	b->base=base;
	XSelectInput(display,base,ExposureMask|KeyPressMask|ButtonPressMask|
				 StructureNotifyMask);
	/* plintf("Browser base: %d \n",base); */
	XStringListToTextProperty(&wname,1,&winname);
	XStringListToTextProperty(&iname,1,&iconname);

	size_hints.flags=PPosition|PSize|PMinSize;
	size_hints.x=0;
	size_hints.y=0;
	/* size_hints.width=width;
 size_hints.height=height; */
	size_hints.min_width=width-15;
	size_hints.min_height=height;
	/* wm_hints.initial_state=IconicState;
 wm_hints.flags=StateHint;
 */
	XClassHint class_hints;
	class_hints.res_name="";
	class_hints.res_class="";

	XSetWMProperties(display,base,&winname,&iconname,NULL,0,&size_hints,NULL,&class_hints);
	make_icon((char*)browse_bits,browse_width,browse_height,base);
	b->upper=make_window(base,0,0,width,ystart+drow*6,1);
	XSetWindowBackground(display,b->upper,MyMainWinColor);
	b->main=make_plain_window(base,0,ystart+drow*6,width,row*drow,1);
	XSetWindowBackground(display,b->main,MyDrawWinColor);
	b->find=br_button(base,0,0,"find",0);
	b->get=br_button(base,1,0,"get ",0);
	b->repl=br_button(base,2,0,"replace",0);
	b->restore=br_button(base,0,1,"restore",0);
	b->write=br_button(base,1,1," write ",0);
	b->load=br_button(base,2,1," load  ",0);
	b->first=br_button(base,0,2,"first",0);
	b->last=br_button(base,1,2,"last ",0);
	b->unrepl=br_button(base,2,2,"unrepl",0);
	b->table=br_button(base,2,3,"table",0);
	b->up=br_button(base,0,3," up ",0);
	b->down=br_button(base,1,3,"down",0);
	b->pgup=br_button(base,0,4,"pgup",0);
	b->pgdn=br_button(base,1,4,"pgdn",0);
	b->left=br_button(base,0,5,"left ",0);
	b->right=br_button(base,1,5,"right",0);
	b->home=br_button(base,0,6,"home",0);
	b->end=br_button(base,1,6,"end ",0);
	b->addcol=br_button(base,2,4,"addcol",0);
	b->delcol=br_button(base,2,5,"delcol",0);
	b->close=br_button(base,2,6,"close",0);
	b->time=br_button(base,5,0,"time ",1);
	b->hint=make_window(base,0,4*drow,width-17,drow-3,1);
	XSelectInput(display,b->time,SIMPMASK);

	for(i=0;i<BMAXCOL;i++) {
		/*  plintf("%d ",i); */
		b->label[i]=br_button_data(base,5,i+1,"1234567890",1);
		XSelectInput(display,b->label[i],SIMPMASK);
		/* plintf(" %d \n",i); */
	}
	if(noicon==0) {
		XIconifyWindow(display,base,screen);
	}
}


static void make_d_table(double xlo, double xhi, int col, char *filename, BROWSER b) {
	int i,npts,ok;
	FILE *fp;
	open_write_file(&fp,filename,&ok);
	if(!ok) {
		return;
	}
	npts=b.iend-b.istart;
	fprintf(fp,"%d\n",npts);
	fprintf(fp,"%g\n%g\n",xlo,xhi);
	for(i=0;i<npts;i++)
		fprintf(fp,"%10.10g\n",b.data[col][i+b.istart]);

	fclose(fp);
	ping();
}


static void redraw_browser(BROWSER b) {
	int i,i0;
	Window w;
	draw_data(b);
	for(i=0;i<BMAXCOL;i++) {
		w=b.label[i];
		i0=i+b.col0-1;
		if(i0<(b.maxcol-1)) {
			XClearWindow(display,w);
			XDrawString(display,w,small_gc,5,CURY_OFFs,uvar_names[i0],
						strlen(uvar_names[i0]));
		}
	}
}


static void replace_column(char *var, char *form, float **dat, int n) {
	int com[200],i,j;
	int intflag=0;
	int dif_var=-1;
	int seq=0;
	double a1,a2,da=0.0;
	float old=0.0,dt,derv=0.0;
	float sum=0.0;
	if(n<2) {
		return;
	}
	dt=NJMP*DELTA_T;
	/* first check for derivative or integral symbol */
	i=0;
	while(i<strlen(form)) {
		if(!isspace(form[i])) {
			break;
		}
		i++;
	}
	if(form[i]=='&') {
		intflag=1;
		form[i]=' ';
	}
	if(form[i]=='@') {
		form[i]=' ';
		find_variable(form,&dif_var);
		if(dif_var<0) {
			err_msg("No such variable");
			return;
		}
	}

	if(dif_var<0) {
		chk_seq(form,&seq,&a1,&a2);
	}
	if(seq==1) {
		if(a1==a2) {
			seq=3;
		} else {
			da=(a2-a1)/((double)(n-1));
		}
	}
	if(seq==2) {
		da=a2;
	}
	if(seq==3) {
		err_msg("Illegal sequence");
		return;
	}

	/*  first compile formula ... */
	if(dif_var<0 && seq==0) {
		if(add_expr(form,com,&i)) {
			NCON=NCON_START;
			NSYM=NSYM_START;
			err_msg("Illegal formula...");
			return;
		}
	}
	/* next check to see if column is known ... */

	find_variable(var,&i);
	if(i<0) {
		err_msg("No such column...");
		NCON=NCON_START;
		NSYM=NSYM_START;
		return;
	}
	R_COL=i;

	/* Okay the formula is cool so lets allocate and replace  */

	wipe_rep();
	old_rep=(float *)malloc(sizeof(float)*n);
	REPLACE=1;
	for(i=0;i<n;i++) {
		old_rep[i]=dat[R_COL][i];
		if(dif_var<0) {
			if(seq==0) {
				for(j=0;j<NODE+1;j++) {
					set_ivar(j,(double)dat[j][i]);
				}
				for(j=NODE;j<NEQ;j++) {
					set_val(uvar_names[j],(double)dat[j+1][i]);
				}
				if(intflag) {
					sum+=(float)evaluate(com);
					dat[R_COL][i]=sum*dt;
				} else {
					dat[R_COL][i]=(float)evaluate(com);
				}
			} else {
				dat[R_COL][i]=(float)(a1+i*da);
			}
		} else {
			if(i==0) {
				derv=(dat[dif_var][1]-dat[dif_var][0])/dt;
			}
			if(i==(n-1)) {
				derv=(dat[dif_var][i]-old)/dt;
			}
			if(i>0 && i<(n-1)) {
				derv=(dat[dif_var][i+1]-dat[dif_var][i])/dt;
			}
			old=dat[dif_var][i];
			dat[R_COL][i]=derv;
		}
	}
	NCON=NCON_START;
	NSYM=NSYM_START;
}


static void resize_browser(Window win, BROWSER *b) {
	unsigned int w,h,hreal;
	int dcol=17*DCURXs,drow=DCURYs+6;
	int i0;
	int newrow,newcol;
	if(my_browser.xflag==0) {
		return;
	}
	if(win!=b->base) {
		return;
	}
	get_new_size(win,&w,&h);
	hreal=h;

	/* first make sure the size is is ok and an integral value
	 * of the proper width and height
	*/
	i0=w/dcol;
	if((w%dcol)>0) {
		i0++;
	}
	if(i0>b->maxcol) {
		i0=b->maxcol;
	}
	w=i0*dcol;
	if(i0<5) {
		w=5*dcol;
	}
	newcol=i0;
	h=hreal-8-5*drow;
	i0=h/drow;
	if((h%drow)>0) {
		i0++;
	}
	if(i0>b->maxrow) {
		i0=b->maxrow;
	}
	h=i0*drow+DCURXs/2;
	newrow=i0;
	/*  Now resize everything   */
	if(b->ncol==newcol && b->nrow==newrow) {
		return;
	}
	b->ncol=newcol;
	b->nrow=newrow;

	XResizeWindow(display,b->base,w-17,hreal);
	XResizeWindow(display,b->upper,w-17,8+drow*3);
	XResizeWindow(display,b->main,w-17,h);

	/* Let the browser know how many rows and columns of data  */
}


static void unreplace_column(void) {
	int i,n=my_browser.maxrow;
	if(!REPLACE) {
		return;
	}
	for(i=0;i<n;i++) {
		my_browser.data[R_COL][i]=old_rep[i];
	}
	wipe_rep();
}


static void write_browser_data(FILE *fp, BROWSER *b) {
	int i,j,l;
	for(i=b->istart;i<b->iend;i++) {
		if(N_plist>0) {
			for(l=0;l<N_plist;l++) {
				j=plotlist[l];
				fprintf(fp,"%.8g ",b->data[j][i]);
			}
		} else {
			for(j=0;j<b->maxcol;j++) {
				fprintf(fp,"%.8g ",b->data[j][i]);
			}
		}
		fprintf(fp,"\n");
	}
}
