#include "graf_par.h"

#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "aniparse.h"
#include "arrayplot.h"
#include "auto_nox.h"
#include "auto_x11.h"
#include "axes2.h"
#include "browse.h"
#include "color.h"
#include "form_ode.h"
#include "ggets.h"
#include "graphics.h"
#include "init_conds.h"
#include "integrate.h"
#include "kinescope.h"
#include "load_eqn.h"
#include "main.h"
#include "many_pops.h"
#include "markov.h"
#include "menu.h"
#include "menudrive.h"
#include "menus.h"
#include "my_ps.h"
#include "my_svg.h"
#include "mykeydef.h"
#include "nullcline.h"
#include "pop_list.h"
#include "rubber.h"
#include "struct.h"
#include "xpplim.h"

/* --- Macros --- */
enum xppLineStyle {
	lsSTABLE_EQUILIBRIUM = 0,
	lsUNSTABLE_EQUILIBRIUM,
	lsSTABLE_PERIODIC = 8,
	lsUNSTABLE_PERIODIC = 9
};

/* --- Types --- */
typedef struct {
	char angle[20];
	char yes[3];
	double start;
	double incr;
	int nclip;
} MOV3D;

typedef struct {
	float *x[MAXBIFCRV],*y[MAXBIFCRV];
	int color[MAXBIFCRV],npts[MAXBIFCRV],nbifcrv;
	Window w;
} BD;


/* --- Forward declarations --- */
static void add_bd_crv(float *x, float *y, int len, int type, int ncrv);
static int alter_curve(char *title, int in_it, int n);
static void axes_opts(void);
static void check_flags(void);
static void check_val(double *x1, double *x2, double *xb, double *xd);
static void corner_cube(double *xlo, double *xhi, double *ylo, double *yhi);
static int create_crv(int ind);
static void create_ps(void);
static void create_svg(void);
static void delete_frz(void);
static void delete_frz_crv(int i);
static void draw_bd(Window w);
static void draw_freeze_key(void);
static void draw_frozen_cline(int index, Window w);
static void edit_curve(void);
static void edit_frz(void);
static void edit_frz_crv(int i);
static void export_graf_data(void);
static void fit_window(void);
static void free_bd(void);
static int freeze_crv(int ind);
static void frz_bd(void);
static void get_2d_view(int ind);
static void get_3d_view(int ind);
static int get_frz_index(Window w);
static void kill_frz(void);
static void movie_rot(double start, double increment, int nclip, int angle);
static void new_curve(void);
static void read_bd(FILE *fp);
static void scroll_window(void);
static void set_key(int x, int y);
static void update_view(float xlo,float xhi, float ylo, float yhi);
static void user_window(void);
static void zoom_in(int i1, int j1, int i2, int j2);
static void zoom_out(int i1, int j1, int i2, int j2);


/* --- Data --- */
static double FreezeKeyX,FreezeKeyY;
static int FreezeKeyFlag;
static int CurrentCurve=0;
static BD my_bd;
static MOV3D mov3d = { "theta","N",45,45,7};

NCLINE nclines[MAXNCLINE];
int AutoFreezeFlag=0;


/* --- Functions --- */
void add_a_curve_com(int c) {
	switch(c) {
	case 0:
		if(MyGraph->nvars>=MAXPERPLOT) {
			err_msg("Too many plots!");
			return;
		}
		new_curve();
		break;
	case 1:
		if(MyGraph->nvars>1) {
			MyGraph->nvars=MyGraph->nvars-1;
		}
		break;
	case 2:
		MyGraph->nvars=1;
		break;
	case 3:
		edit_curve();
		break;
	case 4:
		create_ps();
		break;
	case 5:
		create_svg();
		break;
		/* case 6: freeze();
	  break; */
	case 7:
		axes_opts();
		break;
	case 8:
		export_graf_data();
		break;
		/*  case 9: change_cmap();
	   break; */
	}
	check_flags();
	redraw_the_graph();
}


void auto_freeze_it(void) {
	if(AutoFreezeFlag==0) {
		return;
	}
	create_crv(0);
}


void change_cmap_com(int i) {
	NewColormap(i);
}


void change_view_com(int com) {
	if(com==2) {
		make_my_aplot("Array!");
		edit_aplot();
		return;
	}
	if(com==3) {
		new_vcr();
		return;
	}

	MyGraph->grtype=5*com;
	if(MyGraph->grtype<5) {
		get_2d_view(CurrentCurve);
	} else {
		get_3d_view(CurrentCurve);
	}
	check_flags();
	redraw_the_graph();
}


void check_windows(void) {
	double zip,zap;
	check_val(&MyGraph->xmin,&MyGraph->xmax,&MyGraph->xbar,&MyGraph->dx);
	check_val(&MyGraph->ymin,&MyGraph->ymax,&MyGraph->ybar,&MyGraph->dy);
	check_val(&MyGraph->zmin,&MyGraph->zmax,&MyGraph->zbar,&MyGraph->dz);
	check_val(&MyGraph->xlo,&MyGraph->xhi,&zip,&zap);
	check_val(&MyGraph->ylo,&MyGraph->yhi,&zip,&zap);
}


void default_window(void) {
	if(MyGraph->ThreeDFlag) {
		MyGraph->xmax=x_3d[1];
		MyGraph->ymax=y_3d[1];
		MyGraph->zmax=z_3d[1];
		MyGraph->xmin=x_3d[0];
		MyGraph->ymin=y_3d[0];
		MyGraph->zmin=z_3d[0];

		corner_cube(&(MyGraph->xlo),&(MyGraph->xhi),&(MyGraph->ylo),&(MyGraph->yhi));
		check_windows();
	} else {
		MyGraph->xmax=x_3d[1];
		MyGraph->ymax=y_3d[1];

		MyGraph->xmin=x_3d[0];
		MyGraph->ymin=y_3d[0];

		MyGraph->xlo=MyGraph->xmin;
		MyGraph->ylo=MyGraph->ymin;
		MyGraph->xhi=MyGraph->xmax;
		MyGraph->yhi=MyGraph->ymax;
		check_windows();
	}
	redraw_the_graph();
}


void draw_freeze(Window w) {
	int i,j,type=MyGraph->grtype,lt=0;
	float oldxpl,oldypl,oldzpl=0.0,xpl,ypl,zpl=0.0;
	float *xv,*yv,*zv;
	for(i=0;i<MAXNCLINE;i++) {
		draw_frozen_cline(i,w);
	}
	for(i=0;i<MAXFRZ;i++) {
		if(frz[i].use==1 && frz[i].w==w && frz[i].type==type) {
			if(frz[i].color<0) {
				set_linestyle(-frz[i].color);
				lt=1;
			} else {
				set_linestyle(frz[i].color);
			}
			xv=frz[i].xv;
			yv=frz[i].yv;
			zv=frz[i].zv;
			oldxpl=xv[0];
			oldypl=yv[0];
			if(type>0) {
				oldzpl=zv[0];
			}
			for(j=0;j<frz[i].len;j++) {
				xpl=xv[j];
				ypl=yv[j];
				if(type>0) {
					zpl=zv[j];
				}
				if(lt==0) {
					if(type==0) {
						line_abs(oldxpl,oldypl,xpl,ypl);
					} else {
						line_3d(oldxpl,oldypl,oldzpl,xpl,ypl,zpl);
					}
				} else {
					if(type==0) {
						point_abs(xpl,ypl);
					} else {
						point_3d(xpl,ypl,zpl);
					}
				}
				oldxpl=xpl;
				oldypl=ypl;
				if(type>0) {
					oldzpl=zpl;
				}
			}
		}
	}
	draw_freeze_key();
	draw_bd(w);
}


void freeze_com(int c) {
	switch(c) {
	case 0:
		freeze_crv(0);
		break;
	case 1:
		delete_frz();
		break;
	case 2:
		edit_frz();
		break;
	case 3:
		kill_frz();
		break;
		/*case 4:
   key_frz();
   break; */
	case 5:
		frz_bd();
		break;
	case 6:
		free_bd();
		break;
	case 7:
		AutoFreezeFlag=1-AutoFreezeFlag;
		break;
	}
}


void get_3d_par_com(void) {
	static char *n[]={"Persp (1=On)","ZPlane","ZView","Theta","Phi","Movie(Y/N)",
					  "Vary (theta/phi)","Start angle", "Increment",
					  "Number increments"};
	char values[10][MAX_LEN_SBOX];
	int status;

	int nclip=8,angle=0;
	double start,increment=45;

	if(MyGraph->grtype<5) {
		return;
	}

	sprintf(values[0],"%d",MyGraph->PerspFlag);
	sprintf(values[1],"%g",MyGraph->ZPlane);
	sprintf(values[2],"%g",MyGraph->ZView);
	sprintf(values[3],"%g",MyGraph->Theta);
	sprintf(values[4],"%g",MyGraph->Phi);
	sprintf(values[5],"%s",mov3d.yes);
	sprintf(values[6],"%s",mov3d.angle);
	sprintf(values[7],"%g",mov3d.start);
	sprintf(values[8],"%g",mov3d.incr);
	sprintf(values[9],"%d",mov3d.nclip);

	status=do_string_box(10,5,2,"3D Parameters",n,values,28);
	if(status!=0) {
		MyGraph->PerspFlag=atoi(values[0]);
		MyGraph->ZPlane=atof(values[1]);
		MyGraph->ZView=atof(values[2]);
		MyGraph->Theta=atof(values[3]);
		MyGraph->Phi=atof(values[4]);
		if(values[5][0]=='y' || values[5][0]=='Y') {
			strcpy(mov3d.yes,values[5]);
			strcpy(mov3d.angle,values[6]);
			start=atof(values[7]);
			increment=atof(values[8]);
			nclip=atoi(values[9]);
			mov3d.start=start;
			mov3d.incr=increment;
			mov3d.nclip=nclip;
			angle=0;
			if(mov3d.angle[0]=='p' || mov3d.angle[0]=='P') {
				angle=1;
			}
			/*     XRaiseWindow(display,MyGraph->w); */
			movie_rot(start,increment,nclip,angle);
		}
		make_rot(MyGraph->Theta,MyGraph->Phi);
		/*  Redraw the picture   */
		redraw_the_graph();
	}
}


void get_max(int index, double *vmin, double *vmax) {
	float x0,x1,z;
	double temp;
	int i;
	x0=my_browser.data[index][0];
	x1=x0;
	for(i=0;i<my_browser.maxrow;i++) {
		z=my_browser.data[index][i];
		if(z<x0) {
			x0=z;
		}
		if(z>x1) {
			x1=z;
		}
	}
	*vmin=(double)x0;
	*vmax=(double)x1;
	if(fabs(*vmin-*vmax)<REAL_SMALL) {
		temp=.05*lmax(fabs(*vmin),1.0);
		*vmin=*vmin-temp;
		*vmax=*vmax+temp;
	}
}


void graph_all(int *list, int n, int type) {
	int i;
	if(type==0) {
		for(i=0;i<n;i++) {
			MyGraph->xv[i]=0;
			MyGraph->yv[i]=list[i];
			MyGraph->line[i]=MyGraph->line[0];
			MyGraph->color[i]=i;
		}
		MyGraph->nvars=n;
		MyGraph->grtype=0;
		MyGraph->ThreeDFlag=0;
	}
	if(type==1) {
		MyGraph->nvars=1;
		MyGraph->xv[0]=list[0];
		MyGraph->yv[0]=list[1];
		MyGraph->grtype=0;
		MyGraph->ThreeDFlag=0;
		if(n==3) {
			MyGraph->zv[0]=list[2];
			MyGraph->grtype=5;
			MyGraph->ThreeDFlag=1;
		}
	}
	check_flags();
	fit_window();
}


void ind_to_sym(int ind, char *str) {
	if(ind==0) {
		strcpy(str,"T");
	} else {
		strcpy(str,uvar_names[ind-1]);
	}
}


/*  Bifurcation curve importing */
void init_bd(void) {
	my_bd.nbifcrv=0;
}


void key_frz_com(int c) {
	int x,y;
	switch(c) {
	case 0:
		FreezeKeyFlag=0;
		break;
	case 1:
		MessageBox("Position with mouse");
		if(get_mouse_xy(&x,&y,draw_win)) {
			set_key(x,y);
			draw_freeze_key();
		}
		KillMessageBox();
	}
}


void redraw_the_graph(void) {
	blank_screen(draw_win);
	set_normal_scale();
	do_axes();
	hi_lite(draw_win);
	restore(0,my_browser.maxrow);
	draw_label(draw_win);
	draw_freeze(draw_win);
	redraw_dfield();
	if(MyGraph->Nullrestore) {
		restore_nullclines();
	}
}


void window_zoom_com(int c) {
	int i1,i2,j1,j2;
	switch(c) {
	case 0:
		user_window();
		break;
	case 1:
		if(rubber(&i1,&j1,&i2,&j2,draw_win,RUBBOX)==0) {
			break;
		}
		zoom_in(i1,j1,i2,j2);
		break;
	case 2:
		if(rubber(&i1,&j1,&i2,&j2,draw_win,RUBBOX)==0) {
			break;
		}
		zoom_out(i1,j1,i2,j2);
		break;
	case 3:
		fit_window();
		break;
	case 4:
		default_window();
		break;
	case 5:
		scroll_window();
		break;
	}
	set_normal_scale();
}


/*  a short cut   */
void xi_vs_t(void) {
	char name[20],value[20];
	int i=MyGraph->yv[0];

	ind_to_sym(i,value);
	sprintf(name,"Plot vs t: ");
	new_string(name,value);
	find_variable(value,&i);

	if(i>-1) {
		MyGraph->yv[0]=i;
		MyGraph->grtype=0;
		MyGraph->xv[0]=0;
		if(storind>=2) {
			get_max(MyGraph->xv[0],&(MyGraph->xmin),&(MyGraph->xmax));
			get_max(MyGraph->yv[0],&(MyGraph->ymin),&(MyGraph->ymax));
		} else {
			MyGraph->xmin=T0;
			MyGraph->xmax=TEND;
		}
		MyGraph->xlo=MyGraph->xmin;
		MyGraph->ylo=MyGraph->ymin;
		MyGraph->xhi=MyGraph->xmax;
		MyGraph->yhi=MyGraph->ymax;
		check_windows();
		check_flags();
		set_normal_scale();
		redraw_the_graph();
	}
}


/* --- Static functions --- */
static void add_bd_crv(float *x, float *y, int len, int type, int ncrv) {
	int i;
	if(ncrv>=MAXBIFCRV) {
		return;
	}

	my_bd.x[ncrv]=(float *)malloc(sizeof(float)*len);
	my_bd.y[ncrv]=(float *)malloc(sizeof(float)*len);
	for(i=0;i<len;i++) {
		my_bd.x[ncrv][i]=x[i];
		my_bd.y[ncrv][i]=y[i];
	}
	my_bd.npts[ncrv]=len;
	i=lsSTABLE_EQUILIBRIUM;
	if(type==UNSTABLE_PERIODIC) {
		i=lsUNSTABLE_PERIODIC;
	}
	if(type==STABLE_PERIODIC) {
		i=lsSTABLE_PERIODIC;
	}
	if(type==UNSTABLE_EQUILIBRIUM) {
		i=lsUNSTABLE_EQUILIBRIUM;
	}
	my_bd.color[ncrv]=i;
}


static int alter_curve(char *title, int in_it, int n) {
	static char *nn[]={"*0X-axis","*0Y-axis","*0Z-axis","*4Color","Line type"};
	char values[5][MAX_LEN_SBOX];
	int status,i;
	int i1=MyGraph->xv[in_it],i2=MyGraph->yv[in_it],i3=MyGraph->zv[in_it];
	char n1[15],n2[15],n3[15];

	ind_to_sym(i1,n1);
	ind_to_sym(i2,n2);
	ind_to_sym(i3,n3);
	sprintf(values[0],"%s",n1);
	sprintf(values[1],"%s",n2);
	sprintf(values[2],"%s",n3);
	sprintf(values[3],"%d",MyGraph->color[in_it]);
	sprintf(values[4],"%d",MyGraph->line[in_it]);
	status=do_string_box(5,5,1,title,nn,values,25);
	if(status!=0) {
		find_variable(values[0],&i);
		if(i>-1) {
			MyGraph->xv[n]=i;
		}
		find_variable(values[1],&i);
		if(i>-1) {
			MyGraph->yv[n]=i;
		}
		find_variable(values[2],&i);
		if(i>-1) {
			MyGraph->zv[n]=i;
		}
		MyGraph->line[n]=atoi(values[4]);
		i=atoi(values[3]);
		if(i<0 || i>10) {
			i=0;
		}
		MyGraph->color[n]=i;
		return(1);
	}
	return(0);
}


static void axes_opts(void) {
	static char *n[]={"X-origin","Y-origin","Z-origin",
					  "X-org(1=on)","Y-org(1=on)","Z-org(1=on",
					  "PSFontSize"};
	char values[7][MAX_LEN_SBOX];
	int status;
	sprintf(values[0],"%g",MyGraph->xorg);
	sprintf(values[1],"%g",MyGraph->yorg);
	sprintf(values[2],"%g",MyGraph->zorg);
	sprintf(values[3],"%d",MyGraph->xorgflag);
	sprintf(values[4],"%d",MyGraph->yorgflag);
	sprintf(values[5],"%d",MyGraph->zorgflag);
	sprintf(values[6],"%d",PS_FontSize);
	status=do_string_box(7,7,1,"Axes options",n,values,25);
	if(status!=0) {
		MyGraph->xorg=atof(values[0]);
		MyGraph->yorg=atof(values[1]);
		MyGraph->zorg=atof(values[2]);
		MyGraph->xorgflag=atoi(values[3]);
		MyGraph->yorgflag=atoi(values[4]);
		MyGraph->zorgflag=atoi(values[5]);
		PS_FontSize=atoi(values[6]);
		redraw_the_graph();
	}
}


static void check_flags(void) {
	if(MyGraph->grtype>4) {
		MyGraph->ThreeDFlag=1;
	} else {
		MyGraph->ThreeDFlag=0;
	}
	if((MyGraph->xv[0]==0) || (MyGraph->yv[0]==0) ||
	   ((MyGraph->zv[0]==0) && (MyGraph->ThreeDFlag==1))) {
		MyGraph->TimeFlag=1;
	} else {
		MyGraph->TimeFlag=0;
	}
}


/* see get_max for details */
static void check_val(double *x1, double *x2, double *xb, double *xd) {
	double temp;
	if(*x1==*x2) {
		temp=.05*lmax(fabs(*x1),1.0);
		*x1=*x1-temp;
		*x2=*x2+temp;
	}
	if(*x1>*x2) {
		temp=*x2;
		*x2=*x1;
		*x1=temp;

	}
	*xb=.5*(*x1+*x2);
	*xd=2.0/(*x2-*x1);
}


static void corner_cube(double *xlo, double *xhi, double *ylo, double *yhi) {
	float x,y;
	float x1,x2,y1,y2;
	threedproj(-1.,-1.,-1.,&x,&y);
	x1=x;
	x2=x;
	y1=y;
	y2=y;
	threedproj(-1.,-1.,1.,&x,&y);
	if(x<x1) {
		x1=x;
	}
	if(x>x2) {
		x2=x;
	}
	if(y<y1) {
		y1=y;
	}
	if(y>y2) {
		y2=y;
	}
	threedproj(-1.,1.,-1.,&x,&y);
	if(x<x1) {
		x1=x;
	}
	if(x>x2) {
		x2=x;
	}
	if(y<y1) {
		y1=y;
	}
	if(y>y2) {
		y2=y;
	}
	threedproj(-1.,1.,1.,&x,&y);
	if(x<x1) {
		x1=x;
	}
	if(x>x2) {
		x2=x;
	}
	if(y<y1) {
		y1=y;
	}
	if(y>y2) {
		y2=y;
	}
	threedproj(1.,-1.,-1.,&x,&y);
	if(x<x1) {
		x1=x;
	}
	if(x>x2) {
		x2=x;
	}
	if(y<y1) {
		y1=y;
	}
	if(y>y2) {
		y2=y;
	}
	threedproj(1.,-1.,1.,&x,&y);
	if(x<x1) {
		x1=x;
	}
	if(x>x2) {
		x2=x;
	}
	if(y<y1) {
		y1=y;
	}
	if(y>y2) {
		y2=y;
	}
	threedproj(1.,1.,1.,&x,&y);
	if(x<x1) {
		x1=x;
	}
	if(x>x2) {
		x2=x;
	}
	if(y<y1) {
		y1=y;
	}
	if(y>y2) {
		y2=y;
	}
	threedproj(1.,1.,-1.,&x,&y);
	if(x<x1) {
		x1=x;
	}
	if(x>x2) {
		x2=x;
	}
	if(y<y1) {
		y1=y;
	}
	if(y>y2) {
		y2=y;
	}
	*xlo=x1;
	*ylo=y1;
	*xhi=x2;
	*yhi=y2;
}


static int create_crv(int ind) {
	int i,type,j;
	int ix,iy,iz;

	for(i=0;i<MAXFRZ;i++) {
		if(frz[i].use==0) {
			ix=MyGraph->xv[ind];
			iy=MyGraph->yv[ind];
			iz=MyGraph->zv[ind];
			if(my_browser.maxrow<=2) {
				err_msg("No Curve to freeze");
				return(-1);
			}
			frz[i].xv=(float *) malloc(sizeof(float)*my_browser.maxrow);
			frz[i].yv=(float *) malloc(sizeof(float)*my_browser.maxrow);
			if((type=MyGraph->grtype)>0) {
				frz[i].zv=(float *)malloc(sizeof(float)*my_browser.maxrow);
			}
			if((type>0 && frz[i].zv==NULL) || (type==0 && frz[i].yv==NULL)) {
				err_msg("Cant allocate storage for curve");
				return(-1);
			}
			frz[i].use=1;
			frz[i].len=my_browser.maxrow;
			for(j=0;j<my_browser.maxrow;j++) {
				frz[i].xv[j]=my_browser.data[ix][j];
				frz[i].yv[j]=my_browser.data[iy][j];
				if(type>0) {
					frz[i].zv[j]=my_browser.data[iz][j];
				}
			}
			frz[i].type=type;
			frz[i].w=draw_win;
			sprintf(frz[i].name,"crv%c",'a'+i);
			sprintf(frz[i].key,"%s",frz[i].name);
			return(i);
		}
	}
	err_msg("All curves used");
	return(-1);
}


static void create_ps(void) {
	char filename[XPP_MAX_NAME];
	static char *nn[]={"BW-0/Color-1","Land(0)/Port(1)","Axes fontsize","Font","Linewidth"};
	int status;
	char values[5][MAX_LEN_SBOX];
	sprintf(values[0],"%d",PS_Color);
	sprintf(values[1],"%d",PS_Port);
	sprintf(values[2],"%d",PS_FontSize);
	sprintf(values[3],"%s",PS_FONT);
	sprintf(values[4],"%g",PS_LineWidth);
	status=do_string_box(5,5,1,"Postscript parameters",nn,values,25);
	if(status!=0) {
		PS_Color=atoi(values[0]);
		PS_Port=atoi(values[1]);
		PS_FontSize=atoi(values[2]);
		PS_LineWidth=atof(values[4]);
		sprintf(PS_FONT,"%s",values[3]);
		sprintf(filename,"%s.ps",this_file);
		ping();

		if(!file_selector("Print postscript",filename,"*.ps")) {
			return;
		}
		if(ps_init(filename,PS_Color)) {
			ps_restore();
			ping();
		}
	}
}


static void create_svg(void) {
	char filename[XPP_MAX_NAME];
	strcpy(filename,this_file);
	filename[strlen(filename)-4]='\0';
	strcat(filename,".svg");
	/*sprintf(filename,"%s.svg",tmp);*/
	if(!file_selector("Print svg",filename,"*.svg")) {
		return;
	}
	if(svg_init(filename,PS_Color)) {
		svg_restore();
		ping();
	}
}


static void delete_frz(void) {
	int i;
	i=get_frz_index(draw_win);
	if(i<0) {
		return;
	}
	delete_frz_crv(i);
}


static void delete_frz_crv(int i) {
	if(frz[i].use==0) {
		return;
	}
	frz[i].use=0;
	frz[i].name[0]=0;
	frz[i].key[0]=0;
	free(frz[i].xv);
	free(frz[i].yv);
	if(frz[i].type>0) {
		free(frz[i].zv);
	}
}


static void draw_bd(Window w) {
	int i,j,len;
	float oldxpl,oldypl,xpl,ypl,*x,*y;
	if(w==my_bd.w && my_bd.nbifcrv>0) {
		for(i=0;i<my_bd.nbifcrv;i++) {
			set_linestyle(my_bd.color[i]);
			len=my_bd.npts[i];
			x=my_bd.x[i];
			y=my_bd.y[i];
			xpl=x[0];
			ypl=y[0];
			for(j=0;j<len;j++) {
				oldxpl=xpl;
				oldypl=ypl;
				xpl=x[j];
				ypl=y[j];
				line_abs(oldxpl,oldypl,xpl,ypl);
			}
		}
	}
}


static void draw_freeze_key(void) {
	int ix,iy;
	int i,y0;
	int ix2;
	int dy=2*HChar;
	if(FreezeKeyFlag==SCRNFMT) {
		return;
	}
	if(PltFmtFlag==PSFMT) {
		dy=-dy;
	}
	scale_to_screen((float)FreezeKeyX,(float)FreezeKeyY,&ix,&iy);
	ix2=ix+4*HChar;
	y0=iy;
	for(i=0;i<MAXFRZ;i++) {
		if(frz[i].use==1 && frz[i].w==draw_win && strlen(frz[i].key)>0) {
			set_linestyle(abs(frz[i].color));
			line(ix,y0,ix2,y0);
			set_linestyle(0);
			put_text(ix2+HChar,y0,frz[i].key);
			y0+=dy;
		}
	}
}


static void draw_frozen_cline(int index, Window w) {
	if(nclines[index].use==0 || nclines[index].w!=w) {
		return;
	}
}


static void edit_curve(void) {
	char bob[20];
	int crv=0;
	sprintf(bob,"Edit 0-%d :",MyGraph->nvars-1);
	ping();
	new_int(bob,&crv);
	if(crv>=0 && crv<MyGraph->nvars) {
		sprintf(bob,"Edit curve %d",crv);
		alter_curve(bob,crv,crv);
	}
}


static void edit_frz(void) {
	int i;
	i=get_frz_index(draw_win);
	if(i<0) {
		return;
	}
	edit_frz_crv(i);
}


static void edit_frz_crv(int i) {
	static char *nn[]={"*4Color","Key","Name"};
	char values[3][MAX_LEN_SBOX];
	int status;
	sprintf(values[0],"%d",frz[i].color);
	sprintf(values[1],"%s",frz[i].key);
	sprintf(values[2],"%s",frz[i].name);
	status=do_string_box(3,3,1,"Edit Freeze",nn,values,25);
	if(status!=0) {
		frz[i].color=atoi(values[0]);
		sprintf(frz[i].key,"%s",values[1]);
		sprintf(frz[i].name,"%s",values[2]);
	}
}


static void export_graf_data(void) {
	FILE *fp;
	/*char filename[DEFAULT_STRING_LENGTH];*/
	char filename[XPP_MAX_NAME];
	sprintf(filename,"curve.dat");
	ping();
	if(!file_selector("Export graph data",filename,"*.dat")) {
		return;
	}
	if((fp=fopen(filename,"w"))==NULL) {
		err_msg("Couldn't open file");
		return;
	}
	export_data(fp);
	fclose(fp);
}


static void fit_window(void) {
	double Mx=-1.e25,My=-1.e25,Mz=-1.e25,mx=-Mx,my=-My,mz=-Mz;
	int i,n=MyGraph->nvars;
	if(storind<2) {
		return;
	}
	if(MyGraph->ThreeDFlag) {
		for(i=0;i<n;i++) {
			get_max(MyGraph->xv[i],&(MyGraph->xmin),&(MyGraph->xmax));
			Mx=lmax(MyGraph->xmax,Mx);
			mx=-lmax(-MyGraph->xmin,-mx);

			get_max(MyGraph->yv[i],&(MyGraph->ymin),&(MyGraph->ymax));
			My=lmax(MyGraph->ymax,My);
			my=-lmax(-MyGraph->ymin,-my);

			get_max(MyGraph->zv[i],&(MyGraph->zmin),&(MyGraph->zmax));
			Mz=lmax(MyGraph->zmax,Mz);
			mz=-lmax(-MyGraph->zmin,-mz);
		}
		MyGraph->xmax=Mx;
		MyGraph->ymax=My;
		MyGraph->zmax=Mz;
		MyGraph->xmin=mx;
		MyGraph->ymin=my;
		MyGraph->zmin=mz;

		corner_cube(&(MyGraph->xlo),&(MyGraph->xhi),&(MyGraph->ylo),&(MyGraph->yhi));
		check_windows();
	} else {
		for(i=0;i<n;i++) {
			get_max(MyGraph->xv[i],&(MyGraph->xmin),&(MyGraph->xmax));
			Mx=lmax(MyGraph->xmax,Mx);
			mx=-lmax(-MyGraph->xmin,-mx);

			get_max(MyGraph->yv[i],&(MyGraph->ymin),&(MyGraph->ymax));
			My=lmax(MyGraph->ymax,My);
			my=-lmax(-MyGraph->ymin,-my);
		}
		MyGraph->xmax=Mx;
		MyGraph->ymax=My;

		MyGraph->xmin=mx;
		MyGraph->ymin=my;

		MyGraph->xlo=MyGraph->xmin;
		MyGraph->ylo=MyGraph->ymin;
		MyGraph->xhi=MyGraph->xmax;
		MyGraph->yhi=MyGraph->ymax;
		check_windows();
	}
	redraw_the_graph();
}


static void free_bd(void) {
	int i;
	if(my_bd.nbifcrv>0) {
		for(i=0;i<my_bd.nbifcrv;i++) {
			free(my_bd.x[i]);
			free(my_bd.y[i]);
		}
		my_bd.nbifcrv=0;
	}
}


static int freeze_crv(int ind) {
	int i;
	i=create_crv(ind);
	if(i<0) {
		return(-1);
	}
	edit_frz_crv(i);
	return(1);
}


static void frz_bd(void) {
	FILE *fp;
	/*char filename[DEFAULT_STRING_LENGTH];*/
	char filename[XPP_MAX_NAME];
	sprintf(filename,"diagram.dat");
	ping();
	if(!file_selector("Import Diagram",filename,"*.dat")) {
		return;
	}
	if((fp=fopen(filename,"r"))==NULL) {
		err_msg("Couldn't open file");
		return;
	}
	read_bd(fp);
}


static void get_2d_view(int ind) {
	static char *n[]={"*0X-axis","*0Y-axis","Xmin", "Ymin",
					  "Xmax", "Ymax", "Xlabel","Ylabel"};
	char values[8][MAX_LEN_SBOX];
	int  status,i;
	int i1=MyGraph->xv[ind],i2=MyGraph->yv[ind];
	char n1[15],n2[15];
	ind_to_sym(i1,n1);
	ind_to_sym(i2,n2);
	sprintf(values[0],"%s",n1);
	sprintf(values[1],"%s",n2);
	sprintf(values[2],"%g",MyGraph->xmin);
	sprintf(values[3],"%g",MyGraph->ymin);
	sprintf(values[4],"%g",MyGraph->xmax);
	sprintf(values[5],"%g",MyGraph->ymax);
	sprintf(values[6],"%s",MyGraph->xlabel);
	sprintf(values[7],"%s",MyGraph->ylabel);
	MyGraph->ThreeDFlag=0;
	status=do_string_box(8,4,2,"2D View",n,values,31);
	if(status!=0) {
		/*  get variable names  */
		find_variable(values[0],&i);
		if(i>-1) {
			MyGraph->xv[ind]=i;
		}
		find_variable(values[1],&i);
		if(i>-1) {
			MyGraph->yv[ind]=i;
		}
		MyGraph->xmin=atof(values[2]);
		MyGraph->ymin=atof(values[3]);
		MyGraph->xmax=atof(values[4]);
		MyGraph->ymax=atof(values[5]);
		MyGraph->xlo=MyGraph->xmin;
		MyGraph->ylo=MyGraph->ymin;
		MyGraph->xhi=MyGraph->xmax;
		MyGraph->yhi=MyGraph->ymax;
		sprintf(MyGraph->xlabel,"%s",values[6]);
		sprintf(MyGraph->ylabel,"%s",values[7]);
		check_windows();
	}
}


static void get_3d_view(int ind) {
	static char *n[]={"*0X-axis","*0Y-axis", "*0Z-axis",
					  "Xmin", "Xmax", "Ymin",
					  "Ymax", "Zmin","Zmax",
					  "XLo", "XHi", "YLo", "YHi","Xlabel","Ylabel","Zlabel"};
	char values[16][MAX_LEN_SBOX];
	int  status,i,i1=MyGraph->xv[ind],i2=MyGraph->yv[ind],i3=MyGraph->zv[ind];
	char n1[15],n2[15],n3[15];
	ind_to_sym(i1,n1);
	ind_to_sym(i2,n2);
	ind_to_sym(i3,n3);
	sprintf(values[0],"%s",n1);
	sprintf(values[1],"%s",n2);
	sprintf(values[2],"%s",n3);
	sprintf(values[3],"%g",MyGraph->xmin);
	sprintf(values[5],"%g",MyGraph->ymin);
	sprintf(values[7],"%g",MyGraph->zmin);
	sprintf(values[4],"%g",MyGraph->xmax);
	sprintf(values[6],"%g",MyGraph->ymax);
	sprintf(values[8],"%g",MyGraph->zmax);
	sprintf(values[9],"%g",MyGraph->xlo);
	sprintf(values[11],"%g",MyGraph->ylo);
	sprintf(values[10],"%g",MyGraph->xhi);
	sprintf(values[12],"%g",MyGraph->yhi);
	sprintf(values[13],"%s",MyGraph->xlabel);
	sprintf(values[14],"%s",MyGraph->ylabel);
	sprintf(values[15],"%s",MyGraph->zlabel);
	MyGraph->ThreeDFlag=1;
	status=do_string_box(16,6,3,"3D View",n,values,31);
	if(status!=0) {
		/*  get variable names  */
		find_variable(values[0],&i);
		if(i>-1) {
			MyGraph->xv[ind]=i;
		}
		find_variable(values[1],&i);
		if(i>-1) {
			MyGraph->yv[ind]=i;
		}
		find_variable(values[2],&i);
		if(i>-1) {
			MyGraph->zv[ind]=i;
		}
		sprintf(MyGraph->xlabel,"%s",values[13]);
		sprintf(MyGraph->ylabel,"%s",values[14]);
		sprintf(MyGraph->zlabel,"%s",values[15]);


		MyGraph->xmin=atof(values[3]);
		MyGraph->ymin=atof(values[5]);
		MyGraph->zmin=atof(values[7]);
		MyGraph->xmax=atof(values[4]);
		MyGraph->ymax=atof(values[6]);
		MyGraph->zmax=atof(values[8]);
		MyGraph->xlo=atof(values[9]);
		MyGraph->ylo=atof(values[11]);
		MyGraph->xhi=atof(values[10]);
		MyGraph->yhi=atof(values[12]);
		check_windows();
	}
}


static int get_frz_index(Window w) {
	char *n[MAXFRZ];
	char key[MAXFRZ],ch;

	int i;
	int count=0;
	Window temp=main_win;
	for(i=0;i<MAXFRZ;i++) {
		if(frz[i].use==1 && w==frz[i].w) {
			n[count]=(char *)malloc(20);
			sprintf(n[count],"%s",frz[i].name);
			key[count]='a'+i;

			count++;
		}
	}
	if(count==0) {
		return(-1);
	}
	key[count]=0;
	ch=(char)pop_up_list(&temp,"Curves",n,key,count,12,0,10,8*DCURY+8,
						 no_hint,info_pop,info_message);
	for(i=0;i<count;i++) {
		free(n[i]);
	}
	return((int)(ch-'a'));
}


static void kill_frz(void) {
	int i;
	for(i=0;i<MAXFRZ;i++) {
		if(frz[i].use==1 && frz[i].w==draw_win) {
			delete_frz_crv(i);
		}
	}
}


static void movie_rot(double start, double increment, int nclip, int angle) {
	int i;
	double thetaold=MyGraph->Theta,phiold=MyGraph->Phi;
	reset_film();
	for(i=0;i<=nclip;i++) {
		if(angle==0) {
			make_rot(start+i*increment,phiold);
		} else {
			make_rot(thetaold,start+i*increment);
		}
		redraw_the_graph();
		film_clip();
	}
	MyGraph->Theta=thetaold;
	MyGraph->Phi=phiold;
}


static void new_curve(void) {
	if(alter_curve("New Curve",0,MyGraph->nvars)) {
		MyGraph->nvars=MyGraph->nvars+1;
	}
}


static void read_bd(FILE *fp) {
	int oldtype,type,oldbr,br,ncrv=0,len,f2;
	float x[8000],ylo[8000],yhi[8000];
	len=0;
	if(fscanf(fp,"%g %g %g %d %d %d",&x[len],&ylo[len],&yhi[len],&oldtype,&oldbr,&f2)<6) {
		plintf("Error reading variables from file %s", fp);
	}
	len++;
	while(!feof(fp)) {
		if(fscanf(fp,"%g %g %g %d %d %d",&x[len],&ylo[len],&yhi[len],&type,&br,&f2)<6) {
			plintf("Error reading variables from file %s", fp);
		}
		if(type==oldtype && br==oldbr) {
			len++;
		} else {
			add_bd_crv(x,ylo,len,oldtype,ncrv);
			ncrv++;
			if(oldtype==UNSTABLE_PERIODIC || oldtype==STABLE_PERIODIC) {
				add_bd_crv(x,yhi,len,oldtype,ncrv);
				ncrv++;
			}
			if(oldbr==br) {
				len--;
			}
			x[0]=x[len];
			ylo[0]=ylo[len];
			yhi[0]=yhi[len];

			len=1;
		}
		oldbr=br;
		oldtype=type;
	}
	/*  save this last one */
	if(len>1) {
		add_bd_crv(x,ylo,len,oldtype,ncrv);
		ncrv++;
		if(oldtype==UNSTABLE_PERIODIC || oldtype==STABLE_PERIODIC) {
			add_bd_crv(x,yhi,len,oldtype,ncrv);
			ncrv++;
		}
	}
	plintf( " got %d bifurcation curves\n",ncrv);
	fclose(fp);
	my_bd.nbifcrv=ncrv;
	my_bd.w=draw_win;
}


static void scroll_window(void) {
	XEvent ev;
	int i=0,j=0;
	int state=0;
	float x,y,x0,y0;
	float xlo=MyGraph->xlo;
	float ylo=MyGraph->ylo;
	float xhi=MyGraph->xhi;
	float yhi=MyGraph->yhi;
	float dx=0.0,dy=0.0;
	int alldone=0;
	XSelectInput(display,draw_win,
				 KeyPressMask|ButtonPressMask|ButtonReleaseMask|
				 PointerMotionMask|ButtonMotionMask|ExposureMask);
	while(!alldone) {
		XNextEvent(display,&ev);
		switch(ev.type) {
		case KeyPress:
			alldone=1;
			break;
		case Expose:
			do_expose(ev);
			break;
		case ButtonPress:
			if(state==0) {
				i=ev.xkey.x;
				j=ev.xkey.y;
				scale_to_real(i,j,&x0,&y0);
				state=1;
			}
			break;
		case MotionNotify:
			if(state==1) {
				i=ev.xmotion.x;
				j=ev.xmotion.y;
				scale_to_real(i,j,&x,&y);
				dx=-(x-x0)/2;
				dy=-(y-y0)/2;
				update_view(xlo+dx,xhi+dx,ylo+dy,yhi+dy);
			}
			break;
		case ButtonRelease:
			state=0;
			xlo=xlo+dx;
			xhi=xhi+dx;
			ylo=ylo+dy;
			yhi=yhi+dy;
			break;
		}
	}
}


static void set_key(int x, int y) {
	float xp,yp;
	scale_to_real(x,y,&xp,&yp);
	FreezeKeyX=xp;
	FreezeKeyY=yp;
	FreezeKeyFlag=1;
}


static void update_view(float xlo,float xhi, float ylo, float yhi) {
	MyGraph->xlo=xlo;
	MyGraph->ylo=ylo;
	MyGraph->xhi=xhi;
	MyGraph->yhi=yhi;
	if(MyGraph->grtype<5) {
		MyGraph->xmin=MyGraph->xlo;
		MyGraph->xmax=MyGraph->xhi;
		MyGraph->ymin=MyGraph->ylo;
		MyGraph->ymax=MyGraph->yhi;
	}
	check_windows();
	redraw_the_graph();
}


static void user_window(void) {
	static char *n[]={"X Lo","X Hi","Y Lo","Y Hi"};
	char values[4][MAX_LEN_SBOX];
	int status;
	sprintf(values[0],"%g",MyGraph->xlo);
	sprintf(values[2],"%g",MyGraph->ylo);
	sprintf(values[1],"%g",MyGraph->xhi);
	sprintf(values[3],"%g",MyGraph->yhi);
	status=do_string_box(4,2,2,"Window",n,values,28);
	if(status!=0) {
		MyGraph->xlo=atof(values[0]);
		MyGraph->ylo=atof(values[2]);
		MyGraph->xhi=atof(values[1]);
		MyGraph->yhi=atof(values[3]);
		if(MyGraph->grtype<5) {
			MyGraph->xmin=MyGraph->xlo;
			MyGraph->xmax=MyGraph->xhi;
			MyGraph->ymin=MyGraph->ylo;
			MyGraph->ymax=MyGraph->yhi;
		}
		check_windows();
	}
	redraw_the_graph();
}


static void zoom_in(int i1, int j1, int i2, int j2) {
	float x1,y1,x2,y2;
	float dx=MyGraph->xhi-MyGraph->xlo;
	float dy=MyGraph->yhi-MyGraph->ylo;
	scale_to_real(i1,j1,&x1,&y1);
	scale_to_real(i2,j2,&x2,&y2);
	if(x1==x2 || y1==y2) {
		if(dx < 0) {
			dx=-dx;
		}
		if(dy < 0) {
			dy=-dy;
		}
		dx = dx/2;
		dy = dy/2;

		/*Shrink by thirds and center (track) about the point clicked*/
		MyGraph->xlo=x1-dx/2;
		MyGraph->xhi=x1+dx/2;

		MyGraph->ylo=y1-dy/2;
		MyGraph->yhi=y1+dy/2;
	} else {
		MyGraph->xlo=x1;
		MyGraph->ylo=y1;
		MyGraph->xhi=x2;
		MyGraph->yhi=y2;
	}
	if(MyGraph->grtype<5) {
		MyGraph->xmin=MyGraph->xlo;
		MyGraph->xmax=MyGraph->xhi;
		MyGraph->ymin=MyGraph->ylo;
		MyGraph->ymax=MyGraph->yhi;
	}
	check_windows();
	redraw_the_graph();
	draw_help();
}


static void zoom_out(int i1, int j1, int i2, int j2) {
	float x1,y1,x2,y2;
	float bx,mux,by,muy;
	float dx=MyGraph->xhi-MyGraph->xlo;
	float dy=MyGraph->yhi-MyGraph->ylo;
	scale_to_real(i1,j1,&x1,&y1);
	scale_to_real(i2,j2,&x2,&y2);

	if(x1==x2 || y1==y2) {
		if(dx < 0) {
			dx=-dx;
		}
		if(dy < 0) {
			dy=-dy;
		}
		/*Grow by thirds and center (track) about the point clicked*/
		dx = dx*2;
		dy = dy*2;

		MyGraph->xlo=x1-dx/2;
		MyGraph->xhi=x1+dx/2;

		MyGraph->ylo=y1-dy/2;
		MyGraph->yhi=y1+dy/2;
	} else {
		if(x1>x2) {
			bx=x1;
			x1=x2;
			x2=bx;
		}
		if(y1>y2) {
			by=y1;
			y1=y2;
			y2=by;
		}

		bx=dx*dx/(x2-x1);
		mux=(x1-MyGraph->xlo)/dx;
		MyGraph->xlo=MyGraph->xlo-bx*mux;
		MyGraph->xhi=MyGraph->xlo+bx;

		by=dy*dy/(y2-y1);
		muy=(y1-MyGraph->ylo)/dy;
		MyGraph->ylo=MyGraph->ylo-by*muy;
		MyGraph->yhi=MyGraph->ylo+by;
	}
	if(MyGraph->grtype<5) {
		MyGraph->xmin=MyGraph->xlo;
		MyGraph->xmax=MyGraph->xhi;
		MyGraph->ymin=MyGraph->ylo;
		MyGraph->ymax=MyGraph->yhi;
	}
	check_windows();
	redraw_the_graph();
	draw_help();
}
