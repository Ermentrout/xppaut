#include "nullcline.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "abort.h"
#include "browse.h"
#include "form_ode.h"
#include "ggets.h"
#include "graf_par.h"
#include "graphics.h"
#include "init_conds.h"
#include "integrate.h"
#include "load_eqn.h"
#include "main.h"
#include "many_pops.h"
#include "markov.h"
#include "menudrive.h"
#include "my_rhs.h"
#include "my_ps.h"
#include "my_svg.h"
#include "numerics.h"
#include "parserslow.h"
#include "pop_list.h"
#include "struct.h"
#include "util/timeutil.h"
#include "xpplim.h"

/* --- Macros --- */
#define MAX_NULL 10000


/* --- Types --- */
typedef struct {
	float x,y,z;
} Pt;


typedef struct nclines {
	float *xn,*yn;
	int nmx,nmy;
	int n_ix,n_iy;
	struct nclines *n,*p;
}  NCLINES;


/* --- Forward declarations --- */
static void add_froz_cline(float *xn, int nmx, int n_ix, float *yn, int nmy, int n_iy);
static void clear_froz_cline(void);
static void do_cline(int ngrid, double x1, double y1, double x2, double y2);
static void dump_clines(FILE *fp, float *x, int nx, float *y, int ny);
static float fnull(double x, double y);
static void get_max_dfield(double *y, double *ydot, double u0, double v0, double du, double dv, int n, int inx, int iny, double *mdf);
static int interpolate(Pt p1, Pt p2, double z, float *x, float *y);
static void new_nullcline(int course, double xlo, double ylo, double xhi, double yhi, float *stor, int *npts);
static void quad_contour(Pt p1, Pt p2, Pt p3, Pt p4);
static void redraw_froz_cline(int flag);
static void restor_null(float *v, int n, int d);
static void save_frozen_clines(char *fn);
static void save_the_nullclines(void);
static void start_ncline(void);
static void stor_null(double x1, double y1, double x2, double y2);


/* --- Data --- */
static int DF_IX=-1,DF_IY=-1;
static int DFIELD_TYPE=0;
static int DFSuppress=0;
static int NCSuppress=0;
static int n_nstore=0;
static int ncline_cnt;
static int NullStyle=0; /* 1 is with little vertical/horizontal lines */
static int null_ix,null_iy,WHICH_CRV;
static int num_x_n,num_y_n,num_index;
static float *X_n,*Y_n,*saver,*NTop,*NBot;

static NCLINES *ncperm;
static RANGE_INFO ncrange;

int DOING_DFIELD=0;
int DFBatch=0;
int DF_GRID=16,DF_FLAG=0;
int NCBatch=0;
int NULL_HERE;
int XNullColor=2;
int YNullColor=7;


/* --- Functions --- */
void create_new_cline(void) {
	if(NULL_HERE) {
		new_clines_com(0);
	}
}


void direct_field_com(int c) {
	int i,j,start,k;
	int inx=MyGraph->xv[0]-1;
	int iny=MyGraph->yv[0]-1;
	double y[MAXODE],ydot[MAXODE],xv1,xv2;
	double dtold=DELTA_T;
	float v1[MAXODE],v2[MAXODE];

	double amp,mdf;
	double t;
	double du,dv,u0,v0,dxp,dyp,dz,dup,dvp;
	double oldtrans=TRANS;
	int grid=DF_GRID;

	if(MyGraph->TimeFlag				||
	   MyGraph->xv[0]==MyGraph->yv[0] ||
	   MyGraph->ThreeDFlag) {
		return;
	}

	if(c==2) {
		DF_FLAG=0;
		return;
	}
	if(c==0) {
		DFIELD_TYPE=1;
	}
	if(c==4) {
		DFIELD_TYPE=0;
	}
	new_int("Grid:",&grid);
	if(grid<=1) {
		return;
	}
	DF_GRID=grid;
	du=(MyGraph->xhi-MyGraph->xlo)/(double)grid;
	dv=(MyGraph->yhi-MyGraph->ylo)/(double)grid;

	dup =(double)(DRight-DLeft)/(double)grid;
	dvp=(double)(DTop-DBottom)/(double)grid;
	dz=hypot(dup,dvp)*(.25+.75*DFIELD_TYPE) ;
	u0=MyGraph->xlo;
	v0=MyGraph->ylo;
	set_linestyle(MyGraph->color[0]);
	if(c!=1) {
		DF_FLAG=1;
		if(c==3) {
			DF_FLAG=2;
			du=(MyGraph->xhi-MyGraph->xlo)/(double)(grid+1);
			dv=(MyGraph->yhi-MyGraph->ylo)/(double)(grid+1);
		}
		DF_IX=inx+1;
		DF_IY=iny+1;
		get_ic(2,y);
		get_max_dfield(y,ydot,u0,v0,du,dv,grid,inx,iny,&mdf);
		if (PltFmtFlag==SVGFMT) {
			DOING_DFIELD=1;
			fprintf(svgfile,"<g>\n");
		}
		for(i=0;i<=grid;i++) {
			y[inx]=u0+du*i;
			for(j=0;j<=grid;j++) {
				y[iny]=v0+dv*j;
				my_rhs(0.0,y,ydot,NODE);
				extra(y,0.0,NODE,NEQ);
				if(MyGraph->ColorFlag || DF_FLAG==2) {
					v1[0]=0.0;
					v2[0]=0.0;
					for(k=0;k<NEQ;k++) {
						v1[k+1]=(float)y[k];
						v2[k+1]=v1[k+1]+(float)ydot[k];
					}
					comp_color(v1,v2,NODE,1.0);
				}
				if(DF_FLAG==1) {
					scale_dxdy(ydot[inx],ydot[iny],&dxp,&dyp);
					if(DFIELD_TYPE==0) {
						amp=hypot(dxp,dyp);
						if(amp!=0.0) {
							ydot[inx]/=amp;
							ydot[iny]/=amp;
						}
					} else {
						ydot[inx]/=mdf;
						ydot[iny]/=mdf;
					}
					xv1=y[inx]+ydot[inx]*dz;
					xv2=y[iny]+ydot[iny]*dz;
					bead_abs((float)xv1,(float)xv2);
					line_abs((float)y[inx],(float)y[iny],(float)xv1,(float)xv2);
				}
				if(DF_FLAG==2 && j>0 && i<grid) {
					frect_abs((float)y[inx],(float)y[iny],(float)du,(float)dv);
				}
			}
		}
		TRANS=oldtrans;
		if (PltFmtFlag==SVGFMT) {
			DOING_DFIELD=0;
			fprintf(svgfile,"</g>\n");
		}
		return;
	}
	STORFLAG=0;

	SuppressBounds=1;
	for(k=0;k<2;k++) {
		for(i=0;i<=grid;i++) {
			for(j=0;j<=grid;j++) {
				get_ic(2,y);
				y[inx]=u0+du*i;
				y[iny]=v0+dv*j;
				t=0.0;
				start=1;
				integrate(&t,y,TEND,DELTA_T,1,NJMP,&start);
			}
		}
		DELTA_T=-DELTA_T;
	}
	SuppressBounds=0;
	DELTA_T=dtold;
	if (PltFmtFlag==SVGFMT) {
		DOING_DFIELD=0;
		fprintf(svgfile,"</g>\n");
	}
}


/*  all the nifty 2D stuff here    */
void do_batch_dfield(void) {
	if(!XPPBatch) {
		return;
	}
	switch(DFBatch) {
	case 0:
		return;
	case 1:
		DF_FLAG=1;
		DFIELD_TYPE=1;
		DF_IX=MyGraph->xv[0];
		DF_IY=MyGraph->yv[0];
		redraw_dfield();
		return;
	case 2:
		DF_FLAG=1;
		DFIELD_TYPE=0;
		DF_IX=MyGraph->xv[0];
		DF_IY=MyGraph->yv[0];
		redraw_dfield();
		return;
	case 3:
		DF_FLAG=2;
		DFIELD_TYPE=0;
		DF_IX=MyGraph->xv[0];
		DF_IY=MyGraph->yv[0];
		redraw_dfield();
		return;
	case 4:
		DF_FLAG=1;
		DFIELD_TYPE=1;
		DF_IX=MyGraph->xv[0];
		DF_IY=MyGraph->yv[0];
		redraw_dfield();
		return;
	case 5:
		DF_FLAG=1;
		DFIELD_TYPE=0;
		DF_IX=MyGraph->xv[0];
		DF_IY=MyGraph->yv[0];
		redraw_dfield();
		return;
	}
}


void do_batch_nclines(void) {
	if(!XPPBatch || !NCBatch) {
		return;
	}
	if(NCBatch==1) {
		new_clines_com(0);
		return;
	}

}


void do_range_clines(void) {
	static char *n[]={"*2Range parameter","Steps","Low","High"};
	char values[4][MAX_LEN_SBOX];
	int status,i;
	double z,dz,zold;
	float xmin,xmax,y_tp,y_bot;
	int col1=XNullColor,col2=YNullColor;
	int course=NMESH;
	sprintf(values[0],"%s",ncrange.rv);
	sprintf(values[1],"%d",ncrange.nstep);
	sprintf(values[2],"%g",ncrange.xlo);
	sprintf(values[3],"%g",ncrange.xhi);
	status=do_string_box(4,4,1,"Range Clines",n,values,45);
	if(status!=0) {
		strcpy(ncrange.rv,values[0]);
		ncrange.nstep=atoi(values[1]);
		ncrange.xlo=atof(values[2]);
		ncrange.xhi=atof(values[3]);
		if(ncrange.nstep<=0) {
			return;
		}
		dz=(ncrange.xhi-ncrange.xlo)/(double)ncrange.nstep;
		if(dz<=0.0) {
			return;
		}
		get_val(ncrange.rv,&zold);

		for(i=NODE;i<NODE+NMarkov;i++) {
			set_ivar(i+1+FIX_VAR,last_ic[i]);
		}
		xmin=(float)MyGraph->xmin;
		xmax=(float)MyGraph->xmax;
		y_tp=(float)MyGraph->ymax;
		y_bot=(float)MyGraph->ymin;
		null_ix=MyGraph->xv[0];
		null_iy=MyGraph->yv[0];

		for(i=0;i<=ncrange.nstep;i++) {
			z=(double)i*dz+ncrange.xlo;
			set_val(ncrange.rv,z);
			if(NULL_HERE==0) {
				X_n=(float *)malloc(4*MAX_NULL*sizeof(float));
				Y_n=(float *)malloc(4*MAX_NULL*sizeof(float));
				if(X_n!=NULL  &&  Y_n!=NULL) {
					NULL_HERE=1;
				}
				NTop=(float *)malloc((course+1)*sizeof(float));
				NBot=(float *)malloc((course+1)*sizeof(float));
				if(NTop==NULL || NBot==NULL) {
					NULL_HERE=0;
				}
			} else {
				free(NTop);
				free(NBot);
				NTop=(float *)malloc((course+1)*sizeof(float));
				NBot=(float *)malloc((course+1)*sizeof(float));
				if(NTop==NULL || NBot==NULL) {
					NULL_HERE=0;
					return;
				}
			}
			WHICH_CRV=null_ix;
			set_linestyle(col1);
			new_nullcline(course,xmin,y_bot,xmax,y_tp,X_n,&num_x_n);

			WHICH_CRV=null_iy;
			set_linestyle(col2);
			new_nullcline(course,xmin,y_bot,xmax,y_tp,Y_n,&num_y_n);
			add_froz_cline(X_n,num_x_n,null_ix,Y_n,num_y_n,null_iy);
		}
		set_val(ncrange.rv,zold);
	}
}


void froz_cline_stuff_com(int i) {
	int delay=200;
	if(n_nstore==0) {
		start_ncline();
	}
	switch(i) {
	case 0:
		if(NULL_HERE==0) {
			return;
		}
		add_froz_cline(X_n,num_x_n,null_ix,Y_n,num_y_n,null_iy);
		break;
	case 1:
		clear_froz_cline();
		break;
	case 3:
		new_int("Delay (msec)",&delay);
		if(delay<=0) {
			delay=0;
		}
		redraw_froz_cline(delay);
		break;
	case 2:
		do_range_clines();
		break;
	}
}


/* type=0,1 */
int get_nullcline_floats(float **v,int *n,int who,int type) {
	NCLINES *z;
	int i;
	if(who<0) {
		if(type==0) {
			*v=X_n;
			*n=num_x_n;
		} else {
			*v=Y_n;
			*n=num_y_n;
		}
		if(v==NULL) {
			return 1;
		}
		return 0;
	}
	if(who>ncline_cnt || n_nstore==0) {
		return 1;
	}
	z=ncperm;
	for(i=0;i<who;i++) {
		z=z->n;
	}
	if(z==NULL) {
		return 1;
	}
	if(type==0) {
		*v=z->xn;
		*n=z->nmx;
	} else {
		*v=z->yn;
		*n=z->nmy;
	}
	if(v==NULL) {
		return 1;
	}
	return 0;
}


void new_clines_com(int c) {
	int course=NMESH,i;
	float xmin,xmax,y_tp,y_bot;
	int col1=XNullColor,col2=YNullColor;

	if(MyGraph->ThreeDFlag	||
	   MyGraph->TimeFlag	||
	   MyGraph->xv[0]==MyGraph->yv[0]) {
		return;
	}
	if(c==1) {
		restore_nullclines();
		return;
	}
	if(c==2) {
		MyGraph->Nullrestore=1;
		return;
	}
	if(c==3) {
		MyGraph->Nullrestore=0;
		return;
	}
	if(c==4) {
		froz_cline_stuff();
		return;
	}
	if(c==5) {
		save_the_nullclines();
		return;
	}
	if(c==0) {
		for(i=NODE;i<NODE+NMarkov;i++) {
			set_ivar(i+1+FIX_VAR,last_ic[i]);
		}
		xmin=(float)MyGraph->xmin;
		xmax=(float)MyGraph->xmax;
		y_tp=(float)MyGraph->ymax;
		y_bot=(float)MyGraph->ymin;
		null_ix=MyGraph->xv[0];
		null_iy=MyGraph->yv[0];
		if(NULL_HERE==0) {
			X_n=(float *)malloc(4*MAX_NULL*sizeof(float));
			Y_n=(float *)malloc(4*MAX_NULL*sizeof(float));
			if(X_n!=NULL && Y_n!=NULL) {
				NULL_HERE=1;
			}
			NTop=(float *)malloc((course+1)*sizeof(float));
			NBot=(float *)malloc((course+1)*sizeof(float));
			if(NTop==NULL || NBot==NULL) {
				NULL_HERE=0;
			}
		} else {
			free(NTop);
			free(NBot);
			NTop=(float *)malloc((course+1)*sizeof(float));
			NBot=(float *)malloc((course+1)*sizeof(float));
			if(NTop==NULL || NBot==NULL) {
				NULL_HERE=0;
				return;
			}
		}

		WHICH_CRV=null_ix;
		if(!NCSuppress) {
			set_linestyle(col1);
		}
		new_nullcline(course,xmin,y_bot,xmax,y_tp,X_n,&num_x_n);
		ping();

		WHICH_CRV=null_iy;
		if(!NCSuppress) {
			set_linestyle(col2);
		}
		new_nullcline(course,xmin,y_bot,xmax,y_tp,Y_n,&num_y_n);
		ping();
	}
}


void redraw_dfield(void) {
	int i,j,k;
	int inx=MyGraph->xv[0]-1;
	int iny=MyGraph->yv[0]-1;
	double y[MAXODE],ydot[MAXODE],xv1,xv2;
	float v1[MAXODE],v2[MAXODE];
	FILE *fp;
	double amp,mdf;
	double du,dv,u0,v0,dxp,dyp,dz,dup,dvp;

	int grid=DF_GRID;
	if(DF_FLAG==0 || MyGraph->TimeFlag || MyGraph->xv[0]==MyGraph->yv[0] ||
	   MyGraph->ThreeDFlag || DF_IX!=MyGraph->xv[0] || DF_IY!=MyGraph->yv[0]) {
		return;
	}
	if(DFSuppress==1) {
		fp=fopen("dirfields.dat","w");
		if(fp==NULL) {
			return;
		}
	}

	du=(MyGraph->xhi-MyGraph->xlo)/(double)grid;
	dv=(MyGraph->yhi-MyGraph->ylo)/(double)grid;
	dup =(double)(DRight-DLeft)/(double)grid;
	dvp=(double)(DTop-DBottom)/(double)grid;
	dz=hypot(dup,dvp)*(.25+.75*DFIELD_TYPE);
	u0=MyGraph->xlo;
	v0=MyGraph->ylo;
	if(!DFSuppress) {
		set_linestyle(MyGraph->color[0]);
	}
	get_ic(2,y);
	get_max_dfield(y,ydot,u0,v0,du,dv,grid,inx,iny,&mdf);
	if (PltFmtFlag==SVGFMT) {
		DOING_DFIELD=1;
		fprintf(svgfile,"<g>\n");
	}
	for(i=0;i<=grid;i++) {
		y[inx]=u0+du*i;
		for(j=0;j<=grid;j++) {
			y[iny]=v0+dv*j;
			my_rhs(0.0,y,ydot,NODE);
			extra(y,0.0,NODE,NEQ);
			if(MyGraph->ColorFlag || DF_FLAG==2) {
				v1[0]=0.0;
				v2[0]=0.0;
				for(k=0;k<NEQ;k++) {
					v1[k+1]=(float)y[k];
					v2[k+1]=v1[k+1]+(float)ydot[k];
				}
				if(!DFSuppress) {
					comp_color(v1,v2,NODE,1.0);
				}
			}
			if(DF_FLAG==1 || DF_FLAG==4) {
				scale_dxdy(ydot[inx],ydot[iny],&dxp,&dyp);
				if(DFIELD_TYPE==1) {
					ydot[inx]/=mdf;
					ydot[iny]/=mdf;
				} else {
					amp=hypot(dxp,dyp);
					if(amp!=0.0) {
						ydot[inx]/=amp;
						ydot[iny]/=amp;
					}
				}
				xv1=y[inx]+ydot[inx]*dz;
				xv2=y[iny]+ydot[iny]*dz;
				if(!DFSuppress) {
					bead_abs((float)xv1,(float)xv2);
					line_abs((float)y[inx],(float)y[iny],(float)xv1,(float)xv2);
				} else{
					fprintf(fp,"%g %g %g %g\n",y[inx],y[iny],xv1,xv2);
				}
			}
			if(DF_FLAG==2 && j>0 && i<grid) {
				frect_abs((float)y[inx],(float)y[iny],(float)du,(float)dv);
			}
		}
	}
	if(PltFmtFlag==SVGFMT) {
		DOING_DFIELD=0;
		fprintf(svgfile,"</g>\n");
	}
	if(DFSuppress==1) {
		fclose(fp);
	}
	DFSuppress=0;
}


void restore_nullclines(void) {
	int col1=XNullColor,col2=YNullColor;
	if(NULL_HERE==0) {
		return;
	}
	if(MyGraph->xv[0]==null_ix &&
	   MyGraph->yv[0]==null_iy &&
	   MyGraph->ThreeDFlag==0) {
		set_linestyle(col1);
		restor_null(X_n,num_x_n,1);
		set_linestyle(col2);
		restor_null(Y_n,num_y_n,2);
	}
	redraw_froz_cline(0);
}


void silent_dfields(void) {
	if(DFBatch==5  || DFBatch==4) {
		DFSuppress=1;
		init_ps();
		do_batch_dfield();
		DFSuppress=0;
	}
}


void silent_nullclines(void) {
	FILE *fp;
	if(NCBatch!=2) {
		return;
	}
	NCSuppress=1;
	new_clines_com(0);
	fp=fopen("nullclines.dat","w");
	if(fp==NULL) {
		plintf("Cannot open nullcline file\n");
		return;
	}
	dump_clines(fp,X_n,num_x_n,Y_n,num_y_n);
	fclose(fp);
	NCSuppress=0;
}


/* --- Static functions --- */
static void add_froz_cline(float *xn, int nmx, int n_ix, float *yn, int nmy, int n_iy) {
	NCLINES *z,*znew;
	int i;
	z=ncperm;
	/* move to end */
	while(z->n!=NULL) {
		z=(z->n);
	}
	z->xn=(float *)malloc(4*nmx*sizeof(float));
	for(i=0;i<4*nmx;i++) {
		z->xn[i]=xn[i];
	}
	z->yn=(float *)malloc(4*nmy*sizeof(float));
	for(i=0;i<4*nmy;i++) {
		z->yn[i]=yn[i];
	}
	z->nmx=nmx;
	z->nmy=nmy;
	z->n_ix=n_ix;
	z->n_iy=n_iy;
	z->n=(NCLINES *)malloc(sizeof(NCLINES));
	znew=z->n;
	znew->n=NULL;
	znew->p=z;
	znew->nmx=0;
	znew->nmy=0;
	znew->n_ix=-5;
	znew->n_iy=-5;
	ncline_cnt++;
}


static void clear_froz_cline(void) {
	NCLINES *z,*znew;
	z=ncperm;
	while(z->n!=NULL) {
		z=z->n;
	}
	/*  this is the bottom but there is nothing here that has been stored   */
	znew=z->p;
	if(znew==NULL) {
		return;
	}
	free(z);
	z=znew;
	/* now we are deleting everything */
	while(z->p !=NULL) {
		znew=z->p;
		z->n=NULL;
		z->p=NULL;
		free(z->xn);
		free(z->yn);
		free(z);
		z=znew;
	}
	if(ncperm->nmx>0) {
		free(ncperm->xn);
		ncperm->nmx=0;
	}
	if(ncperm->nmy>0) {
		free(ncperm->yn);
		ncperm->nmy=0;
	}
	ncperm->n=NULL;
	n_nstore=1;
	ncline_cnt=0;
}


static void do_cline(int ngrid, double x1, double y1, double x2, double y2) {
	float dx=(x2-x1)/(float)ngrid;
	float dy=(y2-y1)/(float)ngrid;
	float x,y;
	Pt p[5];
	int i,j;
	int nx=ngrid+1;
	int ny=ngrid+1;

	y=y2;
	for(i=0;i<nx;i++) {
		x=x1+i*dx;
		NBot[i]=fnull(x,y);
	}
	for(j=1;j<ny;j++) {
		y=y2-j*dy;
		NTop[0]=NBot[0];
		NBot[0]=fnull(x1,y);
		for(i=1;i<nx;i++) {
			x=x1+i*dx;
			NTop[i]=NBot[i];
			NBot[i]=fnull(x,y);
			p[0].x=x-dx;
			p[0].y=y+dy;
			p[0].z=NTop[i-1];
			p[1].x=x;
			p[1].y=y+dy;
			p[1].z=NTop[i];
			p[3].x=x-dx;
			p[3].y=y;
			p[3].z=NBot[i-1];
			p[2].x=x;
			p[2].y=y;
			p[2].z=NBot[i];
			/*      Uncomment for triangle contour
	  p[4].x=.25*(p[0].x+p[1].x+p[2].x+p[3].x);
	 p[4].y=.25*(p[0].y+p[1].y+p[2].y+p[3].y);
	 p[4].z=.25*(p[0].z+p[1].z+p[2].z+p[3].z);


	 triangle_contour(p[0],p[1],p[4]);
	 triangle_contour(p[1],p[4],p[2]);
	 triangle_contour(p[4],p[3],p[2]);
	 triangle_contour(p[0],p[4],p[3]); */
			/*   Uncomment for quad contour     */
			quad_contour(p[0],p[1],p[2],p[3]);
			/*     FlushDisplay(); */
		}
	}
}


/* gnuplot format */
static void dump_clines(FILE *fp, float *x, int nx, float *y, int ny) {
	int i;
	fprintf(fp,"# X-nullcline\n");
	for(i=0;i<nx-1;i++) {
		fprintf(fp,"%g %g 1 \n",x[4*i],x[4*i+1]);
		fprintf(fp,"%g %g 1 \n",x[4*i+2],x[4*i+3]);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n# Y-nullcline\n");
	for(i=0;i<ny-1;i++) {
		fprintf(fp,"%g %g 2 \n",y[4*i],y[4*i+1]);
		fprintf(fp,"%g %g 2 \n",y[4*i+2],y[4*i+3]);
		fprintf(fp,"\n");
	}
}


static float fnull(double x, double y) {
	double y1[MAXODE],ydot[MAXODE];
	int i;
	for(i=0;i<NODE;i++) {
		y1[i]=last_ic[i];
	}
	y1[null_ix-1]=(double)x;
	y1[null_iy-1]=(double)y;
	my_rhs(0.0,y1,ydot,NODE);
	return((float)ydot[WHICH_CRV-1]);
}


static void get_max_dfield(double *y, double *ydot, double u0, double v0, double du, double dv,
						   int n, int inx, int iny, double *mdf) {
		   int i,j;
		   double amp,dxp,dyp;
		   *mdf=0.0;
		   for(i=0;i<=n;i++) {
			   y[inx]=u0+du*i;
			   for(j=0;j<=n;j++) {
				   y[iny]=v0+dv*j;
				   my_rhs(0.0,y,ydot,NODE);
				   extra(y,0.0,NODE,NEQ);
				   scale_dxdy(ydot[inx],ydot[iny],&dxp,&dyp);
				   amp=hypot(dxp,dyp);
				   if(amp>*mdf) {
					   *mdf=amp;
				   }
			   }
		   }
	   }


static int interpolate(Pt p1, Pt p2, double z, float *x, float *y) {
	float scale;
	if(p1.z==p2.z) {
		return(0);
	}
	scale=(z-p1.z)/(p2.z-p1.z);
	*x=p1.x+scale*(p2.x-p1.x);
	*y=p1.y+scale*(p2.y-p1.y);
	return(1);
}


static void new_nullcline(int course, double xlo, double ylo, double xhi, double yhi, float *stor, int *npts) {
	num_index=0;
	saver=stor;
	do_cline(course,xlo,ylo,xhi,yhi);
	*npts=num_index;
}


static void quad_contour(Pt p1, Pt p2, Pt p3, Pt p4) {
	float x[4],y[4];
	int count=0;
	if(p1.z*p2.z<=0.0) {
		if(interpolate(p1,p2,0.0,&x[count],&y[count])) {
			count++;
		}
	}
	if(p2.z*p3.z<=0.0) {
		if(interpolate(p3,p2,0.0,&x[count],&y[count])) {
			count++;
		}
	}
	if(p3.z*p4.z<=0.0) {
		if(interpolate(p3,p4,0.0,&x[count],&y[count])) {
			count++;
		}
	}
	if(p1.z*p4.z<=0.0) {
		if(interpolate(p1,p4,0.0,&x[count],&y[count])) {
			count++;
		}
	}
	if(count==2) {
		if(!NCSuppress){
			line_abs(x[0],y[0],x[1],y[1]);
		}
		stor_null(x[0],y[0],x[1],y[1]);
	}
}


static void redraw_froz_cline(int flag) {
	NCLINES *z;
	int col1=XNullColor,col2=YNullColor;
	if(n_nstore==0) {
		return;
	}
	z=ncperm;
	while(1) {
		if(z==NULL || (z->nmx==0 && z->nmy==0)) {
			return;
		}
		if(MyGraph->xv[0]==z->n_ix &&
		   MyGraph->yv[0]==z->n_iy &&
		   MyGraph->ThreeDFlag==0) {
			if(flag>0) {
				waitasec(flag);
				clr_scrn();
			}
			set_linestyle(col1);
			restor_null(z->xn,z->nmx,1);
			set_linestyle(col2);
			restor_null(z->yn,z->nmy,2);
			if(flag>0) {
				FlushDisplay();
			}
		}
		z=z->n;
		if(z==NULL) {
			break;
		}
	}
}


/* d=1 for x and 2 for y  */
static void restor_null(float *v, int n, int d) {
	int i,i4;
	float xm,ym;
	int x1,y1;
	if (PltFmtFlag==SVGFMT) {
		fprintf(svgfile,"<g>\n");
	}
	for(i=0;i<n;i++) {
		i4=4*i;
		line_abs(v[i4],v[i4+1],v[i4+2],v[i4+3]);
		if(NullStyle==1) {
			xm=.5*(v[i4]+v[i4+2]);
			ym=.5*(v[i4+1]+v[i4+3]);
			scale_to_screen(xm,ym,&x1,&y1);
			switch(d) {
			case 1:
				line(x1,y1-4,x1,y1+4);
				break;
			case 2:
				line(x1-4,y1,x1+4,y1);
				break;
			}
		}
	}
	if (PltFmtFlag==SVGFMT) {
		fprintf(svgfile,"</g>\n");
	}
}


static void save_frozen_clines(char *fn) {
	NCLINES *z;
	FILE *fp;
	char fnx[DEFAULT_STRING_LENGTH];
	char ch;
	int i=1;
	if(n_nstore==0) {
		return;
	}
	ch=(char)TwoChoice("YES","NO","Save Frozen Clines?","yn");
	if(ch=='n') {
		return;
	}
	z=ncperm;
	while(1) {
		if(z==NULL || (z->nmx==0 && z->nmy==0)) {
			return;
		}
		sprintf(fnx,"%s.%d",fn,i);
		fp=fopen(fnx,"w");
		if(fp==NULL) {
			err_msg("Cant open file!");
			return;
		}
		dump_clines(fp,z->xn,z->nmx,z->yn,z->nmy);
		fclose(fp);
		i++;
		z=z->n;
		if(z==NULL) {
			break;
		}
	}
}


/* animated nullclines stuff added Aug 31 97 just redraws them.
 * It will allow you to either freeze a range of them or just
 * one at a time freeze
 * - store the current set range_freeze
 * - compute over some range of parameters clear
 * - delete all but the current set animate
 * - replay all frozen ones (not current set )
 */
static void save_the_nullclines(void) {
	FILE *fp;
	char filename[DEFAULT_STRING_LENGTH];
	if(NULL_HERE==0) {
		return;
	}
	sprintf(filename,"nc.dat");
	ping();
	if(!file_selector("Save nullclines",filename,"*.dat")) {
		return;
	}
	fp=fopen(filename,"w");
	if(fp==NULL) {
		err_msg("Cant open file!");
		return;
	}
	dump_clines(fp,X_n,num_x_n,Y_n,num_y_n);
	fclose(fp);
	save_frozen_clines(filename);
}


static void start_ncline(void) {
	n_nstore=1;
	ncperm=(NCLINES *)malloc(sizeof(NCLINES));
	ncperm->p=NULL;
	ncperm->n=NULL;
	ncperm->nmx=0;
	ncperm->nmy=0;
	ncperm->n_ix=-5;
	ncperm->n_iy=-5;
	ncrange.xlo=0;
	ncrange.xhi=1;
	ncrange.nstep=10;
	sprintf(ncrange.rv," ");
}


static void stor_null(double x1, double y1, double x2, double y2) {
	int i;
	if(num_index>=MAX_NULL) {
		return;
	}
	i=4*num_index;
	saver[i]=x1;
	saver[i+1]=y1;
	saver[i+2]=x2;
	saver[i+3]=y2;
	num_index++;
}
