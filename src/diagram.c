#include "diagram.h"

#include <stdio.h>
#include <stdlib.h>

#include "autevd.h"
#include "auto_nox.h"
#include "browse.h"
#include "form_ode.h"
#include "load_eqn.h"
#include "markov.h"
#include "my_ps.h"

#define DALLOC(a) (double *)malloc((a)*sizeof(double))

int NBifs=0;
DIAGRAM *bifd;

void start_diagram(int n) {
	NBifs=1;
	bifd=(DIAGRAM *)malloc(sizeof(DIAGRAM));
	bifd->prev=NULL;
	bifd->next=NULL;
	bifd->index=0;
	bifd->uhi=DALLOC(n);
	bifd->ulo=DALLOC(n);
	bifd->u0=DALLOC(n);
	bifd->ubar=DALLOC(n);
	bifd->evr=DALLOC(n);
	bifd->evi=DALLOC(n);
	bifd->norm=0;
	bifd->lab=0;

	DiagFlag=0;
}

int find_diagram(int irs, int n, int *index, int *ibr, int *ntot,
				 int *itp, int *nfpar, double *a, double *uhi,
				 double * ulo, double *u0, double *par, double *per,
				 int *icp1, int *icp2, int *icp3, int *icp4) {
	int i,found=0;
	DIAGRAM *d;
	d=bifd;

	while(d->next!=NULL) {
		if(d->lab==irs) {
			found=1;
			break;
		}
		d=d->next;
	}
	if(found) {
		*ibr=d->ibr;
		*ntot=d->ntot;
		*index=d->index;
		*itp=d->itp;
		*nfpar=d->nfpar;
		*a=d->norm;
		par=d->par;
		*icp1=d->icp1;
		*icp2=d->icp2;
		*icp3=d->icp3;
		*icp4=d->icp4;
		*per=d->per;
		for(i=0;i<n;i++) {
			u0[i]=d->u0[i];
			ulo[i]=d->ulo[i];
			uhi[i]=d->uhi[i];
		}
		return(1);
	}
	return(0);
}

void edit_start(int ibr, int ntot, int itp, int lab, int nfpar,
				double a, double *uhi, double *ulo, double *u0,
				double *ubar, double *par, double per, int n,
				int icp1, int icp2, int icp3, int icp4,
				double *evr, double *evi) {
	edit_diagram(bifd,ibr,ntot,itp,lab,nfpar,a,uhi,ulo,u0,ubar,
				 par,per,n,icp1,icp2,icp3,icp4,AutoTwoParam,evr,evi,blrtn.torper);
}

void edit_diagram(DIAGRAM *d, int ibr, int ntot, int itp, int lab, int nfpar,
				  double a, double *uhi, double *ulo, double *u0, double *ubar,
				  double *par, double per, int n, int icp1, int icp2, int icp3,
				  int icp4,int flag2, double *evr, double *evi, double tp) {
	int i;
	d->calc=TypeOfCalc;
	d->ibr=ibr;
	d->ntot=ntot;
	d->itp=itp;
	d->lab=lab;
	d->nfpar=nfpar;
	d->norm=a;
	for(i=0;i<8;i++) {
		d->par[i]=par[i];
		/*  printf("%d %g\n",i,par[i]); */
	}
	d->per=per;
	d->icp1=icp1;
	d->icp2=icp2;
	d->icp3=icp3;
	d->icp4=icp4;
	d->flag2=flag2;
	for(i=0;i<n;i++) {
		d->ulo[i]=ulo[i];
		d->uhi[i]=uhi[i];
		d->ubar[i]=ubar[i];
		d->u0[i]=u0[i];
		d->evr[i]=evr[i];
		d->evi[i]=evi[i];
	}
	d->torper=tp;
}

void add_diagram(int ibr, int ntot, int itp, int lab, int nfpar, double a,
				 double *uhi, double *ulo, double *u0, double *ubar, double *par,
				 double per, int n, int icp1, int icp2, int icp3,int icp4,int flag2,
				 double *evr, double *evi) {
	DIAGRAM *d,*dnew;

	d=bifd;
	while(d->next != NULL) {
		d=(d->next);
	}
	d->next=(DIAGRAM *)malloc(sizeof(DIAGRAM));
	dnew=d->next;
	dnew->next=NULL;
	dnew->prev=d;
	dnew->uhi=DALLOC(n);
	dnew->ulo=DALLOC(n);
	dnew->u0=DALLOC(n);
	dnew->ubar=DALLOC(n);
	dnew->evr=DALLOC(n);
	dnew->evi=DALLOC(n);
	dnew->index=NBifs;
	NBifs++;
	edit_diagram(dnew,ibr,ntot,itp,lab,nfpar,a,uhi,ulo,u0,ubar,par,per,n,
				 icp1,icp2,icp3,icp4,flag2,evr,evi,blrtn.torper);
}


void kill_diagrams(void) {
	DIAGRAM *d,*dnew;
	d=bifd;
	while(d->next != NULL) {  /*  Move to the end of the tree  */
		d=d->next;
	}
	while(d->prev != NULL ) {
		dnew=d->prev;
		d->next=NULL;
		d->prev=NULL;
		free(d->uhi);
		free(d->ulo);
		free(d->u0);
		free(d->ubar);
		free(d->evr);
		free(d->evi);
		free(d);
		d=dnew;
	}

	free(bifd->uhi);
	free(bifd->ulo);
	free(bifd->u0);
	free(bifd->ubar);
	free(bifd->evr);
	free(bifd->evi);
	free(bifd);
	start_diagram(NODE);
}

void redraw_diagram(void) {
	DIAGRAM *d;
	int type,flag=0;
	draw_bif_axes();
	d=bifd;
	if(d->next==NULL) {
		return;
	}
	while(1) {
		type=get_bif_type(d->ibr,d->ntot,d->lab);

		if(d->ntot==1) {
			flag=0;
		} else {
			flag=1;
		}
		add_point(d->par,d->per,d->uhi,d->ulo,d->ubar,d->norm,type,flag,
				  d->lab,d->nfpar,d->icp1,d->icp2,d->icp3,d->icp4,d->flag2,d->evr,d->evi);
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
}

void write_info_out(void) {
	char filename[XPP_MAX_NAME];
	DIAGRAM *d;
	int type,i;
	int status;
	int icp1,icp2,icp3,icp4;
	double *par;
	double par1,par2=0,*uhigh,*ulow,per;
	FILE *fp;
	sprintf(filename,"allinfo.dat");
	status=file_selector("Write all info",filename,"*.dat");

	if(status==0) {
		return;
	}
	fp=fopen(filename,"w");
	if(fp==NULL) {
		err_msg("Can't open file");
		return;
	}

	d=bifd;
	if(d->next==NULL) {
		return;
	}
	while(1) {
		type=get_bif_type(d->ibr,d->ntot,d->lab);
		icp1=d->icp1;
		icp2=d->icp2;
		par=d->par;
		per=d->per;
		uhigh=d->uhi;
		ulow=d->ulo;
		par1=par[icp1];
		if(icp2<NAutoPar) {
			par2=par[icp2];
		} else {
			par2=par1;
		}
		fprintf(fp,"%d %d %d %g %g %g ",
				type,d->ibr,d->flag2,par1,par2,per);
		for(i=0;i<NODE;i++) {
			fprintf(fp,"%g ",uhigh[i]);
		}
		for(i=0;i<NODE;i++) {
			fprintf(fp,"%g ",ulow[i]);
		}
		for(i=0;i<NODE;i++) {
			fprintf(fp,"%g %g ",d->evr[i],d->evi[i]);
		}
		fprintf(fp,"\n");
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
	fclose(fp);
}

void load_browser_with_branch(int ibr,int pts,int pte) {
	DIAGRAM *d;
	int type,i,j,pt;
	int icp1,icp2;
	double *par;
	double par1,par2=0,*u0,per;
	int first,last,nrows;
	first=abs(pts);
	last=abs(pte);
	if(last<first) { /* reorder the points so that we will store in right range*/
		i=first;
		first=last;
		last=i;
	}
	nrows=last-first+1;
	d=bifd;
	if(d->next==NULL) {
		return;
	}
	j=0;
	while(1) {
		type=get_bif_type(d->ibr,d->ntot,d->lab);
		pt=abs(d->ntot);
		if((d->ibr==ibr) && (pt>=first) && (pt<=last)) {
			icp1=d->icp1;
			icp2=d->icp2;
			par=d->par;
			per=d->per;
			u0=d->u0;

			par1=par[icp1];
			if(icp2<NAutoPar) {
				par2=par[icp2];
			} else {
				par2=par1;
			}
			storage[0][j]=par1;
			for(i=0;i<NODE;i++) {
				storage[i+1][j]=u0[i];
			}
			j++;
		}
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
	storind=nrows;
	refresh_browser(nrows);
}


void write_init_data_file(void) {
	char filename[XPP_MAX_NAME];
	DIAGRAM *d;
	int type,i;
	int status;
	int icp1,icp2;
	double *par;
	double par1,par2=0,*u0,per;
	FILE *fp;
	sprintf(filename,"initdata.dat");
	status=file_selector("Write init data file",filename,"*.dat");

	if(status==0) {
		return;
	}

	fp=fopen(filename,"w");
	if(fp==NULL) {
		err_msg("Can't open file");
		return;
	}

	d=bifd;
	if(d->next==NULL) {
		return;
	}

	while(1) {
		type=get_bif_type(d->ibr,d->ntot,d->lab);
		icp1=d->icp1;
		icp2=d->icp2;
		par=d->par;
		per=d->per;
		u0=d->u0;
		par1=par[icp1];
		if(icp2<NAutoPar) {
			par2=par[icp2];
		} else {
			par2=par1;
		}
		fprintf(fp,"%d %d %g %g %g ",
				type,d->ibr,par1,par2,per);
		for(i=0;i<NODE;i++) {
			fprintf(fp,"%g ",u0[i]);
		}
		fprintf(fp,"\n");
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
	fclose(fp);
}


void write_pts(void) {
	char filename[XPP_MAX_NAME];
	DIAGRAM *d;
	int type;
	int status;
	int icp1,icp2;
	double *par;
	double x,y1,y2,par1,par2=0,a,*uhigh,*ulow,*ubar,per;
	FILE *fp;
	sprintf(filename,"diagram.dat");
	status=file_selector("Write points",filename,"*.dat");
	if(status==0) {
		return;
	}
	fp=fopen(filename,"w");
	if(fp==NULL) {
		err_msg("Can't open file");
		return;
	}

	d=bifd;
	if(d->next==NULL) {
		return;
	}
	while(1) {
		type=get_bif_type(d->ibr,d->ntot,d->lab);
		icp1=d->icp1;
		icp2=d->icp2;
		par=d->par;
		per=d->per;
		uhigh=d->uhi;
		ulow=d->ulo;
		ubar=d->ubar;
		a=d->norm;
		par1=par[icp1];
		if(icp2<NAutoPar) {
			par2=par[icp2];
		}
		if(check_plot_type(d->flag2,icp1,icp2)==1) {
			auto_xy_plot(&x,&y1,&y2,par1,par2,per,uhigh,ulow,ubar,a);
			fprintf(fp,"%g %g %g %d %d %d\n",
					x,y1,y2,type,abs(d->ibr),d->flag2);
		}
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
	fclose(fp);
}

void post_auto(void) {
	char filename[XPP_MAX_NAME];
	DIAGRAM *d;
	int type,flag=0;
	int status;
	sprintf(filename,"auto.ps");
	status=file_selector("Postscript",filename,"*.ps");

	if(status==0) {
		return;
	}
	if(!ps_init(filename,PS_Color)) {
		return;
	}
	draw_ps_axes();
	d=bifd;
	if(d->next==NULL) {
		return;
	}

	while(1) {
		type=get_bif_type(d->ibr,d->ntot,d->lab);
		if(type < 0) {
			plintf("Unable to get bifurcation type.\n");
		}
		if(d->ntot==1) {
			flag=0;
		} else {
			flag=1;
		}
		add_ps_point(d->par,d->per,d->uhi,d->ulo,d->ubar,d->norm,type,flag,
					 d->lab,d->nfpar,d->icp1,d->icp2,d->flag2,d->evr,d->evi);
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
	ps_end();
	set_normal_scale();
}


void svg_auto(void) {
	char filename[XPP_MAX_NAME];
	DIAGRAM *d;
	int type,flag=0;
	int status;
	sprintf(filename,"auto.svg");
	status=file_selector("SVG",filename,"*.svg");
	if(status==0) {
		return;
	}
	if(!svg_init(filename,PS_Color)) {
		return;
	}
	draw_svg_axes();
	d=bifd;
	if(d->next==NULL) {
		return;
	}
	while(1) {
		type=get_bif_type(d->ibr,d->ntot,d->lab);
		if(type < 0) {
			plintf("Unable to get bifurcation type.\n");
		}
		if(d->ntot==1) {
			flag=0;
		} else {
			flag=1;
		}
		add_ps_point(d->par,d->per,d->uhi,d->ulo,d->ubar,d->norm,type,flag,
					 d->lab,d->nfpar,d->icp1,d->icp2,d->flag2,d->evr,d->evi);
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
	svg_end();

	set_normal_scale();
}


void bound_diagram(double *xlo, double *xhi, double *ylo, double *yhi) {
	DIAGRAM *d;
	int type;

	/*int flag=0;
	  */
	double x,y1,y2,par1,par2=0.0;
	d=bifd;
	if(d->next==NULL) {
		return;
	}
	*xlo=1.e16;
	*ylo=*xlo;
	*xhi=-*xlo;
	*yhi=-*ylo;
	while(1) {
		type=get_bif_type(d->ibr,d->ntot,d->lab);
		if(type <1) {
			plintf("Unable to get bifurcation type.\n");
		}
		par1=d->par[d->icp1];
		if(d->icp2<NAutoPar) {
			par2=d->par[d->icp2];
		}
		auto_xy_plot(&x,&y1,&y2,par1,par2,d->per,d->uhi,d->ulo,d->ubar,d->norm);
		if(x<*xlo) {
			*xlo=x;
		}
		if(x>*xhi) {
			*xhi=x;
		}
		if(y2<*ylo) {
			*ylo=y2;
		}
		if(y1>*yhi) {
			*yhi=y1;
		}
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
}


int save_diagram(FILE *fp, int n) {
	int i;
	DIAGRAM *d;
	fprintf(fp,"%d\n",NBifs-1);
	if(NBifs==1) {
		return(-1);
	}
	d=bifd;
	while(1) {
		fprintf(fp,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
				d->calc,d->ibr,d->ntot,d->itp,d->lab,d->index,d->nfpar,
				d->icp1,d->icp2,d->icp3,d->icp4,d->flag2);
		for(i=0;i<8;i++) {
			fprintf(fp,"%g ",d->par[i]);
		}
		fprintf(fp,"%g %g \n",d->norm,d->per);

		for(i=0;i<n;i++) {
			fprintf(fp,"%f %f %f %f %f %f\n",d->u0[i],d->uhi[i],d->ulo[i],
					d->ubar[i],d->evr[i],d->evi[i]);
		}
		d=d->next;
		if(d==NULL) {
			break;
		}
	}
	return(1);
}


int load_diagram(FILE *fp, int node) {
	double u0[NAUTO],uhi[NAUTO],ulo[NAUTO],ubar[NAUTO],evr[NAUTO],evi[NAUTO],norm,par[8],per;
	int i,flag=0;
	int n;
	int calc,ibr,ntot,itp,lab,index,nfpar,icp1,icp2,icp3,icp4,flag2;
	fscanf(fp,"%d",&n);
	if(n==0) {
		/*    start_diagram(NODE); */
		return(-1);
	}

	while(1) {
		fscanf(fp,"%d %d %d %d %d %d %d %d %d %d %d %d",
			   &calc,&ibr,&ntot,&itp,&lab,&index,&nfpar,
			   &icp1,&icp2,&icp3,&icp4,&flag2);
		for(i=0;i<8;i++) {
			fscanf(fp,"%lg ",&par[i]);
		}
		fscanf(fp,"%lg %lg ",&norm,&per);
		for(i=0;i<node;i++) {
			fscanf(fp,"%lg %lg %lg %lg %lg %lg",&u0[i],&uhi[i],&ulo[i],
				   &ubar[i],&evr[i],&evi[i]);
		}
		if(flag==0) {
			edit_start(ibr,ntot,itp,lab,nfpar,norm,uhi,ulo,u0,ubar,par,per,node,
					   icp1,icp2,icp3,icp4,evr,evi);
			flag=1;
			DiagFlag=1;
		} else {
			add_diagram(ibr,ntot,itp,lab,nfpar,norm,uhi,ulo,u0,ubar,par,per,node,
						icp1,icp2,icp3,icp4,flag2,evr,evi);
		}
		if(index>=n) {
			break;
		}
	}
	return(1);
}
