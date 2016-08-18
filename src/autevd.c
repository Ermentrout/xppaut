#include "autevd.h"
#undef abs

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auto_nox.h"
#include "diagram.h"
#include "solver/gear.h"
#include "xAuto.h"

/* --- Macros --- */
#define MAXDIMHET 12
#define MAXDIMHOM 24


/* --- Types --- */
typedef struct {
	int pt,br;
	double evr[NAUTO],evi[NAUTO];
} EIGVAL;


/* --- Data --- */
int DiagFlag=0;
static EIGVAL my_ev;
XAUTO xAuto;


/* --- Functions --- */
void init_auto(int ndim, int nicp,int nbc, int ips, int irs, int ilp,
			   int ntst, int isp, int isw, int nmx, int npr,
			   double ds, double dsmin, double dsmax, double rl0, double rl1,
			   double a0, double a1, int ip1, int ip2, int ip3, int ip4, int ip5,
			   int nuzr, double epsl, double epsu, double epss, int ncol) {
	/* here are the constants that we do not allow the user to change */
	int nnbc;
	int i;
	xAuto.iad=aauto.iad;
	xAuto.iplt=0;
	xAuto.mxbf=aauto.mxbf;
	xAuto.iid=aauto.iid;
	xAuto.itmx=aauto.itmx;
	xAuto.itnw=aauto.itnw;
	xAuto.nwtn=aauto.nwtn;
	xAuto.jac=0;
	xAuto.iads=aauto.iads;
	xAuto.nthl=1;
	xAuto.ithl[0]=10;
	xAuto.thl[0]=0.0;
	xAuto.nint=0;

	if(ips==4) {
		nnbc=ndim;
	} else {
		nnbc=0;
	}
	xAuto.ndim=ndim;
	xAuto.nbc=nnbc;
	xAuto.ips=ips;
	xAuto.irs=irs;
	xAuto.ilp=ilp;
	xAuto.nicp=nicp;
	xAuto.icp[0]=ip1;
	xAuto.icp[1]=ip2;
	xAuto.icp[2]=ip3;
	xAuto.icp[3]=ip4;
	xAuto.icp[4]=ip5;
	xAuto.ntst=ntst;
	xAuto.ncol=ncol;
	xAuto.isp=isp;
	xAuto.isw=isw;
	xAuto.nmx=nmx;
	xAuto.rl0=rl0;
	xAuto.rl1=rl1;
	xAuto.a0=a0;
	xAuto.a1=a1;
	xAuto.npr=npr;
	xAuto.epsl=epsl;
	xAuto.epss=epss;
	xAuto.epsu=epsu;
	xAuto.ds=ds;
	xAuto.dsmax=dsmax;
	xAuto.dsmin=dsmin;
	xAuto.nuzr=NAutoUzr;
	for(i=0;i<NAutoUzr;i++) {
		xAuto.iuz[i]=UzrPar[i];
		xAuto.vuz[i]=outperiod[i];
	}
}


void send_eigen(int ibr, int ntot, int n, doublecomplex *ev) {
	int i;
	double er,cs,sn;

	my_ev.pt=abs(ntot);
	my_ev.br=abs(ibr);
	for(i=0;i<n;i++) {
		er=exp((ev+i)->r);
		cs=cos((ev+i)->i);
		sn=sin((ev+i)->i);
		my_ev.evr[i]=er*cs;
		my_ev.evi[i]=er*sn;
	}
}

void send_mult(int ibr, int ntot, int n, doublecomplex *ev) {
	int i;

	my_ev.pt=abs(ntot);
	my_ev.br=abs(ibr);
	for(i=0;i<n;i++) {
		my_ev.evr[i]=(ev+i)->r;
		my_ev.evi[i]=(ev+i)->i;
	}
}


/* Only unit 8,3 or q.prb is important; all others are unnecesary */
int get_bif_type(int ibr, int ntot, int lab) {

	int type=STABLE_EQUILIBRIUM;
	if(ibr<0 && ntot<0) {
		type=STABLE_PERIODIC;
	}
	if(ibr<0 && ntot>0) {
		type=UNSTABLE_PERIODIC;
	}
	if(ibr>0 && ntot>0) {
		type=UNSTABLE_EQUILIBRIUM;
	}
	if(ibr>0 && ntot<0) {
		type=STABLE_EQUILIBRIUM;
	}
	return type;
}


void addbif(iap_type *iap, rap_type *rap, integer ntots, integer ibrs,
			double *par,integer *icp,int lab, double *a,
			double *uhigh, double *ulow, double *u0, double *ubar) {
	int type;
	int icp1=icp[0],icp2=icp[1],icp3=icp[2],icp4=icp[3];
	double    per=par[10];

	type=get_bif_type(ibrs,ntots,lab);
	if(iap->ntot==1) {
		add_point(par,per,uhigh,ulow,ubar,*a,type,0,lab,
				  iap->nfpr,icp1,icp2,icp3,icp4,AutoTwoParam,my_ev.evr,my_ev.evi);
	} else {
		add_point(par,per,uhigh,ulow,ubar,*a,type,1,lab,
				  iap->nfpr,icp1,icp2,icp3,icp4,AutoTwoParam,my_ev.evr,my_ev.evi);
	}
	if(DiagFlag==0) {
		/* start_diagram(*ndim); */
		edit_start(ibrs,ntots,iap->itp,lab,iap->nfpr,*a,uhigh,ulow,u0,ubar,
				   par,per,iap->ndim,icp1,icp2,icp3,icp4,my_ev.evr,my_ev.evi);
		DiagFlag=1;
		return;
	}
	add_diagram(ibrs,ntots,iap->itp,lab,iap->nfpr,*a,uhigh,ulow,u0,ubar,
				par,per,iap->ndim,icp1,icp2,icp3,icp4,AutoTwoParam,my_ev.evr,
				my_ev.evi);
}
