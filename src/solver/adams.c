#include "adams.h"

#include "runge_kutta.h"
#include "../delay_handle.h"
#include "../flags.h"
#include "../markov.h"
#include "../my_rhs.h"

/* --- Data --- */
static const double coefp[] = {6.875 / 3.00, -7.375 / 3.00, 4.625 / 3.00,
							   -.375};
static const double coefc[] = {.375, 2.375 / 3.00, -.625 / 3.00, 0.125 / 3.00};
static double *y_s[4], *y_p[4], *ypred;


/* --- Forward declarations --- */
static int abmpc(double *y, double *t, double dt, int neq);


/* --- Functions --- */
int adams(double *y, double *tim, double dt, int nstep, int neq, int *ist, double *work) {
	int istart=*ist,i,istpst,k,ik,n;
	int irk;
	double *work1;
	double x0=*tim,xst=*tim;
	work1=work;
	if(istart==1) {
		for(i=0;i<4;i++) {
			y_p[i]=work+(4+i)*neq;
			y_s[i]=work+(8+i)*neq;
		}
		ypred=work+3*neq;
		goto n20;
	}
	if(istart>1) {
		goto n350;
	}
	istpst=0;
	goto n400;
n20:
	x0=xst;
	my_rhs(x0,y,y_p[3],neq);
	for(k=1;k<4;k++) {
		runge_kutta(y,&x0,dt,1,neq,&irk,work1);
		stor_delay(y);
		for(i=0;i<neq;i++) {
			y_s[3-k][i]=y[i];
		}
		my_rhs(x0,y,y_p[3-k],neq);
	}
	istpst=3;
	if(istpst<=nstep) {
		goto n400;
	}
	ik=4-nstep;
	for(i=0;i<neq;i++) {
		y[i]=y_s[ik-1][i];
	}
	xst=xst+nstep*dt;
	istart=ik;
	goto n1000;
n350:
	ik=istart-nstep;
	if(ik<=1) {
		goto n370;
	}
	for(i=0;i<neq;i++) {
		y[i]=y_s[ik-1][i];
	}
	xst=xst+nstep*dt;
	istart=ik;
	goto n1000;
n370:
	for(i=0;i<neq;i++) {
		y[i]=y_s[0][i];
	}
	if(ik==1) {
		x0=xst+dt*nstep;
		goto n450;
	}
	istpst=istart-1;
n400:
	if(istpst==nstep) {
		goto n450;
	}
	for(n=istpst+1;n<nstep+1;n++) {
		set_wieners(dt,y,x0);
		abmpc(y,&x0,dt,neq);
		stor_delay(y);
	}
n450:
	istart=0;
	xst=x0;
n1000:
	*tim=*tim + nstep*dt;
	*ist=istart;
	return(0);
}


/* --- Static functions --- */
static int abmpc(double *y, double *t, double dt, int neq) {
	double x1,x0=*t;
	int i,k;
	for(i=0;i<neq;i++) {
		ypred[i]=0;
		for(k=0;k<4;k++) {
			ypred[i]=ypred[i]+coefp[k]*y_p[k][i];
		}
		ypred[i]=y[i]+dt*ypred[i];
	}
	for(i=0;i<neq;i++) {
		for(k=3;k>0;k--) {
			y_p[k][i]=y_p[k-1][i];
		}
	}
	x1=x0+dt;
	my_rhs(x1,ypred,y_p[0],neq);
	for(i=0;i<neq;i++) {
		ypred[i]=0;
		for(k=0;k<4;k++) {
			ypred[i]=ypred[i]+coefc[k]*y_p[k][i];
		}
		y[i]=y[i]+dt*ypred[i];
	}
	*t=x1;
	my_rhs(x1,y,y_p[0],neq);
	return(1);
}
