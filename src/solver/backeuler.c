#include "backeuler.h"

#include <math.h>

#include "gear.h"
#include "../delay_handle.h"
#include "../flags.h"
#include "../load_eqn.h"
#include "../markov.h"
#include "../my_rhs.h"
#include "../numerics.h"
#include "../util/matrixalgebra.h"

/* --- Forward declarations --- */
static int one_flag_step_backeuler(double *y, double *t, double dt, int neq, double *yg, double *yp, double *yp2, double *ytemp, double *errvec, double *jac, int *istart);
static int one_step_backeuler(double *y, double *t, double dt, int neq, double *yg, double *yp, double *yp2, double *ytemp, double *errvec, double *jac, int *istart);

/* --- Functions --- */
int backwards_euler(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work) {
	int i,j;
	double *jac,*yg,*yp,*yp2,*ytemp,*errvec;
	yp=work;
	yg=yp+neq;
	ytemp=yg+neq;
	errvec=ytemp+neq;
	yp2=errvec+neq;
	jac=yp2+neq;
	if(NFlags==0) {
		for(i=0;i<nt;i++) {
			if((j=one_step_backeuler(y,tim,dt,neq,yg,yp,yp2,ytemp,errvec,jac,istart))!=0) {
				return(j);
			}
			stor_delay(y);
		}
		return(0);
	}
	for(i=0;i<nt;i++) {
		if((j=one_flag_step_backeuler(y,tim,dt,neq,yg,yp,yp2, ytemp,errvec,jac,istart))!=0) {
			return(j);
		}
		stor_delay(y);
	}
	return(0);
}


/* --- Static functions --- */
static int one_step_backeuler(double *y, double *t, double dt, int neq, double *yg, double *yp, double *yp2,
				 double *ytemp, double *errvec, double *jac, int *istart) {
	int i;
	double err=0.0,err1=0.0;

	int iter=0,info,ipivot[MAXODE1];
	int ml=cv_bandlower,mr=cv_bandupper,mt=ml+mr+1;
	set_wieners(dt,y,*t);
	*t=*t+dt;
	my_rhs(*t,y,yp2,neq);
	for(i=0;i<neq;i++) {
		yg[i]=y[i];
	}
	while(1) {
		err1=0.0;
		err=0.0;
		my_rhs(*t,yg,yp,neq);
		for(i=0;i<neq;i++) {
			errvec[i]=yg[i]-.5*dt*(yp[i]+yp2[i])-y[i];
			err1+=fabs(errvec[i]);
			ytemp[i]=yg[i];
		}
		get_the_jac(*t,yg,yp,ytemp,jac,neq,NEWT_ERR,-.5*dt);
		if(cv_bandflag) {
			for(i=0;i<neq;i++) {
				jac[i*mt+ml]+=1;
			}
			bandfac(jac,ml,mr,neq);
			bandsol(jac,errvec,ml,mr,neq);
		} else {
			for(i=0;i<neq;i++) {
				jac[i*neq+i]+=1.0;
			}
			sgefa(jac,neq,neq,ipivot,&info);
			if(info!=-1) {
				return(-1);
			}
			sgesl(jac,neq,neq,ipivot,errvec,0);
		}
		for(i=0;i<neq;i++) {
			err+=fabs(errvec[i]);
			yg[i]-=errvec[i];
		}
		if(err<EulTol || err1<EulTol) {
			for(i=0;i<neq;i++) {
				y[i]=yg[i];
			}
			return(0);
		}
		iter++;
		if(iter>MaxEulIter) {
			return(-2);
		}
	}
}


static int one_flag_step_backeuler(double *y, double *t, double dt, int neq, double *yg, double *yp,
						  double *yp2, double *ytemp, double *errvec, double *jac, int *istart) {
	double yold[MAXODE],told;
	int i,hit,j;
	double s;
	double dtt=dt;
	int nstep=0;
	while(1) {
		for(i=0;i<neq;i++) {
			yold[i]=y[i];
		}
		told=*t;
		if((j=one_step_backeuler(y,t,dtt,neq,yg,yp,yp2,ytemp,errvec,jac,istart))!=0) {
			return(j);
		}
		if((hit=one_flag_step(yold,y,istart,told,t,neq,&s ))==0) {
			break;
		}
		/* Its a hit !! */
		nstep++;
		dtt=(1-s)*dt;
		if(nstep>(NFlags+2)) {
			printf(" Working too hard?");
			printf("smin=%g\n",s);
			break;
		}
	}
	return 0;
}


