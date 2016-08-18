#include "rosenbrock.h"

#include <math.h>

#include "gear.h"
#include "../delay_handle.h"
#include "../flags.h"
#include "../load_eqn.h"
#include "../markov.h"
#include "../my_rhs.h"
#include "../util/matrixalgebra.h"
#include "../numerics.h"


/* --- Forward declarations --- */
static int one_flag_step_rosen(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr);

/* --- Functions --- */
/* This is rosen  - rosenbock step
 * This uses banded routines as well
 */
int rb23(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr) {
	int out =-1;
	if(NFlags==0) {
		out = rosen(y,tstart,tfinal,istart,n,work,ierr);
	} else {
		out = one_flag_step_rosen(y,tstart,tfinal,istart,n,work,ierr);
	}
	return(out);
}


int rosen(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr) {
	static double htry;
	double epsjac=NEWT_ERR;
	double eps=1e-15,hmin,hmax;
	double tdir=1,t0=*tstart,t=t0;
	double atol=ATOLER,rtol=TOLER;
	double sqrteps=sqrt(eps);
	double thresh=atol/rtol,absh,h;
	double d=1/(2.+sqrt(2.)),e32=6.+sqrt(2.),tnew;
	int i,n2=n*n,done=0,info,ml=cv_bandlower,mr=cv_bandupper,mt=ml+mr+1;
	int ipivot[MAXODE1],nofailed;
	double temp,err,tdel;
	double *ypnew,*k1,*k2,*k3,*f0,*f1,*f2,*dfdt,*ynew,*dfdy;

	*ierr=1;
	ypnew=work;
	k1=ypnew+n;
	k2=k1+n;
	k3=k2+n;
	f0=k3+n;
	f1=f0+n;
	f2=f1+n;
	dfdt=f2+n;
	ynew=dfdt+n;
	dfdy=ynew+n;

	if(t0>tfinal) {
		tdir=-1;
	}
	hmax=fabs(tfinal-t);
	if(*istart==1) {
		htry=hmax;
	}
	my_rhs(t0,y,f0,n);
	hmin=16*eps*fabs(t);
	absh = MIN(hmax, MAX(hmin, htry));
	while(!done) {
		nofailed=1;
		hmin = 16*eps*fabs(t);
		absh = MIN(hmax, MAX(hmin, absh));
		h = tdir * absh;
		if(1.1*absh >= fabs(tfinal - t)) {
			h = tfinal - t;
			absh = fabs(h);
			done = 1;
		}
		get_the_jac(t,y,f0,ypnew,dfdy,n,epsjac,1.0);
		tdel = (t + tdir*MIN(sqrteps*MAX(fabs(t),fabs(t+h)),absh)) - t;
		my_rhs(t+tdel,y,f1,n);
		for(i=0;i<n;i++) {
			dfdt[i]=(f1[i]-f0[i])/tdel;
		}
		while(1) { /* advance a step  */
			for(i=0;i<n2;i++) {
				dfdy[i]=-h*d*dfdy[i];
			}
			for(i=0;i<n;i++) {
				k1[i]=f0[i]+(h*d)*dfdt[i];
			}
			if(cv_bandflag) {
				for(i=0;i<n;i++) {
					dfdy[i*mt+ml]+=1;
				}
				bandfac(dfdy,ml,mr,n);
				bandsol(dfdy,k1,ml,mr,n);
			} else {
				for(i=0;i<n;i++) {
					dfdy[i*n+i]+=1;
				}
				sgefa(dfdy,n,n,ipivot,&info);
				sgesl(dfdy,n,n,ipivot,k1,0);
			}
			for(i=0;i<n;i++) {
				ynew[i]=y[i]+.5*h*k1[i];
			}
			my_rhs(t+.5*h,ynew,f1,n);
			for(i=0;i<n;i++) {
				k2[i]=f1[i]-k1[i];
			}
			if(cv_bandflag) {
				bandsol(dfdy,k2,ml,mr,n);
			} else {
				sgesl(dfdy,n,n,ipivot,k2,0);
			}
			for(i=0;i<n;i++) {
				k2[i]=k2[i]+k1[i];
				ynew[i]=y[i]+h*k2[i];
			}
			tnew=t+h;
			my_rhs(tnew,ynew,f2,n);
			for(i=0;i<n;i++) {
				k3[i]=f2[i] - e32*(k2[i] - f1[i]) - 2*(k1[i] - f0[i]) + (h*d)*dfdt[i];
			}
			if(cv_bandflag) {
				bandsol(dfdy,k3,ml,mr,n);
			} else {
				sgesl(dfdy,n,n,ipivot,k3,0);
			}
			err=0.0;
			for(i=0;i<n;i++) {
				temp=MAX(MAX(fabs(y[i]),fabs(ynew[i])),thresh);
				temp=fabs(k1[i]-2*k2[i]+k3[i])/temp;
				if(err<temp) {
					err=temp;
				}
			}
			err=err*(absh/6);
			if(err>rtol) {
				if(absh<hmin) {
					*ierr=-1;
					return(-1);
				}
				absh = MAX(hmin, absh * MAX(0.1, pow(0.8*(rtol/err),1./3.)));
				h = tdir * absh;
				nofailed=0;
				done=0;
			} else {
				break;
			}
		}
		if(nofailed==1) {
			temp=1.25*pow(err/rtol,1./3.);
			if(temp>0.2) {
				absh=absh/temp;
			} else {
				absh=5*absh;
			}
		}
		t=tnew;
		for(i=0;i<n;i++) {
			y[i]=ynew[i];
			f0[i]=f2[i];
		}
	}
	*tstart=t;
	htry=h;
	*istart=0;
	return(0);
}


/* --- Static functions --- */
static int one_flag_step_rosen(double *y, double *tstart, double tfinal, int *istart, int n, double *work, int *ierr) {
	double yold[MAXODE],told;
	int i,ok,hit;
	double s;
	int nstep=0;
	while(1) {
		for(i=0;i<n;i++) {
			yold[i]=y[i];
		}
		told=*tstart;
		ok=rosen(y,tstart,tfinal,istart,n,work,ierr);
		if(ok==-1) {
			break;
		}
		if((hit=one_flag_step(yold,y,istart,told,tstart,n,&s ))==0) {
			break;
		}
		/* Its a hit !! */
		nstep++;

		if(*tstart==tfinal) {
			break;
		}
		if(nstep>(NFlags+2)) {
			printf(" Working too hard? ");
			printf("smin=%g\n",s);
			*ierr=-2;
			return 1;
			break;
		}
	}
	return 0;
}

