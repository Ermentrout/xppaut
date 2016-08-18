#include "symplect.h"

#include "../delay_handle.h"
#include "../flags.h"
#include "../my_rhs.h"

/* --- Forward declarations --- */
static int one_flag_step_symp(double *y, double dt, double *work, int neq, double *tim, int *istart);
static void one_step_symp(double *y, double h, double *f, int n, double *t);

/* --- Data --- */
static double symp_b[]={7/24.,.75,-1./24};
static double symp_B[]={2/3.,-2./3.,1.0};

/* my first symplectic integrator */
int symplect3(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work) {
	int i;
	if(NFlags==0) {
		for(i=0;i<nt;i++) {
			one_step_symp(y,dt,work,neq,tim);
		}
		stor_delay(y);
		return(0);
	}
	for(i=0;i<nt;i++) {
		one_flag_step_symp(y,dt,work,neq,tim,istart);
		stor_delay(y);
	}
	return(0);
}


static int one_flag_step_symp(double *y, double dt, double *work, int neq, double *tim, int *istart) {
	double yold[MAXODE],told;
	int i,hit;
	double s,dtt=dt;
	int nstep=0;
	while(1) {
		for(i=0;i<neq;i++) {
			yold[i]=y[i];
		}
		told=*tim;
		one_step_symp(y,dtt,work,neq,tim);
		if((hit=one_flag_step(yold,y,istart,told,tim,neq,&s ))==0) {
			break;
		}
		/* Its a hit !! */
		nstep++;
		dtt=(1-s)*dt;
		if(nstep>(NFlags+2)) {
			printf(" Working too hard?? ");
			printf("smin=%g\n",s);
			break;
		}
	}
	return(1);
}


static void one_step_symp(double *y, double h, double *f, int n, double *t) {
	int s,j;
	for(s=0;s<3;s++) {
		for(j=0;j<n;j+=2) {
			y[j]+=(h*symp_b[s]*y[j+1]);
		}
		my_rhs(*t,y,f,n);
		for(j=0;j<n;j+=2) {
			y[j+1]+=(h*symp_B[s]*f[j+1]);
		}
	}
	*t+=h;
}



