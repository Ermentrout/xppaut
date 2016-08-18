#include "euler.h"

#include "../delay_handle.h"
#include "../flags.h"
#include "../markov.h"
#include "../my_rhs.h"
#include "../util/matrixalgebra.h"

/* --- Forward declarations --- */
static int one_flag_step_euler(double *y, double dt, double *work, int neq, double *tim, int *istart);
static void one_step_euler(double *y, double dt, double *yp, int neq, double *t);


/* --- Functions --- */
int euler(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work) {
	int i;
	if(NFlags==0) {
		for(i=0;i<nt;i++) {
			one_step_euler(y,dt,work,neq,tim);
			stor_delay(y);
		}
		return(0);
	}
	for(i=0;i<nt;i++) {
		one_flag_step_euler(y,dt,work,neq,tim,istart);
		stor_delay(y);
	}
	return(0);
}


/* --- Static functions --- */
static int one_flag_step_euler(double *y, double dt, double *work, int neq, double *tim, int *istart) {
	double yold[MAXODE],told;
	int i,hit;
	double s,dtt=dt;
	int nstep=0;
	while(1) {
		for(i=0;i<neq;i++) {
			yold[i]=y[i];
		}
		told=*tim;
		one_step_euler(y,dtt,work,neq,tim);
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


static void one_step_euler(double *y, double dt, double *yp, int neq, double *t) {
	int j;
	set_wieners(dt,y,*t);
	my_rhs(*t,y,yp,neq);
	*t+=dt;
	for(j=0;j<neq;j++) {
		y[j]=y[j]+dt*yp[j];
	}
}
