#include "runge_kutta.h"

#include "../delay_handle.h"
#include "../flags.h"
#include "../markov.h"
#include "../my_rhs.h"


/* --- Forward declarations --- */
static int one_flag_step_rk4(double *y, double dt, double *yval[3], int neq, double *tim, int *istart);
static void one_step_rk4(double *y, double dt, double *yval[3], int neq, double *tim);


/* --- Functions --- */
int runge_kutta(double *y, double *tim, double dt, int nt, int neq, int *istart, double *work) {
	register int j;
	double *yval[3];

	yval[0]=work;
	yval[1]=work+neq;
	yval[2]=work+neq+neq;

	if(NFlags==0) {
		for(j=0;j<nt;j++) {
			one_step_rk4(y,dt,yval,neq,tim);
			stor_delay(y);
		}
		return(0);
	}
	for(j=0;j<nt;j++) {
		one_flag_step_rk4(y,dt,yval,neq,tim,istart);
		stor_delay(y);
	}
	return(0);
}


/* --- Static functions --- */
static int one_flag_step_rk4(double *y, double dt, double *yval[3], int neq, double *tim, int *istart) {
	double yold[MAXODE],told;
	int i,hit;
	double s,dtt=dt;
	int nstep=0;
	while(1) {
		for(i=0;i<neq;i++) {
			yold[i]=y[i];
		}
		told=*tim;
		one_step_rk4(y,dtt,yval,neq,tim);
		if((hit=one_flag_step(yold,y,istart,told,tim,neq,&s ))==0) {
			break;
		}
		/* Its a hit !! */
		nstep++;
		dtt=(1-s)*dt;
		if(nstep>(NFlags+2)) {
			printf(" Working too hard?");
			printf("smin=%g\n",s);
			/* plintflaginfo(); */
			break;
		}
	}
	return(1);
}


static void one_step_rk4(double *y, double dt, double *yval[3], int neq, double *tim) {
	int i;
	double t=*tim,t1,t2;
	set_wieners(dt,y,t);
	my_rhs(t,y,yval[1],neq);
	for(i=0;i<neq;i++) {
		yval[0][i]=y[i]+dt*yval[1][i]/6.00;
		yval[2][i]=y[i]+dt*yval[1][i]*0.5;
	}
	t1=t+.5*dt;
	my_rhs(t1,yval[2],yval[1],neq);
	for(i=0;i<neq;i++) {
		yval[0][i]=yval[0][i]+dt*yval[1][i]/3.00;
		yval[2][i]=y[i]+.5*dt*yval[1][i];
	}
	my_rhs(t1,yval[2],yval[1],neq);
	for(i=0;i<neq;i++) {
		yval[0][i]=yval[0][i]+dt*yval[1][i]/3.000;
		yval[2][i]=y[i]+dt*yval[1][i];
	}
	t2=t+dt;
	my_rhs(t2,yval[2],yval[1],neq);
	for(i=0;i<neq;i++) {
		y[i]=yval[0][i]+dt*yval[1][i]/6.00;
	}
	*tim=t2;
}


