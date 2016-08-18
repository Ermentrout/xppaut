#include "find_fixedPoint.h"

#include <math.h>
#include <stdlib.h>

#include "abort.h"
#include "eig_list.h"
#include "form_ode.h"
#include "ggets.h"
#include "graphics.h"
#include "integrate.h"
#include "load_eqn.h"
#include "markov.h"
#include "menudrive.h"
#include "numerics.h"
#include "util/matrixalgebra.h"

/* --- Forward declarations --- */
static void pr_evec(double *x, double *ev, int n, int pr, double eval,int type);


/* --- Data --- */
double ShootIC[8][MAXODE];
int ShootICFlag;
int ShootIndex;
int ShootType[8];

/* --- Functions --- */
/* main fixed point finder */
void do_sing(double *x, double eps, double err, double big, int maxit, int n, int *ierr, float *stabinfo) {
	int kmem,i,j,ipivot[MAXODE];
	int oldcol,dummy;
	int rp=0,rn=0,cp=0,cn=0,im=0;
	int pose=0,nege=0,pr;
	double *work,*eval,*b,*bp,*oldwork,*ework;
	double temp,oldt=DELTA_T,old_x[MAXODE];

	char ch;
	double real,imag;
	double bigpos=-1e10,bigneg=1e10;
	int bpos=0,bneg=0;
	/* float xl[MAXODE]; */
	kmem=n*(2*n+5)+50;
	work=(double *)malloc(sizeof(double)*kmem);
	if(work==NULL) {
		err_msg("Insufficient core ");
		return;
	}
	ShootICFlag=0;
	ShootIndex=0;
	for(i=0;i<n;i++) {
		old_x[i]=x[i];
	}
	oldwork=work+n*n;
	eval=oldwork+n*n;
	b=eval+2*n;
	bp=b+n;
	ework=bp+n;
	rooter(x,err,eps,big,work,ierr,maxit,n);
	if(*ierr!=0) {
		free(work);
		err_msg("Could not converge to root");
		for(i=0;i<n;i++) {
			x[i]=old_x[i];
		}
		return;
	}
	ping();
	for(i=0;i<n*n;i++) {
		oldwork[i]=work[i];
	}
	/* Transpose for Eigen        */
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			temp=work[i+j*n];
			work[i+j*n]=work[i*n+j];
			work[i*n+j]=temp;
		}
	}
	eigen(n,work,eval,ework,ierr);
	if(*ierr!=0) {
		err_msg("Could not compute eigenvalues");
		free(work);
		return;
	}
	/* succesfully computed evals now lets work with them */
	ch='n';
	if(!PAR_FOL) {
		ch=(char)TwoChoice("YES","NO","Print eigenvalues?","yn");
	}
	pr=0;

	if(ch=='y') {
		printf("\n Eigenvalues:\n");
		pr=1;
	}
	for(i=0;i<n;i++) {
		real=eval[2*i];
		imag=eval[2*i+1];
		if(pr==1) {
			printf(" %f  +  i  %f \n",real,imag);
		}
		if(METHOD == METHOD_DISCRETE) {
			real=real*real+imag*imag-1.00;
		}
		if(fabs(imag)<.00000001) {
			imag=0.0;
		}
		if(real<0.0) {
			if(imag!=0.0) {
				cn++;
				if(real<bigneg) {
					bigneg=real;
					bneg=-1;
				}
			} else {
				rn++;
				nege=i;
				if(real<bigneg) {
					bigneg=real;
					bneg=i;
				}
			}
		}
		if(real>0.0) {
			if(imag!=0.0) {
				cp++;
				if(real>bigpos) {
					bigpos=real;
					bpos=-1;
				}
			} else {
				rp++;
				pose=i;
				if(real>bigpos) {
					bigpos=real;
					bpos=i;
				}
			}
		}
		if((real==0.0) && (imag!=0.0)) {
			im++;
		}
	}     /* eigenvalue count */
	if(((rp+cp)!=0) && ((rn+cn)!=0)) {
		eq_symb(x,1);
	} else {
		if((rp+cp)!=0) {
			eq_symb(x,0);
		} else {
			eq_symb(x,3);
		}
	}
	*stabinfo=(float)(cp+rp)+(float)(cn+rn)/1000.0;

	/* Lets change Work back to transposed oldwork */
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			temp=oldwork[i+j*n];
			work[i+j*n]=oldwork[i*n+j];
			work[i*n+j]=temp;
		}
	}
	create_eq_box(cp,cn,rp,rn,im,x,eval,n);
	if(((rp==1) || (rn==1)) && (n>1)) {
		ch='n';
		if(!PAR_FOL) {
			ch=(char)TwoChoice("YES","NO","Draw Invariant Sets?","yn");
		}
		if((ch=='y') || (PAR_FOL&&SHOOT)) {
			oldt=DELTA_T;
			if(rp==1) {
				/* printf(" One real positive -- pos=%d lam=%g \n",pose,eval[2*pose]); */
				/*     for(i=0;i<n*n;i++)printf(" w=%g o=%g \n",work[i],oldwork[i]); */
				get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,eval[2*pose],ierr);
				if(*ierr==0) {
					change_current_linestyle(UnstableManifoldColor,&oldcol);
					pr_evec(x,b,n,pr,eval[2*pose],1);
					DELTA_T=fabs(DELTA_T);
					shoot(bp,x,b,1);
					shoot(bp,x,b,-1);
					change_current_linestyle(oldcol,&dummy);
				} else {
					err_msg("Failed to compute eigenvector");
				}
			}
			if(rn==1) {
				get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,eval[2*nege],ierr);
				if(*ierr==0) {
					change_current_linestyle(StableManifoldColor,&oldcol);
					pr_evec(x,b,n,pr,eval[2*nege],-1);
					DELTA_T=-fabs(DELTA_T);
					shoot(bp,x,b,1);
					shoot(bp,x,b,-1);
					change_current_linestyle(oldcol,&dummy);
				} else {
					err_msg("Failed to compute eigenvector");
				}
			}
			DELTA_T=oldt;
		}
	}  /* end of normal shooting stuff */

	/* strong (un) stable manifold calculation
	only one-d manifolds calculated */
	/* lets check to see if it is relevant */
	if(((rn>1) && (bneg>=0)) || ((rp>1) && (bpos>=0))) {
		ch='n';
		if(!PAR_FOL) {
			ch=(char)TwoChoice("YES","NO","Draw Strong Sets?","yn");
		}
		if((ch=='y')||(PAR_FOL&&SHOOT)) {
			oldt=DELTA_T;
			if((rp>1) && (bpos>=0)) {/* then there is a strong unstable */
				printf("strong unstable %g \n",bigpos);
				get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,bigpos,ierr);
				if(*ierr==0) {
					change_current_linestyle(UnstableManifoldColor,&oldcol);
					pr_evec(x,b,n,pr,bigpos,1);
					DELTA_T=fabs(DELTA_T);
					shoot(bp,x,b,1);
					shoot(bp,x,b,-1);
					change_current_linestyle(oldcol,&dummy);
				} else {
					err_msg("Failed to compute eigenvector");
				}
			}
			if((rn>1) && (bneg>=0)) {/* then there is a strong stable */
				printf("strong stable %g \n",bigneg);
				get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,bigneg,ierr);
				if(*ierr==0) {
					change_current_linestyle(StableManifoldColor,&oldcol);
					pr_evec(x,b,n,pr,bigneg,-1);
					DELTA_T=-fabs(DELTA_T);
					shoot(bp,x,b,1);
					shoot(bp,x,b,-1);
					change_current_linestyle(oldcol,&dummy);
				} else {
					err_msg("Failed to compute eigenvector");
				}
			}
		}
		DELTA_T=oldt;
	}
	free(work);
	return;
}


/* fixed point with no requests and store manifolds */
void do_sing_info(double *x, double eps, double err, double big,
				  int maxit, int n, double *er, double *em, int *ierr) {
	int kmem,i,j,ipivot[MAXODE];

	int rp=0,rn=0,cp=0,cn=0,im=0;
	int pose=0,nege=0,pr=0;
	double *work,*eval,*b,*bp,*oldwork,*ework;
	double temp,old_x[MAXODE];
	double real,imag;
	double bigpos=-1e10,bigneg=1e10;

	/* float xl[MAXODE]; */
	kmem=n*(2*n+5)+50;
	work=(double *)malloc(sizeof(double)*kmem);
	if(work==NULL) {
		/* err_msg("Insufficient core "); */
		return;
	}
	ShootICFlag=0;
	ShootIndex=0;
	for(i=0;i<n;i++) {
		old_x[i]=x[i];
	}
	oldwork=work+n*n;
	eval=oldwork+n*n;
	b=eval+2*n;
	bp=b+n;
	ework=bp+n;
	rooter(x,err,eps,big,work,ierr,maxit,n);
	if(*ierr!=0) {
		free(work);
		/* err_msg("Could not converge to root"); */
		for(i=0;i<n;i++) {
			x[i]=old_x[i];
		}
		return;
	}

	for(i=0;i<n*n;i++) {
		oldwork[i]=work[i];
	}
	/* Transpose for Eigen        */
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			temp=work[i+j*n];
			work[i+j*n]=work[i*n+j];
			work[i*n+j]=temp;
		}
	}
	eigen(n,work,eval,ework,ierr);
	if(*ierr!=0) {
		free(work);
		return;
	}
	/* succesfully computed evals now lets work with them */

	for(i=0;i<n;i++) {
		real=eval[2*i];
		imag=eval[2*i+1];
		er[i]=real;
		em[i]=imag;

		if(METHOD == METHOD_DISCRETE) {
			real=real*real+imag*imag-1.00;
		}
		if(fabs(imag)<.00000001) {
			imag=0.0;
		}
		if(real<0.0) {
			if(imag!=0.0) {
				cn++;
				if(real<bigneg) {
					bigneg=real;
				}
			} else {
				rn++;
				nege=i;
				if(real<bigneg) {
					bigneg=real;
				}
			}
		}
		if(real>0.0) {
			if(imag!=0.0) {
				cp++;
				if(real>bigpos) {
					bigpos=real;
				}
			} else {
				rp++;
				pose=i;
				if(real>bigpos) {
					bigpos=real;
				}
			}
		}
		if((real==0.0) && (imag!=0.0)) {
			im++;
		}
	}     /* eigenvalue count */
	if(((rp+cp)!=0) && ((rn+cn)!=0)) {
		eq_symb(x,1);
	} else {
		if((rp+cp)!=0) {
			eq_symb(x,0);
		} else {eq_symb(x,3);
		}
	}
	/* Lets change Work back to transposed oldwork */
	for(i=0;i<n;i++) {
		for(j=i+1;j<n;j++) {
			temp=oldwork[i+j*n];
			work[i+j*n]=oldwork[i*n+j];
			work[i*n+j]=temp;
		}
	}
	if((n>1)) {
		if(rp==1) {
			/* printf(" One real positive -- pos=%d lam=%g \n",pose,eval[2*pose]); */
			/*     for(i=0;i<n*n;i++)printf(" w=%g o=%g \n",work[i],oldwork[i]); */
			get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,eval[2*pose],ierr);

			if(*ierr==0) {
				pr_evec(x,b,n,pr,eval[2*pose],1);
			}
		}

		if(rn==1) {
			get_evec(work,oldwork,b,bp,n,maxit,err,ipivot,eval[2*nege],ierr);
			if(*ierr==0) {
				pr_evec(x,b,n,pr,eval[2*nege],-1);
			}
		}
	}

	free(work);
	return;
}


void rooter(double *x, double err, double eps, double big, double *work, int *ierr, int maxit, int n) {
	int i,iter,ipivot[MAXODE],info;
	char ch;
	double *xp,*yp,*y,*xg,*dermat,*dely;
	double r;
	dermat=work;
	xg=dermat+n*n;
	yp=xg+n;
	xp=yp+n;
	y=xp+n;
	dely=y+n;
	iter=0;
	*ierr=0;
	while(1) {
		ch=my_abort();

		if(ch==27) {
			*ierr=1;
			return;
		}
		if(ch=='/') {
			*ierr=1;
			ENDSING=1;
			return;
		}
		if(ch=='p') {
			PAUSER=1;
		}
		getjac(x,y,yp,xp,eps,dermat,n);
		sgefa(dermat,n,n,ipivot,&info);
		if(info!=-1) {
			*ierr=1;
			return;
		}
		for(i=0;i<n;i++) {
			dely[i]=y[i];
		}
		sgesl(dermat,n,n,ipivot,dely,0);
		r=0.0;
		for(i=0;i<n;i++) {
			x[i]=x[i]-dely[i];
			r=r+fabs(dely[i]);
		}
		if(r<err) {
			getjac(x,y,yp,xp,eps,dermat,n);
			if(METHOD == METHOD_DISCRETE) {
				for(i=0;i<n;i++) {
					dermat[i*(n+1)]+=1.0;
				}
			}
			return; /* success !! */
		}
		if((r/(double)n)>big) {
			*ierr=1;
			return;
		}
		iter++;
		if(iter>maxit) {
			*ierr=1;
			return;
		}
	}
}


/* --- Static functions --- */
static void pr_evec(double *x, double *ev, int n, int pr, double eval,int type) {
	int i;
	double d=fabs(DELTA_T)*.1;
	ShootICFlag=1;
	if(ShootIndex<7) {
		for(i=0;i<n;i++) {
			ShootIC[ShootIndex][i]=x[i]+d*ev[i];
			ShootType[ShootIndex]=type;
			ShootIC[ShootIndex+1][i]=x[i]-d*ev[i];
			ShootType[ShootIndex+1]=type;
		}
		ShootIndex+=2;
	}
	if(pr==0) {
		return;
	}
}
