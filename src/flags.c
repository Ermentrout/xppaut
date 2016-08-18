/*  this is a new (Summer 1995) addition to XPP that allows one to
	do things like delta functions and other discontinuous
	stuff.

	The conditions are set up as part of the "ODE" file:

global sign condition {event1;....;eventn}
global sign {condition} {event1;....;eventn}

the {} and ;  are required for the events

condition is anything that when it evaluates to 0 means the flag should be
set.  The sign is like in Poincare maps, thus let C(t1) and C(t2) be
the value of the condition at t1  and t2.
sign = 0 ==>  just find when C(t)=0
sign = 1 ==>  C(t1)<0 C(t2)>0
sign = -1==>  C(t1)>0 C(t2)<0

To get the time of the event, we use linear interpolation:

 t* = t1 + (t2-t1)
		   -------   (0-C(t1))
		  C(t2)-C(t1)
This yields  the variables, etc at that time

Now what are the events:

They are of the form:
   name = expression  the variable  <name> is replaced by the value of
  <expression>

Note that there may be several "conditions" defined and that
these must also be checked to see if they have been switched
and in what order.  This is particularly true for "delta" function
type things.


Here is a simple example -- the kicked cycle:
dx/dt = y
dy/dy = -x -c y

if y=0 and y goes from pos to neg then x=x+b
here is how it would work:

global -1 y {x=x+b}


Here is Tysons model:

global -1 u-.2 {m=.5*m}

*/

/*

type =0 variable
type =1 parameter
type =2 output
type =3 halt
*/

#include "flags.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "derived.h"
#include "form_ode.h"
#include "ggets.h"
#include "init_conds.h"
#include "integrate.h"
#include "newpars.h"
#include "parserslow.h"
#include "pop_list.h"
#include "solver/cv2.h"
#include "solver/dormpri.h"
#include "solver/gear.h"
#include "solver/stiff.h"


/* --- Macros --- */
#define MY_DBL_EPS 5e-16
#define MAX_EVENTS 20 /*  this is the maximum number of events per flag */


/* --- Types --- */
typedef struct {
	double f0,f1;
	double tstar;
	int lhs[MAX_EVENTS];
	double vrhs[MAX_EVENTS];
	char lhsname[MAX_EVENTS][11];
	char *rhs[MAX_EVENTS];
	int *comrhs[MAX_EVENTS];
	char *cond;
	int *comcond;
	int sign,nevents;
	int hit,type[MAX_EVENTS];
	int anypars;
	int nointerp;
} FLAG;


/* --- Data --- */
static FLAG flag[MAXFLAG];
int NFlags=0;
double STOL=1.e-10;


/* --- Functions --- */
int add_global(char *cond, int sign, char *rest) {
	char temp[DEFAULT_STRING_LENGTH];
	int nevents,ii,k,l,lt,j=NFlags;
	char ch;
	if(NFlags>=MAXFLAG) {
		plintf("Too many global conditions\n");
		return(1);
	}
	l=strlen(cond);
	flag[j].cond=(char *) malloc(l+1);
	strcpy(flag[j].cond,cond);
	nevents=0;
	flag[j].lhsname[0][0]=0;
	k=0;
	l=strlen(rest);
	for(ii=0;ii<l;ii++) {
		ch=rest[ii];
		if(ch=='{' || ch==' ') {
			continue;
		}
		if(ch=='}' || ch==';') {
			if(nevents==MAX_EVENTS) {
				printf(" Too many events per flag \n");
				return(1);
			}
			temp[k]=0;
			lt=strlen(temp);
			if(flag[j].lhsname[nevents][0]==0) {
				printf(" No event variable named for %s \n",temp);
				return(1);
			}
			flag[j].rhs[nevents]=(char *)malloc(lt+1);
			strcpy(flag[j].rhs[nevents],temp);
			nevents++;
			k=0;
			if(ch=='}') {
				break;
			}
			continue;
		}
		if(ch=='=') {
			temp[k]=0;
			strcpy(flag[j].lhsname[nevents],temp);

			k=0;
			if(nevents<MAX_EVENTS-1) {
				flag[j].lhsname[nevents+1][0]=0;
			}
			continue;
		}

		temp[k]=ch;
		k++;
	}
	if(nevents==0) {
		plintf(" No events for condition %s \n",cond);
		return(1);
	}
	/*  we now have the condition, the names, and the formulae */
	flag[j].sign=sign;
	flag[j].nevents=nevents;
	NFlags++;
	return(0);
}


void show_flags(void) {
	/* uncomment for debugging */
	/*
 for(i=0;i<NFlags;i++) {
   n=flag[i].nevents;
   plintf(" Flag %d has sign %d and %d events and condition %s \n",
	  i+1,flag[i].sign,n,flag[i].cond);
   for(j=0;j<n;j++)
	 plintf("%d:  %s [%d] = %s \n",j+1,flag[i].lhsname[j],flag[i].lhs[j],
		flag[i].rhs[j]);
 }
 */
}

int compile_flags(void) {
	int j;
	int i,k,index,nc;
	int command[DEFAULT_STRING_LENGTH];
	if(NFlags==0) {
		return 0;
	}
	for(j=0;j<NFlags;j++) {
		if(add_expr(flag[j].cond,command,&nc)) {
			plintf("Illegal global condition:  %s\n",flag[j].cond);
			return(1);
		}
		flag[j].anypars=0;
		flag[j].nointerp=0;
		flag[j].comcond=(int *)malloc(sizeof(int)*(nc+1));
		for(k=0;k<=nc;k++) {
			flag[j].comcond[k]=command[k];
		}
		for(i=0;i<flag[j].nevents;i++) {
			index=find_user_name(INITIAL_CONDITION_BOX,flag[j].lhsname[i]);
			if(index<0) {
				index=find_user_name(PARAMETER_BOX,flag[j].lhsname[i]);
				if(index<0) {
					if(strcasecmp(flag[j].lhsname[i],"out_put")==0)	{
						flag[j].type[i]=2;
						flag[j].lhs[i]=0;
					} else if(strcasecmp(flag[j].lhsname[i],"arret")==0) {
						flag[j].type[i]=3;
						flag[j].lhs[i]=0;
					} else if(strcasecmp(flag[j].lhsname[i],"no_interp")==0) {
						flag[j].nointerp=1;
						flag[j].type[i]=0;
						flag[j].lhs[i]=0;
					} else {
						plintf(" <%s> is not a valid variable/parameter name \n",
							   flag[j].lhsname[i]);
						return(1);
					}
				} else {
					flag[j].lhs[i]=index;
					flag[j].type[i]=1;
					flag[j].anypars=1;
				}
			} else {
				flag[j].lhs[i]=index;
				flag[j].type[i]=0;
			}
			if(add_expr(flag[j].rhs[i],command,&nc)) {
				printf("Illegal event %s for global %s\n",
					   flag[j].rhs[i],flag[j].cond);
				return(1);
			}
			flag[j].comrhs[i]=(int *)malloc(sizeof(int)*(nc+1));
			for(k=0;k<=nc;k++) {
				flag[j].comrhs[i][k]=command[k];
			}
		}
	}
	return(0);
}



/*  here is the shell code for a loop around  integration step  */
int one_flag_step(double *yold, double *ynew, int *istart, double told, double *tnew, int neq, double *s) {
	double dt=*tnew-told;
	double f0,f1,tol,tolmin=1e-10;
	double smin=2;
	int sign,i,j,in,ncycle=0,newhit,nevents;

	if(NFlags==0) {
		return(0);
	}
	for(i=0;i<NFlags;i++) {
		flag[i].tstar=2.0;
		flag[i].hit=0;
	}
	/* If this is the first call, then need f1  */
	if(*istart==1) {
		for(i=0;i<neq;i++) {
			SETVAR(i+1,yold[i]);
		}
		SETVAR(0,told);
		for(i=0;i<NFlags;i++) {
			*istart=0;
		}
	}
	for(i=0;i<NFlags;i++) {
		sign=flag[i].sign;
		flag[i].f0=flag[i].f1;
		f0=flag[i].f0;
		for(j=0;j<neq;j++) {
			SETVAR(j+1,ynew[j]);
		}
		SETVAR(0,*tnew);
		f1=evaluate(flag[i].comcond);
		flag[i].f1=f1;
		tol=fabs(f1-f0);
		/* plintf(" call1 %g %g %g %g\n",told,f0,f1,smin);  */
		switch(sign) {
		case 1:
			if((((f0<0.0) && (f1>0.0)) || ((f0<0.0) && (f1>0.0))) && tol>tolmin) {
				flag[i].hit=ncycle+1;
				flag[i].tstar=f0/(f0-f1);
				/* plintf(" f0=%g, f1=%g tstar=%g at t=%g\n tol=%g",f0,f1,flag[i].tstar,*tnew,tol);  */ /* COMMENT! */
			}
			break;
		case -1:
			if(f0>0.0 && f1<=0.0 && tol>tolmin) {
				flag[i].hit=ncycle+1;
				flag[i].tstar=f0/(f0-f1);
			}
			break;
		case 0:
			if(fabs(f1)<MY_DBL_EPS) {
				flag[i].hit=ncycle+1;
				flag[i].tstar=told;
			}
			break;
		}
		if(flag[i].nointerp==1) {
			flag[i].tstar=1.0;
		}

		if(smin>flag[i].tstar) {
			smin=flag[i].tstar;
		}
	} /* run through flags */

	if(smin<STOL) {
		smin=STOL;
	} else {
		smin=(1+STOL)*smin;
	}
	if(smin>1.0) {
		return(0);
	}
	*tnew=told+dt*smin;
	SETVAR(0,*tnew);
	for(i=0;i<neq;i++) {
		ynew[i]=yold[i]+smin*(ynew[i]-yold[i]);
		SETVAR(i+1,ynew[i]);
	}
	for(i=0;i<NFlags;i++) {
		flag[i].f0=evaluate(flag[i].comcond);
	}
	while(1) { /* run through all possible events  */
		ncycle++;
		newhit=0;
		for(i=0;i<NFlags;i++) {
			nevents=flag[i].nevents;
			if(flag[i].hit==ncycle && flag[i].tstar<=smin) {
				for(j=0;j<nevents;j++) {
					flag[i].vrhs[j]=evaluate(flag[i].comrhs[j]);
					in=flag[i].lhs[j];
					if(flag[i].type[j]==0) {
						SETVAR(in+1,flag[i].vrhs[j]);
					}
				}
			}
		}
		for(i=0;i<NFlags;i++) {
			nevents=flag[i].nevents;
			if(flag[i].hit==ncycle && flag[i].tstar<=smin) {
				for(j=0;j<nevents;j++) {
					in=flag[i].lhs[j];
					if(flag[i].type[j]==0) {
						ynew[in]=flag[i].vrhs[j];
					} else if(flag[i].type[j]==1) {
						set_val(upar_names[in],flag[i].vrhs[j]);
					} else {
						if((flag[i].type[j]==2) && (flag[i].vrhs[j]>0)) {
							send_output(ynew,*tnew);
						}
						if((flag[i].type[j]==3) && (flag[i].vrhs[j]>0)) {
							send_halt(ynew,*tnew);
						}
					}
				}
				if(flag[i].anypars) {
					evaluate_derived();
					redraw_params();
				}
			}
		}
		for(i=0;i<neq;i++) {
			SETVAR(i+1,ynew[i]);
		}
		for(i=0;i<NFlags;i++) {
			flag[i].f1=evaluate(flag[i].comcond);
			if(flag[i].hit>0) {
				continue; /* already hit so dont do anything */
			}
			f1=flag[i].f1;
			sign=flag[i].sign;
			f0=flag[i].f0;
			tol=fabs(f1-f0);
			switch(sign) {
			case 1:
				if(f0<=0.0 && f1>=0.0 && tol>tolmin) {
					flag[i].tstar=smin;
					flag[i].hit=ncycle+1;
					newhit=1;
				}
				break;
			case -1:
				if(f0>=0.0 && f1<=0.0 && tol>tolmin) {
					flag[i].tstar=smin;
					flag[i].hit=ncycle+1;
					newhit=1;
				}
				break;
			case 0:
				if(f0*f1<=0 && (f1!=0 || f0!=0) && tol>tolmin) {
					flag[i].tstar=smin;
					flag[i].hit=ncycle+1;
					newhit=1;
				}
			}
		}
		if(newhit==0) {
			break;
		}
	}
	*s=smin;
	return(1);
}
