#include <stdlib.h>

#include "autlim.h"
#include "auto/auto_f2c.h"
#include "auto_nox.h"
#include "derived.h"
#include "solver/gear.h"
#include "load_eqn.h"
#include "main.h"
#include "my_rhs.h"
#include "numerics.h"
#include "parserslow.h"
#include "pp_shoot.h"
#include "util/matrixalgebra.h"
#include "xAuto.h"

XAUTO xAuto;

/*    Hooks to xpp RHS     */
int func(integer ndim, double *u, integer *icp, double *par, integer ijac,
		 double *f, double *dfdu, double *dfdp) {
	int i,j;
	double zz[NAUTO];
	double  y[NAUTO],yp[NAUTO],xp[NAUTO];
	for(i=0;i<NAutoPar;i++) {
		constants[Auto_index_to_array[i]]=par[i];
	}
	evaluate_derived();
	my_rhs(0.0,u,f,ndim);
	if(ijac==1) {
		getjactrans(u,y,yp,xp,NEWT_ERR,dfdu,ndim);
	}
	if(METHOD != METHOD_DISCRETE || NJMP==1) {
		return 0;
	}

	for(i=1;i<NJMP;i++) {
		for(j=0;j<ndim;j++) {
			zz[j]=f[j];
		}
		my_rhs(0.0,zz,f,ndim);
	}
	return 0;
}


int stpnt(integer ndim, doublereal t, doublereal *u, doublereal *par) {
	int i;
	double p;

	for(i=0;i<NAutoPar;i++) {
		par[i] = constants[Auto_index_to_array[i]];
	}
	if(NewPeriodFlag==0) {
		for(i=0;i<ndim;i++) {
			u[i]=last_ic[i];
		}
		return 0;
	}

	get_start_period(&p);
	par[10]=p;
	if(HomoFlag!=1) {
		get_start_orbit(u,t,p,ndim);
	}
	if(HomoFlag==1) {
		get_shifted_orbit(u,t,p,ndim);
		for(i=0;i<ndim;i++) {
			par[11+i]=homo_l[i];
		}
	}
	if(HomoFlag==2) { /* heteroclinic */
		for(i=0;i<ndim;i++) {
			par[11+i]=homo_l[i];
			par[11+i+ndim]=homo_r[i];
		}
	}
	return 0;
} /* stpnt_ */


/* Subroutine */
int bcnd(integer ndim, double *par, integer *icp, integer nbc, double *u0, double *u1,
		 integer ijac, double *fb, double *dbc) {
	int i;
	/* Hooks to the XPP bc parser!! */
	for(i=0;i<NAutoPar;i++) {
		constants[Auto_index_to_array[i]]=par[i];
	}
	evaluate_derived();
	do_bc(u0,0.0,u1,1.0,fb,nbc);
	return 0;
} /* bcnd_ */


/* Subroutine */
int icnd(integer *ndim, double *par, integer *icp, integer *nint,
		 double *u, double *uold, double *udot, double *upold, double *fi,
		 integer *ijac, double *dint) {
	/*
	int i;
	double dum=0.0;
	for(i=0;i<Homo_n;i++)
	dum+=upold[i]*(u[i]-uold[i]);
  fi[0]=dum;
   */
	return 0;
} /* icnd_ */


/* Subroutine */
int fopt(integer *ndim, double *u, integer *icp, double *par,
		 integer *ijac, double *fs, double *dfdu, double *dfdp) {
	/*     ---------- ---- */
	return 0;
} /* fopt_ */

/*  Not sure what to do here; I think  do nothing  since IEQUIB is always
	-2
*/
int pvls (integer ndim, const doublereal *u, doublereal *par) {
	return 0;
}


