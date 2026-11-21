#include "stdfix.h" 
#include "auto_f2c.h" 
#include "auto_nox.h"
#include "autlim.h"
#include "derived.h"
#include "pp_shoot.h"
#include "xAuto.h"
#include "delay_handle.h"
#include "del_stab.h"

void redo_all_fun_tables();
void getjactrans(double *x,double *y,double *yp,double *xp, double eps, double *d, int n);
extern XAUTO xAuto;

/*    Hooks to xpp RHS     */

extern int (*rhs)();
extern double constants[],last_ic[];

extern int Auto_index_to_array[8];
extern int NewPeriodFlag;
extern int AutoTwoParam,NAutoPar;
extern int HomoFlag;
extern double homo_l[100],homo_r[100];
extern int METHOD,NJMP;
extern double outperiod[];
extern int UzrPar[],NAutoUzr;

extern double NEWT_ERR;
extern double delay_list[];
extern int NDelay;

/* Eigenvalue replacement for steady states of delay equations.  AUTO's
   fnhbae() normally takes the eigenvalues of df/du, which for a DDE is the
   *sum* A_0+A_1+... and says nothing about stability.  Here we return the
   ndm rightmost roots of the characteristic equation
        det( lambda I - A_0 - sum_k A_k exp(-lambda tau_k) ) = 0
   at the steady state u for the parameters par.  Returns 0 on success.   */
int auto_delay_eigen(ndm,u,par,ev,nev)
integer ndm;
double *u,*par;
doublecomplex *ev;
integer *nev;
{
  int i,n=(int)ndm,nun=0,nref=0,nr,want;
  double x[NAUTO],evs[2*(NAUTO+64)],*coef=NULL,cn;
  if(n>NAUTO)return 1;
  for(i=0;i<NAutoPar;i++)
    constants[Auto_index_to_array[i]]=par[i];
  evaluate_derived();
  redo_all_fun_tables();
  for(i=0;i<n;i++)x[i]=u[i];
  cn=delay_char_coef(x,n,NEWT_ERR,&coef);
  if(cn<0.0||coef==NULL)return 1;
  /* look at more roots than there are equations: for a DDE the spectrum is
     infinite and a pair that crosses the imaginary axis need not be among
     the n rightmost ones */
  want=2*n;
  if(want<8)want=8;
  if(want>40)want=40;
  if(want<n)want=n;
  nr=delay_spectrum(coef,delay_list,n,NDelay,want,evs,&nun,&nref,200);
  free(coef);
  if(nr<n)return 1;
  for(i=0;i<nr;i++){
    ev[i].r=evs[2*i];
    ev[i].i=evs[2*i+1];
  }
  *nev=nr;
  return 0;
}

/* Real and imaginary parts of the characteristic matrix of a delay equation
   at a steady state u, for the 2-parameter Hopf continuation.  With
   Delta(i w) = i w I - A_0 - sum_k A_k exp(-i w tau_k) one has
        A_0 + sum_k A_k exp(-i w tau_k) = Cr - i S ,
        Cr = A_0 + sum_k A_k cos(w tau_k),   S = sum_k A_k sin(w tau_k).
   Both are returned in column-major layout M[i+j*n] (row i, column j) like
   dfdu in AUTO.  Parameters must already be in `constants' (funi does that).
   Returns 0 on success.                                                   */
int auto_delay_hopf_matrices(ndm,u,omega,Cr,S)
integer ndm;
double *u,omega,*Cr,*S;
{
  int i,j,k,n=(int)ndm;
  double x[NAUTO],*coef=NULL,cn,c,s;
  if(n>NAUTO)return 1;
  for(i=0;i<n;i++)x[i]=u[i];
  cn=delay_char_coef(x,n,NEWT_ERR,&coef);
  if(cn<0.0||coef==NULL)return 1;
  for(i=0;i<n;i++)for(j=0;j<n;j++){
    Cr[i+j*n]=coef[i*n+j];
    S[i+j*n]=0.0;
  }
  for(k=1;k<=NDelay;k++){
    c=cos(omega*delay_list[k-1]);
    s=sin(omega*delay_list[k-1]);
    for(i=0;i<n;i++)for(j=0;j<n;j++){
      Cr[i+j*n]+=coef[i*n+j+n*n*k]*c;
      S[i+j*n]+=coef[i*n+j+n*n*k]*s;
    }
  }
  free(coef);
  return 0;
}

int func(ndim, u, icp, par, ijac, f, dfdu, dfdp)
integer ndim,*icp,ijac;
double  *u,*par,*f,*dfdu,*dfdp;
{
   int i,j,delay_old;
   double zz[NAUTO];
   double  y[NAUTO],yp[NAUTO],xp[NAUTO];
   for(i=0;i<NAutoPar;i++){
     constants[Auto_index_to_array[i]]=par[i];
     
   }
   evaluate_derived();
   redo_all_fun_tables();
   /* For a delay equation AUTO works on the steady-state problem
      F(x,x,...,x;p)=0, so every delayed variable must equal the
      current one rather than a value from the integration history. */
   delay_old=delay_enter_steady();
   rhs(0.0,u,f,ndim);
   if(ijac==1){
     getjactrans(u,y,yp,xp,NEWT_ERR,dfdu,ndim);
   }
   if(METHOD>0||NJMP==1){delay_leave_steady(delay_old);return 0;}
   for(i=1;i<NJMP;i++){
     for(j=0;j<ndim;j++)
       zz[j]=f[j];
     rhs(0.0,zz,f,ndim);
   }
   delay_leave_steady(delay_old);

   return 0;





} /* func_ */


int stpnt(ndim,t,u,par)
integer ndim;
doublereal *u, *par,t;
{
  int i;

  double p;

  for(i=0;i<NAutoPar;i++)
    par[i] = constants[Auto_index_to_array[i]];

  if(NewPeriodFlag==0){  
    for(i=0;i<ndim;i++)
      u[i]=last_ic[i];
    return 0;
  }

  get_start_period(&p);
  par[10]=p;
  if(HomoFlag!=1)get_start_orbit(u,t,p,ndim);
  /*  printf("%d %d %g %g %g %g \n",ndim,HomoFlag,t,u[0],u[1],p); */
  if(HomoFlag==1){

    get_shifted_orbit(u,t,p,ndim);
    for(i=0;i<ndim;i++){
      par[11+i]=homo_l[i];

    }
  }
  if(HomoFlag==2){ /* heteroclinic */
    for(i=0;i<ndim;i++){
      par[11+i]=homo_l[i];
      par[11+i+ndim]=homo_r[i];

    }

  }
  return 0;

} /* stpnt_ */







/* Subroutine */ int bcnd(ndim, par, icp, nbc, u0, u1, ijac, fb, dbc)
integer ndim;
double *par;
integer *icp, nbc;
double *u0, *u1, *fb;
integer ijac;
double *dbc;
{
 int i;
/* Hooks to the XPP bc parser!! */

 for(i=0;i<NAutoPar;i++){
     constants[Auto_index_to_array[i]]=par[i];
 }

 
 evaluate_derived();
 redo_all_fun_tables();
 do_bc(u0,0.0,u1,1.0,fb,nbc);

    return 0;
} /* bcnd_ */

/* Subroutine */ int icnd(ndim, par, icp, nint, u, uold, udot, upold, fi, 
	ijac, dint)
integer *ndim;
double *par;
integer *icp, *nint;
double *u, *uold, *udot, *upold, *fi;
integer *ijac;
double *dint;
{
   int i;
   double dum=0.0;
   /*
  for(i=0;i<Homo_n;i++)
    dum+=upold[i]*(u[i]-uold[i]);
  fi[0]=dum;
   */
    return 0;
} /* icnd_ */

/* Subroutine */ int fopt(ndim, u, icp, par, ijac, fs, dfdu, dfdp)
integer *ndim;
double *u;
integer *icp;
double *par;
integer *ijac;
double *fs, *dfdu, *dfdp;
{
/*     ---------- ---- */
    return 0;
} /* fopt_ */

/*  Not sure what to do here; I think  do nothing  since IEQUIB is always
    -2 
*/
int pvls (integer ndim, const doublereal *u,
          doublereal *par)
{
  return 0;
}


