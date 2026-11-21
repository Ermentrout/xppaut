#include "stdfix.h"

#include "eig_list.h" 
#include "gear.h"
#include "ggets.h"

#include <math.h>
#include <stdio.h>
#include "xpplim.h"
#include "del_stab.h"

#define Z(a,b) z[(a)+n*(b)]
#define DING ping()
/* this code takes the determinant of a complex valued matrix
*/
int TwoChoice(char *c1,char *c2, char *q,char *key);
extern double variable_shift[2][MAXODE],AlphaMax,OmegaMax;

extern double delay_list[MAXDELAY];
extern int NDelay,del_stab_flag,WhichDelay,DelayGrid,PAR_FOL;
extern int (*rhs)();
void eq_symb(double *x,int flag);
double amax();

/*typedef struct{
  double r,i;
}COMPLEX;
*/

/* The
 code here replaces the do_sing code if the equation is
   a delay differential equation. 
*/

/* Build the coefficient matrices of the linearization about the
   steady state x of
        x'(t) = f( x(t), x(t-tau_1), ..., x(t-tau_m) ) :
   A_0 = df/dx(t), A_k = df/dx(t-tau_k).  The distinct delays found at x are
   left in delay_list[0..NDelay-1].  *coef_out is malloc'ed (caller frees)
   with layout coef[j*n+i+n*n*k] = (A_k)_{j,i}.  Central differences are
   used.  Returns an upper bound for |lambda| of any root with Re>=0
   (sum of the matrix 1-norms), or <0 on failure.

   Note: at a steady state a state-dependent delay tau(x) can be treated as
   the constant tau(x*): the extra term in the linearization is proportional
   to x'(t-tau)=0.                                                        */
double delay_char_coef(double *x,int n,double eps,double **coef_out)
{
  double *coef,yp[MAXODE],ym[MAXODE],y[MAXODE],xp[MAXODE],dx;
  double colmax,colsum,colnorm;
  int i,j,k,oldflag=del_stab_flag;
  double h0=1e-6;     /* central-difference step (relative) */

  del_stab_flag=0;     /* scan: collect the distinct delays */
  NDelay=0;
  rhs(0.0,x,y,n);
  for(i=0;i<n;i++){
    variable_shift[0][i]=x[i];
    variable_shift[1][i]=x[i];
  }
  coef=(double *)malloc(sizeof(double)*n*n*(NDelay+1));
  if(coef==NULL){del_stab_flag=oldflag;return -1.0;}

  del_stab_flag=-1;
  WhichDelay=-1;
  colnorm=0.0;
  colmax=0.0;
  for(i=0;i<n;i++){
    colsum=0.0;
    dx=h0*amax(1.0,fabs(x[i]));
    for(j=0;j<n;j++)xp[j]=x[j];
    xp[i]=x[i]+dx;
    rhs(0.0,xp,yp,n);
    xp[i]=x[i]-dx;
    rhs(0.0,xp,ym,n);
    for(j=0;j<n;j++){
      coef[j*n+i]=(yp[j]-ym[j])/(2.0*dx);
      colsum+=fabs(coef[j*n+i]);
    }
    if(colsum>colmax)colmax=colsum;
  }
  colnorm=colmax;
  for(k=0;k<NDelay;k++){
    WhichDelay=k;
    colmax=0.0;
    for(i=0;i<n;i++){
      colsum=0.0;
      dx=h0*amax(1.0,fabs(x[i]));
      for(j=0;j<n;j++)variable_shift[1][j]=variable_shift[0][j];
      variable_shift[1][i]=x[i]+dx;
      rhs(0.0,x,yp,n);
      variable_shift[1][i]=x[i]-dx;
      rhs(0.0,x,ym,n);
      variable_shift[1][i]=x[i];
      for(j=0;j<n;j++){
        coef[j*n+i+n*n*(k+1)]=(yp[j]-ym[j])/(2.0*dx);
        colsum+=fabs(coef[j*n+i+n*n*(k+1)]);
      }
      if(colsum>colmax)colmax=colsum;
    }
    colnorm+=colmax;
  }
  WhichDelay=-1;
  del_stab_flag=oldflag;
  *coef_out=coef;
  return colnorm;
}

void do_delay_sing(x,eps,err,big,maxit,n,ierr,stabinfo)
     double *x,eps,err,big;
     int *ierr,n,maxit;
     float *stabinfo;
{
      double rr[2];

 double colnorm=0,colmax,colsum;
 double *work,old_x[MAXODE],sign;
 double *coef,yp[MAXODE],y[MAXODE],xp[MAXODE],dx;
 int kmem=n*(2*n+5)+50,i,j,k,okroot;

 double *ev;
 ev=(double *)malloc(2*n*sizeof(double));
 for(i=0;i<(2*n);i++)ev[i]=0.0;
 /* first we establish how many delays there are */
 del_stab_flag=0;
 for(i=0;i<n;i++)old_x[i]=x[i];
 work=(double *)malloc(kmem*sizeof(double));
 rooter(x,err,eps,big,work,ierr,maxit,n);
 if(*ierr!=0)
   {
     del_stab_flag=1;
     free(work);
     err_msg("Could not converge to root");
     for(i=0;i<n;i++)x[i]=old_x[i];
     return;
   }
 /* OKAY -- we have the root */
 free(work);
 colnorm=delay_char_coef(x,n,eps,&coef);
 if(colnorm<0.0){
   del_stab_flag=1;
   free(ev);
   err_msg("Out of memory in delay stability");
   return;
 }
 del_stab_flag=-1;
 /* plintf("Norm= %g \n",colnorm); */
 /* Proper spectrum: Chebyshev collocation of the generator, polished by
    Newton on the characteristic equation (del_spec.c).  The old contour
    count below is kept as an independent cross-check.                    */
 {
   int nun=0,nref=0,nr,cp=0,cn=0,rp=0,rn=0,im=0,ii,tot;
   int nreq=(n>12)?n:12;          /* roots to compute: at least n, and 12 */
   double *evs;
   double re,imag;
   char ch='n';
   if(nreq>40&&n<=40)nreq=40;
   evs=(double *)malloc(2*nreq*sizeof(double));
   if(evs!=NULL)for(ii=0;ii<2*nreq;ii++)evs[ii]=0.0;
   nr=(evs==NULL)?-1:delay_spectrum(coef,delay_list,n,NDelay,nreq,evs,&nun,&nref,600);
   if(nr>0){
     for(ii=nr;ii<n;ii++){evs[2*ii]=evs[2*(nr-1)];evs[2*ii+1]=evs[2*(nr-1)+1];}
     for(ii=0;ii<n;ii++){
       re=evs[2*ii];imag=evs[2*ii+1];
       if(fabs(imag)<1e-8)imag=0.0;
       if(fabs(re)<1e-9*(1.0+fabs(imag)))re=0.0;
       if(re<0.0){ if(imag!=0.0)cn++; else rn++; }
       else if(re>0.0){ if(imag!=0.0)cp++; else rp++; }
       else if(imag!=0.0)im++;
     }
     tot=cp+rp;
     if(nun>tot)cp+=nun-tot;   /* more unstable roots than slots shown */
     sign=plot_args(coef,delay_list,n,NDelay,DelayGrid,colnorm,colnorm);
     if(abs((int)sign)!=nun)
       plintf(" Note: %d unstable roots from Newton/collocation but argument-principle count is %d -- check delays/DelayGrid\n",nun,abs((int)sign));
     if(!PAR_FOL)
       ch=(char)TwoChoice("YES","NO","Print eigenvalues?","yn");
     if(ch=='y'){
       plintf("\n Characteristic roots (delay equation) -- the %d rightmost:\n",nr);
       if(NDelay>0){
         plintf(" Delays at this equilibrium:");
         for(ii=0;ii<NDelay;ii++)plintf(" %g",delay_list[ii]);
         plintf("\n");
       }
       for(ii=0;ii<nr;ii++)
         plintf(" %3d:  %.10g  %c  i  %.10g \n",ii+1,evs[2*ii],
                (evs[2*ii+1]<0.0)?'-':'+',fabs(evs[2*ii+1]));
       plintf(" Roots with Re>0: %d     (argument-principle check: %d)\n",nun,abs((int)sign));
     }
     if(cp+rp!=0&&(rn+cn)!=0)eq_symb(x,1);
     else if(cp+rp!=0)eq_symb(x,0);
     else eq_symb(x,3);
     *stabinfo=(float)(cp+rp)+(float)(cn+rn)/1000.0;
     create_eq_box(cp,cn,rp,rn,im,x,evs,n);
     free(evs);
     free(coef);
     free(ev);
     del_stab_flag=1;
     return;
   }
   free(evs);
 }
 /* sign=plot_args(coef,delay_list,n,NDelay,DelayGrid,AlphaMax,OmegaMax); */
 sign=plot_args(coef,delay_list,n,NDelay,DelayGrid,colnorm,colnorm);

 okroot=find_positive_root(coef,delay_list,n,NDelay,colnorm,err,eps,big,maxit,rr);
 if(okroot>0){
   ev[0]=rr[0];
   ev[1]=rr[1];
 }
 free(coef);  
 *stabinfo=(float)fabs(sign);
 /* if(*stabinfo>0) */
 i=(int)sign;
if(i==0&&okroot==1&&AlphaMax>0)
  i=2;

 create_eq_box(abs(i),2,0,0,0,x,ev,n);
 /* DING; */
 del_stab_flag=1;
 free(ev);
 if(okroot==1)*stabinfo=AlphaMax;
}





COMPLEX csum(z,w)
     COMPLEX z,w;
{
  COMPLEX sum;
  sum.r=z.r+w.r;
  sum.i=z.i+w.i;
  return sum;
}

COMPLEX cdif(z,w)
     COMPLEX z,w;
{
   COMPLEX sum;
  sum.r=z.r-w.r;
  sum.i=z.i-w.i;
  return sum;
 }

COMPLEX cmlt(z,w)
     COMPLEX z,w;
{
   COMPLEX sum;
  sum.r=z.r*w.r-z.i*w.i;
  sum.i=z.r*w.i+z.i*w.r;
  return sum;
}

COMPLEX cdivv(z,w)
     COMPLEX z,w;
{
  COMPLEX sum;
  double amp=w.r*w.r+w.i*w.i;
  sum.r=(z.r*w.r+z.i*w.i)/amp;
  sum.i=(z.i*w.r-z.r*w.i)/amp;
  return sum;
}

COMPLEX cexp2(z)
     COMPLEX z;
{
  COMPLEX sum;
  double ex=exp(z.r);
  sum.r=ex*cos(z.i);
  sum.i=ex*sin(z.i);
  return sum;
}

void switch_rows(z,i1,i2,n)
     COMPLEX *z;
     int i1,i2,n;
{
  COMPLEX zt;
  int j;
  for(j=0;j<n;j++){
    zt=Z(i1,j);
    Z(i1,j)=Z(i2,j);
    Z(i2,j)=zt;
  }
}

COMPLEX rtoc(x,y)
     double x,y;
{
  COMPLEX sum;
  sum.i=y;
  sum.r=x;
  return sum;
}

void cprintn(z)
     COMPLEX z;
{
  plintf(" %g + i %g \n",z.r,z.i);
}

void cprint(z)
     COMPLEX z;
{
printf("(%g,%g) ",z.r,z.i);
}

void cprintarr(z,n,m)
     COMPLEX *z;
     int n,m;
{
  int i,j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++)
      cprint(z[i+m*j]);
    plintf("\n");
  }
}

double c_abs(z)
     COMPLEX z;
{
 return(sqrt(z.i*z.i+z.r*z.r));
}

COMPLEX cdeterm(z,n)
     COMPLEX *z;
     int n;
{
  int i,j,imax=0,k;
  double q,qmax;
  COMPLEX sign=rtoc(1.0,0.0),mult,sum,zd;
  for(j=0;j<n;j++){
    qmax=0.0;
    for(i=j;i<n;i++){
      q=c_abs(Z(i,j));
      if(q>qmax){
	qmax=q;
	imax=i;
      }
    }
    if(qmax==0.0)return(rtoc(0.0,0.0));
    switch_rows(z,imax,j,n);
    if(imax>j)sign=cmlt(rtoc(-1.0,0.0),sign);
    zd=Z(j,j);
    for(i=j+1;i<n;i++){
      mult=cdivv(Z(i,j),zd);
      for(k=j+1;k<n;k++){
	Z(i,k)=cdif(Z(i,k),cmlt(mult,Z(j,k)));
      }
    }
  }
  sum=sign;
  for(j=0;j<n;j++)
    sum=cmlt(sum,Z(j,j));
  return sum;
}
COMPLEX cxdeterm(z,n)
     COMPLEX *z;
     int n;
{
  int i,j,k;
  COMPLEX ajj,sum,mult;
  for(j=0;j<n;j++){
    ajj=Z(j,j);
    for(i=j+1;i<n;i++){
      mult=cdivv(Z(i,j),ajj);
      for(k=j+1;k<n;k++){
        Z(i,k)=cdif(Z(i,k),cmlt(mult,Z(j,k)));
      }
    }
  }
 /* now it should be diagonalized */
  sum=rtoc(1.0,0.0);
  for(j=0;j<n;j++){
    sum=cmlt(sum,Z(j,j));
  }
  return sum;
}
	 
void make_z(z,delay,n,m,coef,lambda)
     COMPLEX lambda;
     COMPLEX *z;
     double *coef,*delay;
     int n,m;
{
  int i,j,k,km;
  COMPLEX temp,eld;
  
 for(j=0;j<n;j++)
    for(i=0;i<n;i++){
      if(i==j)temp=lambda;
      else temp=rtoc(0.0,0.0);
      /* cprintn(temp); */
      z[i+j*n]=cdif(temp,rtoc(coef[i+j*n],0.0)); /* initialize the array */
    }
  for(k=0;k<m;k++){
    km=(k+1)*n*n;
    temp=rtoc(-delay[k],0.0); /* convert delay to complex number */
    eld=cexp2(cmlt(temp,lambda)); /* compute exp(-lambda*tau) */
    /* cprintn(eld); */
    for(j=0;j<n;j++)
      for(i=0;i<n;i++)
	z[i+j*n]=cdif(z[i+j*n],cmlt(eld,rtoc(coef[km+i+n*j],0.0)));
  }
}

int find_positive_root(coef,delay,n,m,rad,err,eps,big,maxit,rr)
     double *coef,*delay,*rr;
     int n,m,maxit;
     double rad;
     double eps,err,big;
{
  COMPLEX lambda,lambdap;
  COMPLEX det,*z,detp;
  double jac[4];
  double xl,yl,r,xlp,ylp;

  int k;
  
    lambda.r=AlphaMax;
    lambda.i=OmegaMax;
 
   z=(COMPLEX *)malloc(sizeof(COMPLEX)*n*n); 
 
  /* now Newtons Method for maxit times */
  for(k=0;k<maxit;k++){
      
    make_z(z,delay,n,m,coef,lambda);
    det=cdeterm(z,n);

    r=c_abs(det);
    if(r<err){ /* within the tolerance */
      process_root(lambda.r,lambda.i);
      AlphaMax=lambda.r;
      OmegaMax=lambda.i;
      return 1;
    }
    xl=lambda.r;
    yl=lambda.i;
   
    /* compute the Jacobian */
    if(fabs(xl)>eps)
      r=eps*fabs(xl);
    else
      r=eps*eps;
    xlp=xl+r;
    lambdap=rtoc(xlp,yl);
    make_z(z,delay,n,m,coef,lambdap);
    detp=cdeterm(z,n);
   jac[0]=(detp.r-det.r)/r;
   jac[2]=(detp.i-det.i)/r;
    if(fabs(yl)>eps)
      r=eps*fabs(yl);
    else
      r=eps*eps;
    ylp=yl+r;
    lambdap=rtoc(xl,ylp);
    make_z(z,delay,n,m,coef,lambdap);
    detp=cdeterm(z,n);
    jac[1]=(detp.r-det.r)/r;
    jac[3]=(detp.i-det.i)/r;
    r=jac[0]*jac[3]-jac[1]*jac[2];
    if(r==0){
      plintf(" singular jacobian \n");
      return -1;
    }
   xlp=(jac[3]*det.r-jac[1]*det.i)/r;
   ylp=(-jac[2]*det.r+jac[0]*det.i)/r;
    xl=xl-xlp;
    yl=yl-ylp;
    r=fabs(xlp)+fabs(ylp);
    lambda.r=xl;
    lambda.i=yl;
    if(r<err)
    { /* within the tolerance */
      process_root(lambda.r,lambda.i);
      AlphaMax=lambda.r;
      OmegaMax=lambda.i;
      rr[0]=AlphaMax;
      rr[1]=OmegaMax;
      return 1;
    }
    if(r>big){
      plintf("Failed to converge \n");
      return -1;
    }
  }
      
  plintf("Max iterates exceeded \n");
  return -1;
}
void process_root(real,im)
     double real,im;
{
  plintf("Root: %g + I %g \n",real,im); 
}
double get_arg(delay,coef,m,n,lambda)  
     COMPLEX lambda;
     double *coef;
     double *delay;
     int m,n;
{
  int i,j,k,km;
  COMPLEX *z;
  COMPLEX temp,eld;
  double arg;
  if(m==0)return(0);  /* no delays so don't use this! */
  z=(COMPLEX *)malloc(sizeof(COMPLEX)*n*n); 
  for(j=0;j<n;j++)
    for(i=0;i<n;i++){
      if(i==j)temp=lambda;
      else temp=rtoc(0.0,0.0);
      /* cprintn(temp); */
      z[i+j*n]=cdif(temp,rtoc(coef[i+j*n],0.0)); /* initialize the array */
    }
  for(k=0;k<m;k++){
    km=(k+1)*n*n;
    temp=rtoc(-delay[k],0.0); /* convert delay to complex number */
    eld=cexp2(cmlt(temp,lambda)); /* compute exp(-lambda*tau) */
    /* cprintn(eld); */
    for(j=0;j<n;j++)
      for(i=0;i<n;i++)
	z[i+j*n]=cdif(z[i+j*n],cmlt(eld,rtoc(coef[km+i+n*j],0.0)));
  }
  /*  the array is done  */
 /* cprintarr(z,n,n); */ 
  temp=cdeterm(z,n);
  /* cprint(lambda); 
  cprint(temp); 
  plintf(" \n"); */
   free(z); 
  arg=atan2(temp.i,temp.r);
  /*   plintf("%g %g %g \n",lambda.r,lambda.i,arg); */ 
  return(arg);
}   

int test_sign(old,new)
     double old,new;
{
  if(old>0.0&&new<0.0){
    if(old>2.9&&new<-2.9)return 1;
    return(0); /* doesnt pass threshold */
  }
  if(old<0.0&&new>0.0){
    if(old<-2.9&&new>2.9)return -1;
    return 0;
  }
  return 0;
}

/* code for establishing delay stability
   sign=plot_args(coef,delay,n,m,npts,amax,wmax)
    coef is a real array of length  (m+1)*n^2
    each n^2 block is the jacobian with respect to the mth delay
    m total delays
    n is size of system
    npts is number of pts on each part of contour
    contour is
      i wmax -----<---------    amax+i wmax
       |                            |
       V                            ^
       |                            |
     -i wmax ----->-----------  amax-i wmax
  
     sign is the number of roots in the contour using the argument
     principle
*/  

int plot_args(coef,delay,n,m,npts,almax,wmax)
     double *coef;
     int n,m,npts;
     double almax,wmax,*delay;
{
  int i;
  int sign=0;
  COMPLEX lambda;
  double x,y,arg,oldarg=0.0;
  double ds;  /* steplength */
  /* first the contour from i wmax -- -i wmax */
  ds=2*wmax/npts;
   x=0.0;
  for(i=0;i<npts;i++){
    y=wmax-i*ds;
    lambda=rtoc(x,y);
    arg=get_arg(delay,coef,m,n,lambda);
   /* plintf(" %d %g \n",i,arg); */
    sign=sign+test_sign(oldarg,arg);
    oldarg=arg;
 
  }
 /* lower contour   */
  y=-wmax;
  ds=almax/npts;
  for(i=0;i<npts;i++){
    x=i*ds;
    lambda=rtoc(x,y);
    arg=get_arg(delay,coef,m,n,lambda);
/*        plintf(" %d %g \n",i+npts,arg); */
       sign=sign+test_sign(oldarg,arg);
    oldarg=arg;
 
  }
/* right contour */
 x=almax;
 ds=2*wmax/npts;
  for(i=0;i<npts;i++){
    y=-wmax+i*ds;
    lambda=rtoc(x,y);
    arg=get_arg(delay,coef,m,n,lambda);
/*     plintf(" %d %g \n",i+2*npts,arg); */
      sign=sign+test_sign(oldarg,arg);
    oldarg=arg;
 
  }
 
/* top contour */
  y=wmax;
  ds=almax/npts;
  for(i=0;i<npts;i++){
    x=almax-i*ds;
    lambda=rtoc(x,y);
    arg=get_arg(delay,coef,m,n,lambda);
/*         plintf(" %d %g \n",i+3*npts,arg); */
    sign=sign+test_sign(oldarg,arg);
    oldarg=arg;
  
  }
  return sign;
}
 















