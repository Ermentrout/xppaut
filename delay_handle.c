#include "delay_handle.h"
#include "dde_po.h"
#include "parserslow.h"
#include "ggets.h"
#include "integrate.h"
void free(void *ptr);
void *calloc();
/* #include "stdfix.h"   */
/*   This handles the delay stuff    */

#include <stdio.h>
#include <math.h>
#include "xpplim.h"
#include "getvar.h"

double AlphaMax=2,OmegaMax=2;

double *DelayWork;
int LatestDelay;
int MaxDelay;
int DelayFlag=0;

int NDelay,del_stab_flag,WhichDelay,DelayGrid=1000;

/* History bookkeeping.  Samples live on a uniform grid of spacing DelayDD
   (=|DELTA_T| at allocation time) in a ring buffer.  Sample k steps back
   from the newest one is stored at ring index (LatestDelay+k)%MaxDelay.
   DelayTLatest is the *time stamp* of the newest sample.  Time stamps make
   two things possible:
     (a) get_delay() can be anchored at the actual stage time of the
         integrator (RK stages, DP stages, ...) rather than at the last
         stored step, which restores the order of the integrator, and
     (b) storing a sample whose time equals the newest time stamp
         overwrites it instead of pushing a duplicate. */
double DelayTLatest=0.0,DelayDD=1.0;
int DelayHaveTime=0;
double variable_shift[2][MAXODE];
double delay_list[MAXDELAY];

double evaluate();
extern double DELTA_T,T0,DELAY;
extern int NODE,NCON,NSYM,NSYM_START,NCON_START,NMarkov;
 
extern char delay_string[MAXODE][80];
extern double variables[];
extern int NVAR;



/* del_stab_flag modes (see parserslow2.c do_delay / do_delay_shift):
     > 0   normal: read the history buffer
       0   scanning: collect the distinct delays, return current variable
      -1   Jacobian pieces for the characteristic matrix (variable_shift)
      -2   steady state: every delayed variable is just the current variable
           (this is what AUTO and any equilibrium solver must see)
      -3   periodic orbit solver (dde_po.c): the delayed value is read off
           the periodic orbit being computed                              */

#define DELAY_SAME(a,b) (fabs((a)-(b))<=1e-12*(1.0+fabs(a)))

double delay_stab_eval(double delay,int var)
/* this returns appropriate values for delay jacobian */
{
  int i;

  if(del_stab_flag==-3){      /* periodic orbit of a DDE (dde_po.c): x(t-tau) read off the orbit */
    if(var<1||var>NODE)return(GETVAR(var));
    return(ddepo_delay(delay,var));
  }

  if(del_stab_flag==-2)       /* steady state:  x(t-tau)=x(t) */
    return(GETVAR(var));

  if(del_stab_flag==0) /* search for all delays  */
    {
      for(i=0;i<NDelay;i++){
        if(DELAY_SAME(delay,delay_list[i]))
          return(GETVAR(var));
      }
      if(NDelay<MAXDELAY){
        delay_list[NDelay]=delay;
        NDelay++;
      }
      return(GETVAR(var));
    }
 /*  now we must determine the value to return  */
 /*  del_stab_flag =-1 */
  for(i=0;i<NDelay;i++){
    if(DELAY_SAME(delay,delay_list[i]))
      if(i==WhichDelay)
        return variable_shift[1][var-1];
  }
  return variable_shift[0][var-1];
}

/* Run f() with all delayed variables replaced by the current variable.
   Used by AUTO so that F(x,x,...,x;p)=0 is what is continued. */
int delay_enter_steady(void)
{
  int old=del_stab_flag;
  if(DelayFlag)del_stab_flag=-2;
  return old;
}
void delay_leave_steady(int old)
{
  del_stab_flag=old;
}

int alloc_delay(double big)
{
 int n,i;
 double dd=fabs(DELTA_T);

 if(dd==0.0)dd=1.0;
 /* enough samples to cover [0,big] plus guard points for the cubic stencil
    (the old code could be one or two samples short, and the stencil then
    wrapped around the ring buffer and picked up the newest data) */
 n=(int)ceil(big/dd-1e-9)+4;
 if(n<8)n=8;

 MaxDelay=n;
 LatestDelay=1;
 DelayFlag=0;
 DelayHaveTime=0;
 DelayDD=dd;
 DelayWork=(double *)calloc((size_t)n*(size_t)NODE,sizeof(double));
 if(DelayWork==NULL){
  err_msg("Could not allocate memory for Delay");
  return(0);
 }
 DelayFlag=1;
 NDelay=0;
 WhichDelay=-1;
 del_stab_flag=1;
 for(i=0;i<n*(NODE );i++)DelayWork[i]=0.0;
 return(1);
}

void free_delay()
{
 if(DelayFlag)free(DelayWork);
 DelayFlag=0;
}

/* push (or, if t equals the newest time stamp, overwrite) a sample */
void stor_delay(double *y,double t)
{
 int i,in;
 int nodes=NODE;
 if(DelayFlag==0)return;
 if(DelayHaveTime&&fabs(t-DelayTLatest)<=1e-6*DelayDD){
   in=LatestDelay*nodes;
   for(i=0;i<nodes;i++)DelayWork[i+in]=y[i];
   return;
 }
 --LatestDelay;
 if(LatestDelay<0)LatestDelay+=MaxDelay;
 in=LatestDelay*nodes;
 for(i=0;i<nodes;i++)DelayWork[i+in]=y[i];
 DelayTLatest=t;
 DelayHaveTime=1;
}

double get_delay_old(int in,double tau)
{
 double x=tau/DelayDD;
 int n1=(int)x;
 int n2=n1+1;
 int nodes=NODE;
 int i1,i2;
 double x1,x2;
 if(tau<0.0||tau>DELAY){
   err_msg("Delay negative or too large");
   stop_integration();
   return(0.0);
 }
 i1=(n1+LatestDelay)%MaxDelay;
 i2=(n2+LatestDelay)%MaxDelay;
 x1=DelayWork[in+(nodes )*i1];
 x2=DelayWork[in+(nodes )*i2];
 return(x1+(x-n1)*(x2-x1));
}

void polint(xa,ya,n,x,y,dy)
     double *ya,*xa,*y,*dy,x;
     int n;
{
  int i,m,ns=1;
  double den,dif,dift,h0,hp,w;
  double c[10],d[10];
  dif=fabs(x-xa[0]);
  for(i=1;i<=n;i++){
    if( (dift=fabs(x-xa[i-1]))<dif){
      ns=i;
      dif=dift;
    }
    c[i-1]=ya[i-1];
    d[i-1]=ya[i-1];
  }
  *y=ya[(ns--) -1];
  for(m=1;m<n;m++){
    for(i=1;i<=n-m;i++){
      h0=xa[i-1]-x;
      hp=xa[i+m-1]-x;
      w=c[i]-d[i-1];
      if((den=h0-hp)==0.0)return ;
      den=w/den;
      d[i-1]=hp*den;
      c[i-1]=h0*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns]:d[ns-- -1]));
  }
}

/* Delayed value x_in(t_now - tau), t_now being the time at which the
   right-hand side is currently being evaluated (variables[0]).

   Samples are on a grid of spacing DelayDD ending at DelayTLatest.  The
   integrator is always somewhere in [DelayTLatest, DelayTLatest+DelayDD]
   while it works on the next step, so t_now-tau is located at

         s = (tau - (t_now-DelayTLatest)) / DelayDD     (grid units back)

   and we interpolate with the cubic through the four nearest samples.
   When tau is shorter than the offset into the step (s<0) the value is
   extrapolated from the four newest samples (at most one grid step).

   If t_now is not in the current step -- e.g. rhs() called with t=0 from a
   routine that does not care about time -- we fall back to the old
   behaviour of anchoring at the newest sample.                          */
double get_delay(int in,double tau)
{
 double dd=DelayDD,s,off,u,y;
 double sgn=(DELTA_T<0.0)?-1.0:1.0;
 int nodes=NODE;
 int m,k0,i0,i1,i2,i3;

 if(tau<0.0||tau>DELAY*(1.0+1e-12)+1e-12){
   err_msg("Delay negative or too large");
   stop_integration();
   return(0.0);
 }
 if(tau==0.0) /* check for zero delay and ignore the rest */
   return DelayWork[in+nodes*(LatestDelay%MaxDelay)];

 off=0.0;
 if(DelayHaveTime){
   off=(variables[0]-DelayTLatest)*sgn;
   if(off<-1e-6*dd||off>dd*(1.0+1e-6))off=0.0;
   else if(off<0.0)off=0.0;
 }
 s=(tau-off)/dd;
 m=(int)floor(s);
 k0=m-1;
 if(k0<0)k0=0;
 if(k0>MaxDelay-4)k0=MaxDelay-4;
 u=s-k0;       /* position within stencil nodes 0,1,2,3 */

 i0=(k0+LatestDelay)%MaxDelay;
 i1=(k0+1+LatestDelay)%MaxDelay;
 i2=(k0+2+LatestDelay)%MaxDelay;
 i3=(k0+3+LatestDelay)%MaxDelay;
 y= -(u-1.0)*(u-2.0)*(u-3.0)/6.0 * DelayWork[in+nodes*i0]
    + u*(u-2.0)*(u-3.0)/2.0      * DelayWork[in+nodes*i1]
    - u*(u-1.0)*(u-3.0)/2.0      * DelayWork[in+nodes*i2]
    + u*(u-1.0)*(u-2.0)/6.0      * DelayWork[in+nodes*i3];
 return(y);
}

/*  Handling of the initial data  */
int do_init_delay(double big)
{
 double t=T0,old_t,y[MAXODE];
 int i,nt,j;
 int len;
 int *del_form[MAXODE];

 /* if dt changed since allocation the grid is stale: rebuild it */
 if(DelayFlag&&fabs(fabs(DELTA_T)-DelayDD)>1e-12*DelayDD){
   free_delay();
   if(alloc_delay(big)==0)return(0);
 }
 nt=MaxDelay-1;
 NCON=NCON_START;
 NSYM=NSYM_START;
 for(i=0;i<(NODE );i++){
   del_form[i]=(int *)calloc(200,sizeof(int));
   if(del_form[i]==NULL){
     err_msg("Failed to allocate delay formula ...");
     for(j=0;j<i;j++)free(del_form[j]);
     NCON=NCON_START;
     NSYM=NSYM_START;
     return(0);
   }
   if(add_expr(delay_string[i],del_form[i],&len)){
     err_msg("Illegal delay expression");
     for(j=0;j<=i;j++)free(del_form[j]);
     NCON=NCON_START;
     NSYM=NSYM_START;
     return(0);
   }
 }        /*  Okay all formulas are cool... */
 LatestDelay=1;
 DelayHaveTime=0;

 get_val("t",&old_t);

 for(i=nt;i>=0;i--){
   t=T0-DelayDD*i;
   set_val("t",t);
   for(j=0;j<(NODE );j++)
     y[j]=evaluate(del_form[j]);
   stor_delay(y,t);
 }
 for(j=0;j<(NODE );j++)free(del_form[j]);
 NCON=NCON_START;
 NSYM=NSYM_START;
 set_val("t",old_t);
 return(1);
}
