/*  Glue between the DDE periodic orbit solver (dde_po.c) and XPP's AUTO window.

    do_auto() calls dde_po_run() for  Run -> Start -> Periodic  when the model
    has delays.  The starting cycle is the last simulation, which should be
    exactly one period of a stable oscillation (as for AUTO); the period is
    the length of the stored time interval.  The continuation parameter,
    range, step sizes, NTST, NCOL, NMX and NPR come from the AUTO windows.
    Every accepted point goes to the bifurcation diagram like an AUTO point.

    The mesh is refined automatically (NTST doubles) when the defect indicator says the
    orbit is not resolved, so NTST is only the starting value.

    A grabbed Hopf (HB) point of a steady-state branch also works: Grab it, then
    Run -> Periodic.  The first-order orbit from the critical eigenvector of the
    characteristic matrix is corrected with a small fixed amplitude and the branch is
    followed with growing amplitude (dde_po_hopf_run).

    Every plotted point is labeled and its orbit is kept in memory, so a point of a
    delay periodic branch can be grabbed and extended (Run -> Extend): the sign of Ds
    gives the direction in the parameter.  The stored orbits live only in this session.

    Not done yet: stability (Floquet multipliers), bifurcation detection,
    two-parameter continuation.  All points are drawn
    as open circles ("stability unknown").
*/
#include "dde_po.h"
#include "auto_f2c.h"
#include "auto_nox.h"
#include "diagram.h"
#include "derived.h"
#include "tabular.h"
#include "ggets.h"
#include "xpplim.h"
#include "autlim.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SPER 3
#define UPER 4
#define DD_MAXN 64

extern BIFUR Auto;
extern char upar_names[MAXPAR][11];
extern int AutoPar[8];
extern GRABPT grabpt;
int auto_delay_hopf_matrices(integer ndm,double *u,double omega,double *Cr,double *S);
extern int NODE,NAutoPar,storind,DiagFlag,AutoTwoParam,Auto_index_to_array[8];
extern int del_stab_flag;
extern double constants[];
extern float **storage;
extern int (*rhs)();
extern DIAGRAM *bifd;
void refreshdisplay(void);
int redraw_auto_menus(void);

/* ---------------- model through XPP ---------------- */
static void xp_rhs(void *ctx,double t,const double *x,double *f)
{
  double y[DD_MAXN];
  int i,old;
  (void)ctx;
  for(i=0;i<NODE;i++)y[i]=x[i];
  old=del_stab_flag;
  del_stab_flag=-3;                 /* delay(x,tau) reads the orbit being solved */
  rhs(t,y,f,NODE);
  del_stab_flag=old;
}

static void xp_setpar(void *ctx,int ip,double v)
{
  (void)ctx;
  constants[Auto_index_to_array[ip]]=v;
  evaluate_derived();
  redo_all_fun_tables();
}

static double xp_getpar(void *ctx,int ip)
{
  (void)ctx;
  return constants[Auto_index_to_array[ip]];
}

/* ---------------- starting cycle: the last simulation ---------------- */
typedef struct { double t0,T; } SimOrbit;

static void sim_orbit(void *ctx,double s,double *x)
{
  SimOrbit *o=(SimOrbit *)ctx;
  double t=o->t0+s*o->T,lam;
  int lo=0,hi=storind-1,j;
  while(hi-lo>1){
    int mid=(lo+hi)/2;
    if(storage[0][mid]<=t)lo=mid;else hi=mid;
  }
  lam=(t-storage[0][lo])/(storage[0][hi]-storage[0][lo]);
  if(lam<0.0)lam=0.0;
  if(lam>1.0)lam=1.0;
  for(j=0;j<NODE;j++)
    x[j]=(1.0-lam)*storage[j+1][lo]+lam*storage[j+1][hi];
}

/* ---------------- sending points to the diagram ---------------- */
typedef struct {
  DDEPO *P;
  int ibr,ntot,npr,nplot,lastpt;
  double evr[NAUTO],evi[NAUTO];
  int flags;
  double par0,lastpar,par20,lastpar2;
  int lab0,startitp;
  /* torus curves: the rotation angle and the strong resonances it passes */
  int have_th,nres;
  double th_prev,p1_prev,p2_prev,th0;
  struct { int kind; double p1,p2; } res[16];
  /* stability (Floquet multipliers) */
  int nstable,nunstable,nfloqfail,have_prev,prev_unst,nev;
  double prev_par,prev_dr,prev_di,prev_mod,worst_triv;
  struct { double plo,phi,pest; int kind; int to_unstable; } ev[16];
  int nsp;
  struct { int kind,lab; double par,T; } sp[16];
} RunInfo;

static int next_branch(void)
{
  DIAGRAM *d=bifd;
  int mx=0;
  if(DiagFlag==0)return 1;
  while(d!=NULL){
    if(abs(d->ibr)>mx)mx=abs(d->ibr);
    d=d->next;
  }
  return mx+1;
}

/* orbits of the plotted points, for restarting (Run -> Extend) */
typedef struct {
  int ibr,ntot,lab,n,ntst,ncol;
  double T,par[8];
  double *X;
} SavedPt;
#define MAXSAVED 4000
static SavedPt *saved=NULL;
static int nsaved=0;

static void save_point(int ibr,int ntot,int lab,const DDEPO *P,double T,const double *aux)
{
  SavedPt *sp;
  int i;
  if(saved==NULL){
    saved=(SavedPt *)calloc(MAXSAVED,sizeof(SavedPt));
    if(saved==NULL)return;
  }
  if(nsaved>=MAXSAVED){             /* drop the oldest */
    free(saved[0].X);
    memmove(saved,saved+1,sizeof(SavedPt)*(MAXSAVED-1));
    nsaved=MAXSAVED-1;
  }
  sp=&saved[nsaved];
  sp->X=(double *)malloc(sizeof(double)*(size_t)P->nmesh*P->n);
  if(sp->X==NULL)return;
  memcpy(sp->X,P->X,sizeof(double)*(size_t)P->nmesh*P->n);
  sp->ibr=ibr;sp->ntot=ntot;sp->lab=lab;sp->n=P->n;sp->ntst=P->ntst;sp->ncol=P->ncol;sp->T=T;
  for(i=0;i<8;i++)sp->par[i]=aux[i];
  nsaved++;
}

/* Forget every stored orbit.  Called when the AUTO diagram is destroyed
   (File -> Reset diagram): the orbits are keyed by branch and label, and both
   numberings start again from 1, so a stale entry could otherwise be taken for
   an orbit of the new diagram, and the memory would only ever grow. */
void dde_po_clear_saved(void)
{
  int i;
  if(saved==NULL)return;
  for(i=0;i<nsaved;i++){free(saved[i].X);saved[i].X=NULL;}
  free(saved);
  saved=NULL;
  nsaved=0;
}

static int max_label(void)
{
  DIAGRAM *d=bifd;
  int mx=0;
  if(DiagFlag==0)return 0;
  while(d!=NULL){
    if(d->lab>mx)mx=d->lab;
    d=d->next;
  }
  return mx;
}

/* kind of a stability change from the multiplier that is outside the unit circle:
   1 = real, negative (period doubling), 2 = real, positive (fold / branch point), 3 = complex pair (torus) */
static int mult_kind(double mr,double mi)
{
  if(fabs(mi)<=1e-6*(1.0+fabs(mr)))return mr<0.0?1:2;
  return 3;
}

static int xp_report(void *rctx,int istep,double par,double T,const double *X,int status)
{
  RunInfo *R=(RunInfo *)rctx;
  DDEPO *P=R->P;
  int n=P->n,c,g,i,lab,itp,labeled,stable=0,special=(status&16)!=0;
  double uhi[NAUTO],ulo[NAUTO],ubar[NAUTO],u0[NAUTO],aux[20],norm=0.0;
  if(istep==0)R->par0=par;
  R->lastpar=par;
  R->flags|=status;
  R->lastpt=istep;
  /* As in AUTO: every point is plotted, but only the first and the last point (EP) and every
     NPR-th point are labeled; only labeled points keep their orbit (for Import and Extend). */
  for(c=0;c<n;c++){
    double lo=X[c],hi=X[c],sum=0.0;
    for(g=0;g<P->nmesh;g++){
      double v=X[g*n+c];
      if(v<lo)lo=v;
      if(v>hi)hi=v;
      sum+=v;
      norm+=v*v;
    }
    uhi[c]=hi;ulo[c]=lo;ubar[c]=sum/P->nmesh;u0[c]=X[c];
  }
  norm=sqrt(norm/P->nmesh);
  for(i=0;i<20;i++)aux[i]=0.0;
  for(i=0;i<NAutoPar&&i<8;i++)aux[i]=constants[Auto_index_to_array[i]];
  aux[10]=T;
  if(special&&P->event==2){
    /* a torus needs a complex pair on the unit circle; a real pair with product 1 is a neutral saddle */
    int k,cplx=0;
    for(k=0;k<P->nmu;k++)
      if(fabs(P->mu_i[k])>1e-8&&fabs(hypot(P->mu_r[k],P->mu_i[k])-1.0)<1e-3)cplx=1;
    if(!cplx)return 0;
  }
  R->ntot++;
  labeled=(special||istep==0||(status&2)||(R->npr>0&&R->ntot%R->npr==0));
  if(labeled)R->nplot++;
  lab=labeled?R->lab0+R->nplot:0;
  /* 9 = end point (EP), 7 = period doubling (PD), 8 = torus (TR), 5 = fold of cycles (LP), 4 = labeled point, 0 = plain */
  itp=special?(P->event==1?7:(P->event==2?8:5)):((istep==0||(status&2))?9:(labeled?4:0));
  /* Floquet multipliers: stability, and the multiplier window */
  {
    double mr[70],mi[70],te=0.0;
    int nm,k,unst=0;
    if(P->locate&&P->mu_ok){          /* computed by the solver for its test functions */
      nm=P->nmu;te=P->mu_triv;
      for(k=0;k<nm&&k<70;k++){mr[k]=P->mu_r[k];mi[k]=P->mu_i[k];}
    }
    else nm=ddepo_floquet(P,mr,mi,(n<6?6:n),&te);
    double dr=0.0,di=0.0,dmod=0.0;
    for(k=0;k<n&&k<NAUTO;k++){R->evr[k]=0.0;R->evi[k]=0.0;}
    stable=0;
    if(nm>0){
      for(k=0;k<nm;k++){
	double md=sqrt(mr[k]*mr[k]+mi[k]*mi[k]);
	if(md>1.0)unst++;
	if(k==0){dr=mr[0];di=mi[0];dmod=md;}
      }
      for(k=0;k<n&&k<nm;k++){R->evr[k]=mr[k];R->evi[k]=mi[k];}
      stable=(unst==0);
      if(te>R->worst_triv)R->worst_triv=te;
      if(!special){if(stable)R->nstable++;else R->nunstable++;}
      if(!special&&R->have_prev&&unst!=R->prev_unst&&R->nev<16&&fabs(R->prev_mod-1.0)>1e-6){   /* not if the previous point is itself on the boundary */
	/* a change: classify by the dominant multiplier on the unstable side */
	int to_unst=(unst>R->prev_unst);
	double cr=to_unst?dr:R->prev_dr,ci=to_unst?di:R->prev_di;
	double m0=R->prev_mod,m1=dmod,est=par;
	if(m1!=m0&&(m0-1.0)*(m1-1.0)<0.0)est=R->prev_par+(par-R->prev_par)*(1.0-m0)/(m1-m0);
	else est=0.5*(R->prev_par+par);
	R->ev[R->nev].plo=R->prev_par;R->ev[R->nev].phi=par;R->ev[R->nev].pest=est;
	R->ev[R->nev].kind=mult_kind(cr,ci);R->ev[R->nev].to_unstable=to_unst;
	R->nev++;
      }
      if(!special){R->have_prev=1;R->prev_unst=unst;R->prev_par=par;R->prev_dr=dr;R->prev_di=di;R->prev_mod=dmod;}
    }
    else if(!special)R->nfloqfail++;
  }
  /* sign convention of AUTO's diagram: ibr<0 periodic; ntot<0 stable (filled circle), ntot>0 unstable (open) */
  add_point(aux,T,uhi,ulo,ubar,norm,stable?SPER:UPER,istep==0?0:1,lab,1,Auto.icp1,Auto.icp2,Auto.icp3,Auto.icp4,
	    AutoTwoParam,R->evr,R->evi);
  if(DiagFlag==0){
    edit_start(-R->ibr,stable?-R->ntot:R->ntot,itp,lab,1,norm,uhi,ulo,u0,ubar,aux,T,n,Auto.icp1,Auto.icp2,
	       Auto.icp3,Auto.icp4,R->evr,R->evi);
    DiagFlag=1;
  }
  else
    add_diagram(-R->ibr,stable?-R->ntot:R->ntot,itp,lab,1,norm,uhi,ulo,u0,ubar,aux,T,n,Auto.icp1,Auto.icp2,
		Auto.icp3,Auto.icp4,AutoTwoParam,R->evr,R->evi);
  if(labeled)save_point(-R->ibr,R->ntot,lab,P,T,aux);
  if(special&&R->nsp<16){
    R->sp[R->nsp].kind=P->event;R->sp[R->nsp].par=par;R->sp[R->nsp].lab=lab;R->sp[R->nsp].T=T;
    R->nsp++;
  }
  if(status&(2|4))refreshdisplay();
  return 0;
}

/* Esc or the ABORT button in the AUTO window stops a run (the same check AUTO's own engine makes; it
   also lets the window redraw while a run is going) */
extern int byeauto_(int *iflag);
static int xp_poll(void *ctx)
{
  int flag=0;
  (void)ctx;
  byeauto_(&flag);
  return flag;
}

/* One summary line for the terminal: what was computed, where and why it ended */
static void run_summary(const char *what,int np,double T,int ntst0,const DDEPO *P,const RunInfo *R)
{
  static const char *why[]={
    "reached the end of the parameter range",
    "reached Nmax steps (raise Nmax and Extend from the last point to go on)",
    "amplitude collapsed: a Hopf point is near",
    "the step size fell below Dsmin without Newton converging: the orbit may be too sharp, try a larger NTST or NCOL, or a smaller Dsmax",
    "stopped",
    "stopped by the user (Esc or ABORT); the last point is labeled, Extend continues from it"};
  int st=P->stop;
  if(st<0||st>5)st=4;
  printf("%s: %d points, last period %g, mesh NTST=%d (started with %d, refined %d times)%s\n",
	 what,np,T,P->ntst,ntst0,P->nrefine,(R->flags&1)?", passed a turning point":"");
  printf("  ended at %s=%g: %s\n",upar_names[AutoPar[Auto.icp1]],R->lastpar,why[st]);
  if(P->nunres>0)
    printf("  WARNING: at %d points the mesh limit (3000 unknowns) was reached before the orbit was resolved\n"
	   "  (defect up to %.1e, tolerance %.1e): those points are less accurate\n",P->nunres,P->maxdefect,P->dtol);
  {
    int k;
    static const char *kind[]={"","a real multiplier crosses -1: period doubling",
			       "a real multiplier crosses +1: fold or branch point of cycles",
			       "a complex pair crosses the unit circle: torus (Neimark-Sacker)"};
    printf("  Stability (Floquet multipliers): stable at %d points, unstable at %d points",R->nstable,R->nunstable);
    if(R->nfloqfail>0)printf(", not computed at %d points",R->nfloqfail);
    printf("\n");
    for(k=0;k<R->nsp;k++)
      printf("  %s at %s=%.8g (label %d, period %g)\n",
	     R->sp[k].kind==1?"PERIOD DOUBLING (PD)":(R->sp[k].kind==2?"TORUS (Neimark-Sacker, TR)":"FOLD OF CYCLES (LP)"),
	     upar_names[AutoPar[Auto.icp1]],R->sp[k].par,R->sp[k].lab,R->sp[k].T);
    for(k=0;k<R->nev;k++){
      int j,covered=0;
      for(j=0;j<R->nsp;j++){   /* already located: do not print the rough estimate as well */
	double lo=R->ev[k].plo<R->ev[k].phi?R->ev[k].plo:R->ev[k].phi,hi=R->ev[k].plo<R->ev[k].phi?R->ev[k].phi:R->ev[k].plo;
	if(R->sp[j].par>=lo&&R->sp[j].par<=hi&&((R->sp[j].kind==1&&R->ev[k].kind==1)||(R->sp[j].kind==2&&R->ev[k].kind==3)))covered=1;
      }
      if(covered)continue;
      printf("  stability change between %s=%g and %g (near %g): the cycle becomes %s; %s\n",
	     upar_names[AutoPar[Auto.icp1]],R->ev[k].plo,R->ev[k].phi,R->ev[k].pest,
	     R->ev[k].to_unstable?"unstable":"stable",kind[R->ev[k].kind]);
    }
    if(R->worst_triv>1e-3)
      printf("  WARNING: the trivial multiplier is off by up to %.1e: the multipliers are not accurate, try more NTST\n",R->worst_triv);
  }
}

/* ---------------- the run ---------------- */
int dde_po_run(void)
{
  DDEPO P;
  SimOrbit so;
  RunInfo R;
  double p0,T,dsmin;
  int n=NODE,ntst,ncol,np,i;
  char msg[300];

  if(storind<10){
    err_msg("Run a simulation over exactly one period of a stable oscillation first");
    return 0;
  }
  if(n>DD_MAXN){err_msg("Too many variables for the delay periodic orbit solver");return 0;}
  if(Auto.icp1<0||Auto.icp1>=NAutoPar){err_msg("Select the continuation parameter in the Parameter window");return 0;}
  so.t0=storage[0][0];
  so.T=storage[0][storind-1]-storage[0][0];
  if(so.T<=0.0){err_msg("The last simulation has zero length");return 0;}
  p0=xp_getpar(NULL,Auto.icp1);
  if(p0<Auto.rl0||p0>Auto.rl1){
    sprintf(msg,"Parameter value %g is outside the range [%g,%g] set in the Numerics window",p0,Auto.rl0,Auto.rl1);
    err_msg(msg);
    return 0;
  }
  ntst=Auto.ntst;ncol=Auto.ncol;
  if(ntst<8)ntst=8;
  if(ncol<2)ncol=2;
  if(ncol>8)ncol=8;
  memset(&P,0,sizeof(P));
  P.rhs=xp_rhs;P.setpar=xp_setpar;P.getpar=xp_getpar;P.poll=xp_poll;P.pollctx=NULL;P.ctx=NULL;P.ip=Auto.icp1;
  P.locate=1;
  if(!ddepo_alloc(&P,n,ntst,ncol)){err_msg("Could not set up the periodic orbit solver");return 0;}
  if(n*ntst*ncol>3000){
    ddepo_free(&P);
    sprintf(msg,"Too large: %d variables x NTST %d x NCOL %d. Reduce NTST or NCOL (limit 3000 unknowns)",n,ntst,ncol);
    err_msg(msg);
    return 0;
  }
  ddepo_set_guess(&P,so.T,sim_orbit,&so);

  memset(&R,0,sizeof(R));
  R.P=&P;
  R.ibr=next_branch();
  R.lab0=max_label();
  R.npr=Auto.npr;
  redraw_auto_menus();
  /* AUTO's Dsmin (default 0.001) is far larger than the steps needed near a Hopf point, where
     the step must stay below the shrinking amplitude: use a much smaller floor */
  dsmin=fabs(Auto.dsmin)<1e-6?fabs(Auto.dsmin):1e-6;
  np=ddepo_continue(&P,Auto.rl0,Auto.rl1,(Auto.ds>=0.0)?1:-1,fabs(Auto.ds),dsmin,
		    fabs(Auto.dsmax),Auto.nmx,xp_report,&R);
  T=P.T;
  xp_setpar(NULL,Auto.icp1,p0);        /* leave the model at its starting parameter */
  ddepo_free(&P);
  for(i=0;i<1;i++)refreshdisplay();
  if(np==-2){
    err_msg("Newton did not converge on the starting orbit. Integrate over exactly one period of a STABLE oscillation (Initial conditions -> last, then run for one period), or raise NTST");
    return 0;
  }
  if(np<0){err_msg("Delay periodic orbit continuation failed to start");return 0;}
  run_summary("Delay periodic orbit continuation",np,T,ntst,&P,&R);
  ping();
  return 1;
}

/* ---------------- start at a grabbed Hopf point ---------------- */
int dde_po_hopf_run(void)
{
  static double Cr[DD_MAXN*DD_MAXN],S[DD_MAXN*DD_MAXN];
  DDEPO P;
  RunInfo R;
  double va[DD_MAXN],vb[DD_MAXN],u0[DD_MAXN],saved[8];
  double eps_list[5]={0.05,0.02,0.1,0.01,0.2};
  double per=grabpt.per,omega,res,pH,scale=1.0,dsmin,T;
  int n=NODE,ntst,ncol,np=-2,i,k,ibr;
  char msg[300];

  if(n>DD_MAXN){err_msg("Too many variables for the delay periodic orbit solver");return 0;}
  if(Auto.icp1<0||Auto.icp1>=NAutoPar){err_msg("Select the continuation parameter in the Parameter window");return 0;}
  if(!(per>0.0)||per>1e30){err_msg("The grabbed Hopf point has no finite period (frequency 0): cannot start a cycle here");return 0;}
  omega=2.0*M_PI/per;
  for(i=0;i<NAutoPar&&i<8;i++){
    saved[i]=constants[Auto_index_to_array[i]];
    constants[Auto_index_to_array[i]]=grabpt.par[i];
  }
  evaluate_derived();
  redo_all_fun_tables();
  for(i=0;i<n;i++){
    u0[i]=grabpt.u0[i];
    if(fabs(u0[i])>scale)scale=fabs(u0[i]);
  }
  pH=xp_getpar(NULL,Auto.icp1);
  if(pH<Auto.rl0||pH>Auto.rl1){
    sprintf(msg,"The Hopf point (parameter %g) is outside the range [%g,%g] set in the Numerics window",pH,Auto.rl0,Auto.rl1);
    err_msg(msg);
    goto restore_fail;
  }
  if(auto_delay_hopf_matrices((integer)n,u0,omega,Cr,S)){
    err_msg("Could not form the characteristic matrix at the Hopf point");
    goto restore_fail;
  }
  res=ddepo_hopf_vector(n,Cr,S,omega,va,vb);
  if(res<0.0||res>0.1){
    err_msg("This does not look like a Hopf point: the characteristic matrix is not singular at the stored frequency");
    goto restore_fail;
  }
  ntst=Auto.ntst;ncol=Auto.ncol;
  if(ntst<8)ntst=8;
  if(ncol<2)ncol=2;
  if(ncol>8)ncol=8;
  memset(&P,0,sizeof(P));
  P.rhs=xp_rhs;P.setpar=xp_setpar;P.getpar=xp_getpar;P.poll=xp_poll;P.pollctx=NULL;P.ctx=NULL;P.ip=Auto.icp1;
  P.locate=1;
  if(n*ntst*ncol>3000||!ddepo_alloc(&P,n,ntst,ncol)){
    err_msg("Too large: reduce NTST or NCOL (limit 3000 unknowns = variables x NTST x NCOL)");
    goto restore_fail;
  }
  ibr=next_branch();
  redraw_auto_menus();
  dsmin=fabs(Auto.dsmin)<1e-6?fabs(Auto.dsmin):1e-6;
  /* the first-order orbit is corrected with a small fixed amplitude; try a few sizes */
  for(k=0;k<5;k++){
    memset(&R,0,sizeof(R));
    R.P=&P;R.ibr=ibr;R.lab0=max_label();R.npr=Auto.npr;
    ddepo_hopf_guess(&P,u0,omega,va,vb,eps_list[k]*scale);
    np=ddepo_continue(&P,Auto.rl0,Auto.rl1,1,fabs(Auto.ds),dsmin,fabs(Auto.dsmax),Auto.nmx,xp_report,&R);
    if(np!=-2)break;
  }
  T=P.T;
  ddepo_free(&P);
  for(i=0;i<NAutoPar&&i<8;i++)constants[Auto_index_to_array[i]]=saved[i];
  evaluate_derived();
  redo_all_fun_tables();
  refreshdisplay();
  if(np==-5){err_msg("The first cycle lies outside [Par Min, Par Max]: widen the parameter range on the side where the cycles exist");return 0;}
  if(np==-2){err_msg("Newton did not converge on the cycle near the Hopf point (tried several small amplitudes); try more NTST/NCOL");return 0;}
  if(np<0){err_msg("Delay periodic orbit continuation from the Hopf point failed");return 0;}
  printf("Hopf start: omega=%g (period %g), cycles bifurcate towards %s parameter values (first cycle at %g)\n",
	 omega,per,(R.par0>pH)?"larger":"smaller",R.par0);
  run_summary("Delay periodic orbit continuation",np,T,ntst,&P,&R);
  ping();
  return 1;

restore_fail:
  for(i=0;i<NAutoPar&&i<8;i++)constants[Auto_index_to_array[i]]=saved[i];
  evaluate_derived();
  redo_all_fun_tables();
  return 0;
}

/* ---------------- restart from a point of a delay periodic branch ---------------- */
int dde_po_extend_run(void)
{
  DDEPO P;
  RunInfo R;
  SavedPt *sp=NULL;
  double keep[8],p0,T,dsmin;
  int n=NODE,i,np,ntst0;
  char msg[300];

  if(Auto.icp1<0||Auto.icp1>=NAutoPar){err_msg("Select the continuation parameter in the Parameter window");return 0;}
  for(i=nsaved-1;i>=0;i--){
    SavedPt *c=&saved[i];
    if(c->ibr==grabpt.ibr&&c->ntot==abs(grabpt.ntot)&&c->n==n&&
       fabs(c->T-grabpt.per)<=1e-6*fabs(c->T)&&
       fabs(c->par[Auto.icp1]-grabpt.par[Auto.icp1])<=1e-9*(1.0+fabs(c->par[Auto.icp1]))){
      sp=c;break;
    }
  }
  if(sp==NULL){
    err_msg("The orbit of this point is not stored: only points computed in this session by the delay periodic orbit solver can be extended");
    return 0;
  }
  for(i=0;i<NAutoPar&&i<8;i++){
    keep[i]=constants[Auto_index_to_array[i]];
    constants[Auto_index_to_array[i]]=sp->par[i];
  }
  evaluate_derived();
  redo_all_fun_tables();
  p0=xp_getpar(NULL,Auto.icp1);
  memset(&P,0,sizeof(P));
  P.rhs=xp_rhs;P.setpar=xp_setpar;P.getpar=xp_getpar;P.poll=xp_poll;P.pollctx=NULL;P.ctx=NULL;P.ip=Auto.icp1;
  ntst0=sp->ntst;
  P.locate=1;
  if(!ddepo_alloc(&P,n,sp->ntst,sp->ncol)){
    err_msg("Could not set up the periodic orbit solver");
    goto restore;
  }
  memcpy(P.X,sp->X,sizeof(double)*(size_t)P.nmesh*n);
  P.T=sp->T;
  memset(&R,0,sizeof(R));
  R.P=&P;
  R.ibr=next_branch();
  R.lab0=max_label();
  R.npr=Auto.npr;
  redraw_auto_menus();
  dsmin=fabs(Auto.dsmin)<1e-6?fabs(Auto.dsmin):1e-6;
  np=ddepo_continue(&P,Auto.rl0,Auto.rl1,(Auto.ds>=0.0)?1:-1,fabs(Auto.ds),dsmin,
		    fabs(Auto.dsmax),Auto.nmx,xp_report,&R);
  T=P.T;
  ddepo_free(&P);
  for(i=0;i<NAutoPar&&i<8;i++)constants[Auto_index_to_array[i]]=keep[i];
  evaluate_derived();
  redo_all_fun_tables();
  refreshdisplay();
  if(np==-5){
    sprintf(msg,"Parameter value %g is outside the range [%g,%g] set in the Numerics window",p0,Auto.rl0,Auto.rl1);
    err_msg(msg);
    return 0;
  }
  if(np<0){err_msg("Could not restart the delay periodic orbit continuation from this point");return 0;}
  {
    char what[160];
    sprintf(what,"Extend from point %d of branch %d (%s parameter values)",abs(grabpt.ntot),abs(grabpt.ibr),
	    (Auto.ds>=0.0)?"increasing":"decreasing");
    run_summary(what,np,T,ntst0,&P,&R);
  }
  ping();
  return 1;

restore:
  for(i=0;i<NAutoPar&&i<8;i++)constants[Auto_index_to_array[i]]=keep[i];
  evaluate_derived();
  redo_all_fun_tables();
  return 0;
}

/* ---------------- File -> Import orbit, and "draw labeled orbits" ---------------- */
/* The orbit of a labeled point of a delay periodic branch is not in AUTO's .s file (which
   is indexed by label and would give another orbit with the same label number): it is
   given here from the stored solution.  Returns 1 and malloc'd arrays ts (fraction of the
   period, nrow values from 0 to 1) and us (nrow x NODE values), or 0 if the point is not
   known (e.g. it is from a diagram loaded from a file). */
int dde_po_get_orbit(int ibr,int lab,double per,int maxrows,int *nrow,double **ts,double **us)
{
  int i,g,n,nm,st,nr;
  SavedPt *sp=NULL;
  for(i=nsaved-1;i>=0;i--)
    if(saved[i].ibr==ibr&&saved[i].lab==lab&&fabs(saved[i].T-per)<=1e-6*fabs(saved[i].T)){sp=&saved[i];break;}
  if(sp==NULL||sp->n!=NODE||maxrows<3)return 0;
  n=sp->n;nm=sp->ntst*sp->ncol;
  st=(nm+maxrows-2)/(maxrows-1);          /* thin the nodes if the data storage is small */
  if(st<1)st=1;
  nr=(nm+st-1)/st+1;
  *ts=(double *)malloc(sizeof(double)*nr);
  *us=(double *)malloc(sizeof(double)*(size_t)nr*n);
  if(*ts==NULL||*us==NULL){free(*ts);free(*us);return 0;}
  for(i=0;i<nr-1;i++){
    g=i*st;
    (*ts)[i]=(double)g/nm;
    memcpy(*us+(size_t)i*n,sp->X+(size_t)g*n,sizeof(double)*n);
  }
  (*ts)[nr-1]=1.0;                         /* closing point = first point */
  memcpy(*us+(size_t)(nr-1)*n,sp->X,sizeof(double)*n);
  *nrow=nr;
  printf("Delay periodic orbit of label %d: %d points over one period %g loaded\n",lab,nr,per);
  return 1;
}

/* ---------------- switch to the period-doubled branch at a grabbed PD point ---------------- */
int dde_po_pd_run(void)
{
  DDEPO P1,P2;
  RunInfo R;
  SavedPt *sp=NULL;
  double keep[8],p0,*v=NULL,mu=0.0,mm=0.0,scale=0.0,T,dsmin;
  double eps_list[5]={0.01,0.005,0.02,0.002,0.04};
  int n=NODE,i,k,np=-2,ibr,ntst2;
  char msg[300];

  if(Auto.icp1<0||Auto.icp1>=NAutoPar){err_msg("Select the continuation parameter in the Parameter window");return 0;}
  if(grabpt.itp%10!=7){err_msg("This is not a period doubling point");return 0;}
  for(i=nsaved-1;i>=0;i--){
    SavedPt *c=&saved[i];
    if(c->ibr==grabpt.ibr&&c->ntot==abs(grabpt.ntot)&&c->n==n&&
       fabs(c->T-grabpt.per)<=1e-6*fabs(c->T)&&
       fabs(c->par[Auto.icp1]-grabpt.par[Auto.icp1])<=1e-9*(1.0+fabs(c->par[Auto.icp1]))){
      sp=c;break;
    }
  }
  if(sp==NULL){
    err_msg("The orbit of this point is not stored: only points computed in this session by the delay periodic orbit solver can be used");
    return 0;
  }
  ntst2=2*sp->ntst;
  if(n*ntst2*sp->ncol>3000){
    err_msg("The doubled orbit is too large (limit 3000 unknowns = variables x NTST x NCOL): use a smaller NTST or NCOL for the first branch");
    return 0;
  }
  for(i=0;i<NAutoPar&&i<8;i++){
    keep[i]=constants[Auto_index_to_array[i]];
    constants[Auto_index_to_array[i]]=sp->par[i];
  }
  evaluate_derived();
  redo_all_fun_tables();
  p0=xp_getpar(NULL,Auto.icp1);
  memset(&P1,0,sizeof(P1));memset(&P2,0,sizeof(P2));
  P1.rhs=xp_rhs;P1.setpar=xp_setpar;P1.getpar=xp_getpar;P1.poll=xp_poll;P1.pollctx=NULL;P1.ctx=NULL;P1.ip=Auto.icp1;
  P2=P1;
  P2.locate=1;
  if(!ddepo_alloc(&P1,n,sp->ntst,sp->ncol)){err_msg("Could not set up the periodic orbit solver");goto restore;}
  memcpy(P1.X,sp->X,sizeof(double)*(size_t)P1.nmesh*n);
  P1.T=sp->T;
  v=(double *)malloc(sizeof(double)*(size_t)P1.nmesh*n);
  if(v==NULL||ddepo_floquet_mode(&P1,-1.0,v,&mu,&mm)||fabs(mu+1.0)>0.05||mm>1e-3){
    err_msg("Could not compute the critical Floquet eigenfunction at this point (is it a period doubling point?)");
    ddepo_free(&P1);free(v);goto restore;
  }
  for(i=0;i<P1.nmesh*n;i++)if(fabs(P1.X[i])>scale)scale=fabs(P1.X[i]);
  if(scale<=0.0)scale=1.0;
  if(!ddepo_alloc(&P2,n,ntst2,sp->ncol)){err_msg("Could not set up the periodic orbit solver");ddepo_free(&P1);free(v);goto restore;}
  ibr=next_branch();
  redraw_auto_menus();
  dsmin=fabs(Auto.dsmin)<1e-6?fabs(Auto.dsmin):1e-6;
  /* the first-order doubled orbit is corrected with a small fixed amplitude; try a few sizes */
  for(k=0;k<5;k++){
    memset(&R,0,sizeof(R));
    R.P=&P2;R.ibr=ibr;R.lab0=max_label();R.npr=Auto.npr;
    if(ddepo_pd_guess(&P2,&P1,v,eps_list[k]*scale)){np=-1;break;}
    np=ddepo_continue(&P2,Auto.rl0,Auto.rl1,1,fabs(Auto.ds),dsmin,fabs(Auto.dsmax),Auto.nmx,xp_report,&R);
    if(np!=-2)break;
  }
  T=P2.T;
  ddepo_free(&P1);free(v);
  for(i=0;i<NAutoPar&&i<8;i++)constants[Auto_index_to_array[i]]=keep[i];
  evaluate_derived();
  redo_all_fun_tables();
  refreshdisplay();
  if(np==-5){err_msg("The first doubled cycle lies outside [Par Min, Par Max]: widen the parameter range on the side where the doubled cycles exist");ddepo_free(&P2);return 0;}
  if(np==-2){err_msg("Newton did not converge on the doubled cycle (tried several small amplitudes); try more NTST/NCOL");ddepo_free(&P2);return 0;}
  if(np<0){err_msg("Continuation of the period-doubled branch failed");ddepo_free(&P2);return 0;}
  printf("Period-doubling switch at %s=%g: critical multiplier %.9f, doubled cycles exist for %s parameter values (first at %g)\n",
	 upar_names[AutoPar[Auto.icp1]],p0,mu,(R.par0>p0)?"larger":"smaller",R.par0);
  {
    char what[160];
    sprintf(what,"Period-doubled branch from label %d",grabpt.lab);
    run_summary(what,np,T,ntst2,&P2,&R);
  }
  ddepo_free(&P2);
  ping();
  (void)msg;
  return 1;

restore:
  for(i=0;i<NAutoPar&&i<8;i++)constants[Auto_index_to_array[i]]=keep[i];
  evaluate_derived();
  redo_all_fun_tables();
  return 0;
}

/* ---------------- two-parameter curve of a fold of cycles ---------------- */
#define XP_TWOPAR 4   /* Auto.plot value of the Two par axes (P_P in auto_nox.c) */
static int xp_report2(void *rctx,int istep,double p1,double p2,double T,const double *X,int status)
{
  RunInfo *R=(RunInfo *)rctx;
  DDEPO *P=R->P;
  int n=P->n,c,g,i,lab,itp,labeled;
  double uhi[NAUTO],ulo[NAUTO],ubar[NAUTO],u0[NAUTO],aux[20],norm=0.0;
  if(istep==0){R->par0=p1;R->par20=p2;}
  if(R->startitp==8){
    /* theta is not reduced modulo 2 pi.  R1/R2: it crosses a multiple of pi (multipliers at +1 / -1);
       R4/R3: the folded angle crosses pi/2 / 2pi/3 */
    double th=P->theta,fold_prev,fold_cur;
    if(R->have_th){
      double a=R->th_prev,b=th;
      int ka=(int)floor(a/M_PI),kb=(int)floor(b/M_PI);
      if(ka!=kb&&R->nres<16){
	int kk=(ka<kb)?kb:ka;
	double fr=(kk*M_PI-a)/(b-a);
	R->res[R->nres].kind=(kk%2==0)?1:2;
	R->res[R->nres].p1=R->p1_prev+fr*(p1-R->p1_prev);R->res[R->nres].p2=R->p2_prev+fr*(p2-R->p2_prev);
	R->nres++;
      }
      fold_prev=fmod(fabs(a),2.0*M_PI);if(fold_prev>M_PI)fold_prev=2.0*M_PI-fold_prev;
      fold_cur=fmod(fabs(b),2.0*M_PI);if(fold_cur>M_PI)fold_cur=2.0*M_PI-fold_cur;
      {
	static const double tgt[2]={2.0*M_PI/3.0,0.5*M_PI};
	int q;
	for(q=0;q<2;q++)
	  if((fold_prev-tgt[q])*(fold_cur-tgt[q])<0.0&&R->nres<16){
	    double fr=(tgt[q]-fold_prev)/(fold_cur-fold_prev);
	    R->res[R->nres].kind=3+q;                      /* 3: R3, 4: R4 */
	    R->res[R->nres].p1=R->p1_prev+fr*(p1-R->p1_prev);R->res[R->nres].p2=R->p2_prev+fr*(p2-R->p2_prev);
	    R->nres++;
	  }
      }
    }
    else R->th0=th;
    R->have_th=1;R->th_prev=th;R->p1_prev=p1;R->p2_prev=p2;
  }
  R->lastpar=p1;R->lastpar2=p2;
  R->flags|=status;
  R->lastpt=istep;
  for(c=0;c<n;c++){
    double lo=X[c],hi=X[c],sum=0.0;
    for(g=0;g<P->nmesh;g++){
      double v=X[g*n+c];
      if(v<lo)lo=v;
      if(v>hi)hi=v;
      sum+=v;
      norm+=v*v;
    }
    uhi[c]=hi;ulo[c]=lo;ubar[c]=sum/P->nmesh;u0[c]=X[c];
  }
  norm=sqrt(norm/P->nmesh);
  for(i=0;i<20;i++)aux[i]=0.0;
  for(i=0;i<NAutoPar&&i<8;i++)aux[i]=constants[Auto_index_to_array[i]];
  aux[Auto.icp1]=p1;aux[Auto.icp2]=p2;
  aux[10]=T;
  for(i=0;i<NAUTO;i++){R->evr[i]=0.0;R->evi[i]=0.0;}
  R->ntot++;
  labeled=(istep==0||(status&2)||(R->npr>0&&R->ntot%R->npr==0));
  if(labeled)R->nplot++;
  lab=labeled?R->lab0+R->nplot:0;
  itp=(istep==0)?R->startitp:((status&2)?9:(labeled?4:0));   /* the start is the LP or PD point, the end an EP */
  add_point(aux,T,uhi,ulo,ubar,norm,UPER,istep==0?0:1,lab,1,Auto.icp1,Auto.icp2,Auto.icp3,Auto.icp4,
	    AutoTwoParam,R->evr,R->evi);
  if(DiagFlag==0){
    edit_start(-R->ibr,R->ntot,itp,lab,1,norm,uhi,ulo,u0,ubar,aux,T,n,Auto.icp1,Auto.icp2,
	       Auto.icp3,Auto.icp4,R->evr,R->evi);
    DiagFlag=1;
  }
  else
    add_diagram(-R->ibr,R->ntot,itp,lab,1,norm,uhi,ulo,u0,ubar,aux,T,n,Auto.icp1,Auto.icp2,
		Auto.icp3,Auto.icp4,AutoTwoParam,R->evr,R->evi);
  if(status&2)refreshdisplay();
  return 0;
}

/* Run -> (two parameters) on a grabbed LP (kind 0), PD (kind 1) or TR (kind 2) point of a delay periodic branch */
static int curve2_run(int kind)
{
  DDEPO P;
  RunInfo R;
  SavedPt *sp=NULL;
  double keep[8],pmin,pmax,qmin,qmax,dsmin,p1s,p2s;
  int n=NODE,i,np,iplot;
  static const char *why[]={"left the box (the plot window of the Two par axes, or Par Min/Max)",
			    "reached Nmax steps (raise Nmax)","","Newton failed even with the smallest step (Dsmin)",
			    "stopped","stopped by the user (Esc or ABORT)"};

  if(Auto.icp1<0||Auto.icp1>=NAutoPar||Auto.icp2<0||Auto.icp2>=NAutoPar||Auto.icp1==Auto.icp2){
    err_msg("Choose two different parameters (Par1, Par2) in the Parameter window");return 0;}
  if(grabpt.itp%10!=(kind==0?5:(kind==1?7:8))){err_msg(kind==0?"This is not a fold (LP) point":(kind==1?"This is not a period doubling (PD) point":"This is not a torus (TR) point"));return 0;}
  for(i=nsaved-1;i>=0;i--){
    SavedPt *c=&saved[i];
    if(c->ibr==grabpt.ibr&&c->ntot==abs(grabpt.ntot)&&c->n==n&&
       fabs(c->T-grabpt.per)<=1e-6*fabs(c->T)&&
       fabs(c->par[Auto.icp1]-grabpt.par[Auto.icp1])<=1e-9*(1.0+fabs(c->par[Auto.icp1]))){
      sp=c;break;
    }
  }
  if(sp==NULL){
    err_msg("The orbit of this point is not stored: only points computed in this session by the delay periodic orbit solver can be used");
    return 0;
  }
  if(n*sp->ntst*sp->ncol>(kind==2?480:800)){
    err_msg(kind==2?"The torus curve needs about three times the unknowns of the orbit: use NTST x NCOL x variables <= 480":"The curve needs about twice the unknowns of the orbit: use NTST x NCOL x variables <= 800");
    return 0;
  }
  for(i=0;i<NAutoPar&&i<8;i++){
    keep[i]=constants[Auto_index_to_array[i]];
    constants[Auto_index_to_array[i]]=sp->par[i];
  }
  evaluate_derived();
  redo_all_fun_tables();
  p1s=sp->par[Auto.icp1];p2s=sp->par[Auto.icp2];
  /* the box in which the curve is followed */
  if(Auto.plot==XP_TWOPAR&&Auto.xmax>Auto.xmin&&Auto.ymax>Auto.ymin){
    pmin=Auto.xmin;pmax=Auto.xmax;qmin=Auto.ymin;qmax=Auto.ymax;
  }
  else{pmin=Auto.rl0;pmax=Auto.rl1;qmin=-1e30;qmax=1e30;}
  if(p1s<pmin||p1s>pmax||p2s<qmin||p2s>qmax){
    err_msg("The point lies outside the box (set the Two par axes, or Par Min/Max, around it)");
    goto restore;
  }
  memset(&P,0,sizeof(P));
  P.rhs=xp_rhs;P.setpar=xp_setpar;P.getpar=xp_getpar;P.poll=xp_poll;P.pollctx=NULL;P.ctx=NULL;P.ip=Auto.icp1;
  if(!ddepo_alloc(&P,n,sp->ntst,sp->ncol)){err_msg("Could not set up the periodic orbit solver");goto restore;}
  memcpy(P.X,sp->X,sizeof(double)*(size_t)P.nmesh*n);
  P.T=sp->T;
  dsmin=fabs(Auto.dsmin)<1e-6?fabs(Auto.dsmin):1e-6;
  memset(&R,0,sizeof(R));
  R.P=&P;R.ibr=next_branch();R.lab0=max_label();R.npr=Auto.npr;R.startitp=(kind==0)?5:(kind==1?7:8);
  redraw_auto_menus();
  np=(kind==0?ddepo_fold_curve:(kind==1?ddepo_pd_curve:ddepo_tr_curve))(&P,Auto.icp2,pmin,pmax,qmin,qmax,Auto.ds>=0.0?1:-1,fabs(Auto.ds),dsmin,fabs(Auto.dsmax),
		      Auto.nmx,xp_report2,&R);
  for(i=0;i<NAutoPar&&i<8;i++)constants[Auto_index_to_array[i]]=keep[i];
  evaluate_derived();
  redo_all_fun_tables();
  refreshdisplay();
  if(np==-2)err_msg(kind==0?"Newton did not converge on the fold itself: the LP point may be inaccurate (try a larger NTST) or the second parameter has no effect":(kind==1?"Newton did not converge on the period doubling itself: the PD point may be inaccurate (try a larger NTST) or the second parameter has no effect":"Newton did not converge on the torus point itself: the TR point may be inaccurate (try a larger NTST) or the second parameter has no effect"));
  else if(np==-3||np==-4)err_msg(kind==0?"Could not start the fold curve (no regular tangent or null vector at this point)":(kind==1?"Could not start the period-doubling curve (no multiplier near -1, or no regular tangent)":"Could not start the torus curve (no complex multiplier near the unit circle, or no regular tangent)"));
  else if(np<0)err_msg("The curve could not be set up");
  else{
    iplot=(Auto.plot==XP_TWOPAR);
    printf("%s from label %d: %d points in (%s,%s), from (%g,%g) to (%g,%g)\n",kind==0?"Fold curve of cycles":(kind==1?"Period-doubling curve of cycles":"Torus curve of cycles"),grabpt.lab,np,
	   upar_names[AutoPar[Auto.icp1]],upar_names[AutoPar[Auto.icp2]],R.par0,R.par20,R.lastpar,R.lastpar2);
    printf("  ended: %s\n",(P.stop>=0&&P.stop<=5&&why[P.stop][0])?why[P.stop]:"stopped");
    if(kind==2){
      int q;
      printf("  rotation angle theta: from %.4f to %.4f (continuous, not reduced modulo 2 pi)\n",R.th0,P.theta);
      for(q=0;q<R.nres;q++)
	printf("  strong resonance %s (%s) near (%s,%s)=(%.5g,%.5g)\n",
	       R.res[q].kind==1?"R1":(R.res[q].kind==2?"R2":(R.res[q].kind==3?"R3":"R4")),
	       R.res[q].kind==1?"multiplier +1 twice":(R.res[q].kind==2?"multiplier -1 twice":(R.res[q].kind==3?"theta = 2pi/3":"theta = pi/2")),
	       upar_names[AutoPar[Auto.icp1]],upar_names[AutoPar[Auto.icp2]],R.res[q].p1,R.res[q].p2);
      if(P.stop==3&&fabs(sin(P.theta))<0.05)
	printf("  the curve ended near theta = 0 or pi (a strong resonance, where Newton fails)\n");
    }
    if(!iplot)printf("  (choose Axes -> Two par to see the curve)\n");
    ping();
  }
  ddepo_free(&P);
  return 1;
restore:
  for(i=0;i<NAutoPar&&i<8;i++)constants[Auto_index_to_array[i]]=keep[i];
  evaluate_derived();
  redo_all_fun_tables();
  return 0;
}

int dde_po_fold2_run(void){return curve2_run(0);}
int dde_po_pd2_run(void){return curve2_run(1);}
int dde_po_tr2_run(void){return curve2_run(2);}
