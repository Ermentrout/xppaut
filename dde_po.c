/*  Periodic orbits of delay differential equations: collocation + Newton +
    pseudo-arclength continuation.  See dde_po.h.

    Unknowns   X (nmesh*n values), T (period), and, in a continuation, the
               parameter p.
    Equations  ntst*ncol*n collocation equations
                   dx/ds - T f(x(s), x(s-tau/T), ...; p) = 0
               one integral phase condition  int x . xref' ds = 0
               (and the arclength equation in a continuation).

    The Jacobian is assembled from the n x n derivatives of f with respect
    to the current state and to each delayed state (finite differences of
    the model, through the delay hook), times the Lagrange weights of the
    collocation point and of the point where the delayed value is read.  The
    columns for T and p are finite differences of the residual.  Dense LU.
*/
#include "dde_po.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DD_MAXSLOT 16
#define DD_MAXUSED 64

/* ---------------- the delay hook ---------------- */
static const double *H_X=NULL;
static const DDEPO *H_P=NULL;
static double H_s=0.0,H_T=1.0;
static int H_nseen=0;
static double H_tau[DD_MAXSLOT];
static int H_record=0,H_nused=0;
static int H_uvar[DD_MAXUSED],H_uslot[DD_MAXUSED];
static double H_uval[DD_MAXUSED];
static int H_pvar=-1,H_pslot=-1;
static double H_ph=0.0;
static int ASM_CENTRAL=0;      /* assemble(): central differences for the derivative blocks (accurate, slower) */

static double lag(int m,int k,double s)
{
  double v=1.0;
  int q;
  for(q=0;q<=m;q++)if(q!=k)v*=(s-(double)q/m)/(((double)k-q)/m);
  return v;
}

static double lagd(int m,int k,double s)
{
  double sum=0.0;
  int q,r;
  for(r=0;r<=m;r++){
    double t;
    if(r==k)continue;
    t=1.0/(((double)k-r)/m);
    for(q=0;q<=m;q++)if(q!=k&&q!=r)t*=(s-(double)q/m)/(((double)k-q)/m);
    sum+=t;
  }
  return sum;
}

/* periodic interpolant of component c at s (any real s) for node values X */
static double interp(const DDEPO *P,const double *X,int c,double s)
{
  int j,k,m=P->ncol,n=P->n;
  double sg,v=0.0;
  s-=floor(s);
  sg=s*P->ntst;
  j=(int)floor(sg);
  if(j>=P->ntst)j=P->ntst-1;
  if(j<0)j=0;
  sg-=j;
  for(k=0;k<=m;k++){
    int g=(j*m+k)%P->nmesh;
    v+=lag(m,k,sg)*X[g*n+c];
  }
  return v;
}

double ddepo_eval(const DDEPO *P,int c,double s)
{
  return interp(P,P->X,c,s);
}

double ddepo_delay(double tau,int var)
{
  int q,slot;
  double val;
  if(H_P==NULL||var<1||var>H_P->n)return 0.0;
  for(q=0;q<H_nseen;q++)
    if(fabs(tau-H_tau[q])<=1e-12*(1.0+fabs(tau)))break;
  if(q==H_nseen){
    if(H_nseen<DD_MAXSLOT)H_tau[H_nseen++]=tau;
    else q=DD_MAXSLOT-1;
  }
  slot=q;
  val=interp(H_P,H_X,var-1,H_s-tau/H_T);
  if(H_record){
    int u;
    for(u=0;u<H_nused;u++)if(H_uvar[u]==var-1&&H_uslot[u]==slot)break;
    if(u==H_nused&&H_nused<DD_MAXUSED){
      H_uvar[u]=var-1;H_uslot[u]=slot;H_uval[u]=val;H_nused++;
    }
  }
  if(var-1==H_pvar&&slot==H_pslot)val+=H_ph;
  return val;
}

static void call_rhs(DDEPO *P,const double *X,double s,double T,const double *x,double *f)
{
  H_P=P;H_X=X;H_s=s;H_T=T;H_nseen=0;
  P->rhs(P->ctx,s*T,x,f);
}

/* ---------------- setup ---------------- */
static void gauss01(int m,double *z,double *w)
{
  int i,k,it;
  for(i=0;i<m;i++){
    double x=cos(M_PI*(i+0.75)/(m+0.5)),dp=1.0,p0,p1,p2,dx;
    for(it=0;it<100;it++){
      p0=1.0;p1=x;
      for(k=2;k<=m;k++){p2=((2.0*k-1.0)*x*p1-(k-1.0)*p0)/k;p0=p1;p1=p2;}
      dp=m*(x*p1-p0)/(x*x-1.0);
      dx=p1/dp;
      x-=dx;
      if(fabs(dx)<1e-15)break;
    }
    p0=1.0;p1=x;
    for(k=2;k<=m;k++){p2=((2.0*k-1.0)*x*p1-(k-1.0)*p0)/k;p0=p1;p1=p2;}
    dp=m*(x*p1-p0)/(x*x-1.0);
    z[i]=0.5*(1.0-x);
    w[i]=1.0/((1.0-x*x)*dp*dp);
  }
}

int ddepo_alloc(DDEPO *P,int n,int ntst,int ncol)
{
  int i,k,m=ncol;
  if(n>64||ncol<1||ntst<2)return 0;
  P->n=n;P->ntst=ntst;P->ncol=ncol;
  P->nmesh=ntst*ncol;
  P->hopf=0;P->phi=NULL;
  if(P->tol<=0.0)P->tol=1e-8;
  if(P->amin==0.0)P->amin=1e-3;
  if(P->dtol==0.0)P->dtol=3e-4;
  if(P->maxunk<=0)P->maxunk=3000;
  if(P->maxit<=0)P->maxit=8;
  P->X=(double *)calloc((size_t)P->nmesh*n,sizeof(double));
  P->Xd=(double *)calloc((size_t)P->nmesh*n,sizeof(double));
  P->gz=(double *)malloc(sizeof(double)*m);
  P->gw=(double *)malloc(sizeof(double)*m);
  P->L=(double *)malloc(sizeof(double)*m*(m+1));
  P->Lp=(double *)malloc(sizeof(double)*m*(m+1));
  if(!P->X||!P->Xd||!P->gz||!P->gw||!P->L||!P->Lp){ddepo_free(P);return 0;}
  gauss01(m,P->gz,P->gw);
  for(i=0;i<m;i++)
    for(k=0;k<=m;k++){
      P->L[i*(m+1)+k]=lag(m,k,P->gz[i]);
      P->Lp[i*(m+1)+k]=lagd(m,k,P->gz[i]);
    }
  return 1;
}

void ddepo_free(DDEPO *P)
{
  free(P->X);free(P->Xd);free(P->gz);free(P->gw);free(P->L);free(P->Lp);free(P->phi);
  P->X=P->Xd=P->gz=P->gw=P->L=P->Lp=P->phi=NULL;
  P->hopf=0;
}

void ddepo_set_guess(DDEPO *P,double T,void (*orbit)(void *,double,double *),void *octx)
{
  int j,k,m=P->ncol,n=P->n;
  for(j=0;j<P->ntst;j++)
    for(k=0;k<m;k++)
      orbit(octx,((double)j+(double)k/m)/P->ntst,P->X+(j*m+k)*n);
  P->T=T;
}

/* reference derivative at the collocation points, for the phase condition */
static void update_ref(DDEPO *P)
{
  int j,i,k,b,m=P->ncol,n=P->n;
  for(j=0;j<P->ntst;j++)
    for(i=0;i<m;i++)
      for(b=0;b<n;b++){
	double s=0.0;
	for(k=0;k<=m;k++)s+=P->Lp[i*(m+1)+k]*P->X[((j*m+k)%P->nmesh)*n+b];
	P->Xd[(j*m+i)*n+b]=s*P->ntst;
      }
}

/* ---------------- residual and Jacobian ---------------- */
/* collocation residual, nmesh*n entries */
static void colloc_residual(DDEPO *P,const double *X,double T,double *R)
{
  int j,i,k,a,m=P->ncol,n=P->n;
  double xc[64],f[64];
  for(j=0;j<P->ntst;j++)
    for(i=0;i<m;i++){
      double s=((double)j+P->gz[i])/P->ntst;
      int row=(j*m+i)*n;
      for(a=0;a<n;a++){
	double v=0.0;
	for(k=0;k<=m;k++)v+=P->L[i*(m+1)+k]*X[((j*m+k)%P->nmesh)*n+a];
	xc[a]=v;
      }
      H_record=0;H_pvar=-1;
      call_rhs(P,X,s,T,xc,f);
      for(a=0;a<n;a++){
	double dx=0.0;
	for(k=0;k<=m;k++)dx+=P->Lp[i*(m+1)+k]*X[((j*m+k)%P->nmesh)*n+a];
	R[row+a]=dx*P->ntst-T*f[a];
      }
    }
}

static double phase_residual(DDEPO *P,const double *X)
{
  int j,i,k,a,m=P->ncol,n=P->n;
  double sum=0.0;
  for(j=0;j<P->ntst;j++)
    for(i=0;i<m;i++){
      double wq=P->gw[i]/P->ntst;
      for(a=0;a<n;a++){
	double v=0.0;
	for(k=0;k<=m;k++)v+=P->L[i*(m+1)+k]*X[((j*m+k)%P->nmesh)*n+a];
	sum+=wq*v*P->Xd[(j*m+i)*n+a];
      }
    }
  return sum;
}

double ddepo_residual(DDEPO *P)
{
  double *R=(double *)malloc(sizeof(double)*P->nmesh*P->n);
  double mx=0.0;
  int i;
  if(R==NULL)return -1.0;
  colloc_residual(P,P->X,P->T,R);
  for(i=0;i<P->nmesh*P->n;i++)if(fabs(R[i])>mx)mx=fabs(R[i]);
  free(R);
  return mx;
}

/* Rows 0..N0-1 collocation, row N0 phase condition; columns 0..N0-1 X,
   N0 for T, N0+1 for the parameter.  M has leading dimension ld=N0+2.
   R gets the N0+1 residuals. */
static void assemble(DDEPO *P,const double *X,double T,double par,double *M,int ld,double *R)
{
  int j,i,k,a,b,u,m=P->ncol,n=P->n,N0=P->nmesh*n;
  double xc[64],f[64],f1[64],xp[64];
  double *R1;
  memset(M,0,sizeof(double)*(size_t)ld*(N0+1));
  for(j=0;j<P->ntst;j++)
    for(i=0;i<m;i++){
      double s=((double)j+P->gz[i])/P->ntst;
      int row=(j*m+i)*n,nused,slots;
      int uvar[DD_MAXUSED],uslot[DD_MAXUSED];
      double uval[DD_MAXUSED],tau[DD_MAXSLOT];
      for(a=0;a<n;a++){
	double v=0.0;
	for(k=0;k<=m;k++)v+=P->L[i*(m+1)+k]*X[((j*m+k)%P->nmesh)*n+a];
	xc[a]=v;
      }
      /* base evaluation, recording which delayed values are used */
      H_record=1;H_nused=0;H_pvar=-1;H_pslot=-1;H_ph=0.0;
      call_rhs(P,X,s,T,xc,f);
      H_record=0;
      nused=H_nused;slots=H_nseen;
      for(u=0;u<nused;u++){uvar[u]=H_uvar[u];uslot[u]=H_uslot[u];uval[u]=H_uval[u];}
      for(u=0;u<slots;u++)tau[u]=H_tau[u];
      /* derivative terms and residual */
      for(a=0;a<n;a++){
	double dx=0.0;
	for(k=0;k<=m;k++){
	  int g=(j*m+k)%P->nmesh;
	  double lp=P->Lp[i*(m+1)+k]*P->ntst;
	  dx+=lp*X[g*n+a];
	  M[(size_t)(row+a)*ld+g*n+a]+=lp;
	}
	R[row+a]=dx-T*f[a];
      }
      /* derivative of f wrt the current state */
      for(b=0;b<n;b++){
	double h=1e-6*(1.0+fabs(xc[b])),f2[64];
	memcpy(xp,xc,sizeof(double)*n);
	if(ASM_CENTRAL)h=1e-5*(1.0+fabs(xc[b]));
	xp[b]+=h;
	call_rhs(P,X,s,T,xp,f1);
	if(ASM_CENTRAL){xp[b]-=2.0*h;call_rhs(P,X,s,T,xp,f2);}
	for(a=0;a<n;a++){
	  double A=ASM_CENTRAL?(f1[a]-f2[a])/(2.0*h):(f1[a]-f[a])/h;
	  for(k=0;k<=m;k++){
	    int g=(j*m+k)%P->nmesh;
	    M[(size_t)(row+a)*ld+g*n+b]-=T*A*P->L[i*(m+1)+k];
	  }
	}
      }
      /* derivative of f wrt each delayed value that is used */
      for(u=0;u<nused;u++){
	double h=1e-6*(1.0+fabs(uval[u])),sd,sg,f2[64];
	int jd,kd;
	if(ASM_CENTRAL)h=1e-5*(1.0+fabs(uval[u]));
	H_pvar=uvar[u];H_pslot=uslot[u];H_ph=h;
	call_rhs(P,X,s,T,xc,f1);
	if(ASM_CENTRAL){H_ph=-h;call_rhs(P,X,s,T,xc,f2);}
	H_pvar=-1;H_ph=0.0;
	sd=s-tau[uslot[u]]/T;
	sd-=floor(sd);
	sg=sd*P->ntst;
	jd=(int)floor(sg);
	if(jd>=P->ntst)jd=P->ntst-1;
	if(jd<0)jd=0;
	sg-=jd;
	for(a=0;a<n;a++){
	  double B=ASM_CENTRAL?(f1[a]-f2[a])/(2.0*h):(f1[a]-f[a])/h;
	  for(kd=0;kd<=m;kd++){
	    int g=(jd*m+kd)%P->nmesh;
	    M[(size_t)(row+a)*ld+g*n+uvar[u]]-=T*B*lag(m,kd,sg);
	  }
	}
      }
    }
  /* phase condition */
  for(j=0;j<P->ntst;j++)
    for(i=0;i<m;i++){
      double wq=P->gw[i]/P->ntst;
      for(a=0;a<n;a++){
	double xd=P->Xd[(j*m+i)*n+a];
	for(k=0;k<=m;k++){
	  int g=(j*m+k)%P->nmesh;
	  M[(size_t)N0*ld+g*n+a]+=wq*P->L[i*(m+1)+k]*xd;
	}
      }
    }
  R[N0]=phase_residual(P,X);
  /* columns for T and the parameter: finite differences of the residual */
  R1=(double *)malloc(sizeof(double)*N0);
  if(R1){
    double dT=1e-6*(1.0+fabs(T)),dp=1e-6*(1.0+fabs(par));
    colloc_residual(P,X,T+dT,R1);
    for(i=0;i<N0;i++)M[(size_t)i*ld+N0]=(R1[i]-R[i])/dT;
    P->setpar(P->ctx,P->ip,par+dp);
    colloc_residual(P,X,T,R1);
    P->setpar(P->ctx,P->ip,par);
    for(i=0;i<N0;i++)M[(size_t)i*ld+N0+1]=(R1[i]-R[i])/dp;
    free(R1);
  }
}

/* ---------------- error indicator ---------------- */
/* Residual of  dx/ds = T f  at points that are NOT collocation points (the middles
   between equally spaced nodes), relative to the local size of dx/ds (largest value
   on this and the two neighbouring intervals, so that zero crossings of dx/ds do not
   count).  Small values mean the mesh resolves the orbit; this works for orbits that
   sit close to zero for a long time and then spike. */
double ddepo_defect(DDEPO *P)
{
  int j,q,k,a,m=P->ncol,n=P->n,nt=P->ntst;
  double xc[64],f[64],gmax=0.0,worst=0.0;
  double *def=(double *)malloc(sizeof(double)*(size_t)nt*m*n);
  double *mag=(double *)malloc(sizeof(double)*(size_t)nt*m*n);
  double *dm=(double *)malloc(sizeof(double)*(size_t)nt*n);
  if(!def||!mag||!dm){free(def);free(mag);free(dm);return -1.0;}
  for(j=0;j<nt;j++){
    for(a=0;a<n;a++)dm[j*n+a]=0.0;
    for(q=0;q<m;q++){
      double sg=((double)q+0.5)/m,s=((double)j+sg)/nt;
      for(a=0;a<n;a++){
	double v=0.0;
	for(k=0;k<=m;k++)v+=lag(m,k,sg)*P->X[((j*m+k)%P->nmesh)*n+a];
	xc[a]=v;
      }
      H_record=0;H_pvar=-1;
      call_rhs(P,P->X,s,P->T,xc,f);
      for(a=0;a<n;a++){
	double dx=0.0,tf=P->T*f[a],mg;
	for(k=0;k<=m;k++)dx+=lagd(m,k,sg)*P->X[((j*m+k)%P->nmesh)*n+a];
	dx*=nt;
	def[(j*m+q)*n+a]=dx-tf;
	mg=fabs(dx)>fabs(tf)?fabs(dx):fabs(tf);
	mag[(j*m+q)*n+a]=mg;
	if(mg>dm[j*n+a])dm[j*n+a]=mg;
	if(mg>gmax)gmax=mg;
      }
    }
  }
  for(j=0;j<nt;j++){
    int jl=(j+nt-1)%nt,jr=(j+1)%nt;
    for(a=0;a<n;a++){
      double sc=dm[j*n+a];
      if(dm[jl*n+a]>sc)sc=dm[jl*n+a];
      if(dm[jr*n+a]>sc)sc=dm[jr*n+a];
      sc+=1e-30*gmax;
      if(sc<=0.0)continue;
      for(q=0;q<m;q++){
	double e=fabs(def[(j*m+q)*n+a])/sc;
	if(e>worst)worst=e;
      }
    }
  }
  free(def);free(mag);free(dm);
  return worst;
}

/* ---------------- dense LU ---------------- */
/* solves A x = b in place (b overwritten); A is destroyed; returns 0 if ok */
static int lu_solve(int N,double *A,int ld,double *b)
{
  int i,j,k;
  for(k=0;k<N;k++){
    int piv=k;
    double mx=fabs(A[(size_t)k*ld+k]),t;
    for(i=k+1;i<N;i++){
      t=fabs(A[(size_t)i*ld+k]);
      if(t>mx){mx=t;piv=i;}
    }
    if(mx<1e-300)return 1;
    if(piv!=k){
      for(j=0;j<N;j++){t=A[(size_t)k*ld+j];A[(size_t)k*ld+j]=A[(size_t)piv*ld+j];A[(size_t)piv*ld+j]=t;}
      t=b[k];b[k]=b[piv];b[piv]=t;
    }
    for(i=k+1;i<N;i++){
      double fct=A[(size_t)i*ld+k]/A[(size_t)k*ld+k];
      if(fct!=0.0){
	for(j=k+1;j<N;j++)A[(size_t)i*ld+j]-=fct*A[(size_t)k*ld+j];
	b[i]-=fct*b[k];
      }
    }
  }
  for(i=N-1;i>=0;i--){
    double s=b[i];
    for(j=i+1;j<N;j++)s-=A[(size_t)i*ld+j]*b[j];
    b[i]=s/A[(size_t)i*ld+i];
  }
  return 0;
}

/* ---------------- start at a Hopf point ---------------- */
double ddepo_hopf_vector(int n,const double *Cr,const double *S,double omega,double *va,double *vb)
{
  int N=2*n,i,j,it,mi=0;
  double *K=(double *)malloc(sizeof(double)*(size_t)N*N);
  double *A=(double *)malloc(sizeof(double)*(size_t)N*N);
  double *x=(double *)malloc(sizeof(double)*N);
  double *y=(double *)malloc(sizeof(double)*N);
  double knorm=0.0,res=-1.0,mx,th,c,s;
  if(!K||!A||!x||!y){free(K);free(A);free(x);free(y);return -1.0;}
  /* Delta = Gr + i Gi, Gr=-Cr, Gi=omega I+S;  (Gr+iGi)(a+ib) = (Gr a-Gi b)+i(Gi a+Gr b) */
  for(i=0;i<n;i++)
    for(j=0;j<n;j++){
      double gr=-Cr[i+j*n],gi=S[i+j*n]+((i==j)?omega:0.0);
      K[(size_t)i*N+j]=gr;       K[(size_t)i*N+n+j]=-gi;
      K[(size_t)(n+i)*N+j]=gi;   K[(size_t)(n+i)*N+n+j]=gr;
      if(fabs(gr)>knorm)knorm=fabs(gr);
      if(fabs(gi)>knorm)knorm=fabs(gi);
    }
  /* the size of Delta is that of its terms: the sum omega I + S can cancel exactly */
  for(i=0;i<n*n;i++){
    if(fabs(Cr[i])>knorm)knorm=fabs(Cr[i]);
    if(fabs(S[i])>knorm)knorm=fabs(S[i]);
  }
  if(fabs(omega)>knorm)knorm=fabs(omega);
  if(knorm<=0.0)knorm=1.0;
  for(i=0;i<N;i++)x[i]=1.0+0.37*i;
  /* inverse iteration: K is (nearly) singular, so the shifted solves blow up along the null space */
  for(it=0;it<5;it++){
    for(i=0;i<N;i++){
      for(j=0;j<N;j++)A[(size_t)i*N+j]=K[(size_t)i*N+j];
      A[(size_t)i*N+i]+=1e-9*knorm;
      y[i]=x[i];
    }
    if(lu_solve(N,A,N,y))break;
    mx=0.0;
    for(i=0;i<N;i++)if(fabs(y[i])>mx)mx=fabs(y[i]);
    if(mx<=0.0||mx!=mx)break;
    for(i=0;i<N;i++)x[i]=y[i]/mx;
  }
  for(i=0;i<n;i++){va[i]=x[i];vb[i]=x[n+i];}
  /* rotate the phase so that the largest component is real and equal to 1 */
  mx=0.0;
  for(i=0;i<n;i++){
    double m=sqrt(va[i]*va[i]+vb[i]*vb[i]);
    if(m>mx){mx=m;mi=i;}
  }
  if(mx>0.0){
    th=atan2(vb[mi],va[mi]);
    c=cos(th);s=sin(th);
    for(i=0;i<n;i++){
      double a=va[i],b=vb[i];
      va[i]=(a*c+b*s)/mx;
      vb[i]=(b*c-a*s)/mx;
    }
  }
  /* residual of the null vector, relative to the size of Delta */
  {
    double r=0.0;
    for(i=0;i<N;i++){
      double sum=0.0;
      for(j=0;j<n;j++)sum+=K[(size_t)i*N+j]*va[j]+K[(size_t)i*N+n+j]*vb[j];
      if(fabs(sum)>r)r=fabs(sum);
    }
    res=r/knorm;
  }
  free(K);free(A);free(x);free(y);
  return res;
}

void ddepo_hopf_guess(DDEPO *P,const double *u0,double omega,const double *va,const double *vb,double eps)
{
  int g,c,n=P->n;
  free(P->phi);
  P->phi=(double *)malloc(sizeof(double)*(size_t)P->nmesh*n);
  if(P->phi==NULL)return;
  for(g=0;g<P->nmesh;g++){
    double s=2.0*M_PI*g/P->nmesh;
    for(c=0;c<n;c++){
      double ph=va[c]*cos(s)-vb[c]*sin(s);
      P->phi[g*n+c]=ph;
      P->X[g*n+c]=u0[c]+eps*ph;
    }
  }
  P->T=2.0*M_PI/omega;
  P->hopf=1;
  update_ref(P);
}

/* ---------------- Newton at a fixed parameter ---------------- */
int ddepo_correct(DDEPO *P)
{
  int N0=P->nmesh*P->n,ld=N0+2,N=N0+1,it,i,j;
  double par=P->getpar(P->ctx,P->ip);
  double *M=(double *)malloc(sizeof(double)*(size_t)ld*ld);
  double *S=(double *)malloc(sizeof(double)*(size_t)N*N);
  double *R=(double *)malloc(sizeof(double)*ld);
  int result=-1;
  if(!M||!S||!R){free(M);free(S);free(R);return -1;}
  update_ref(P);
  for(it=1;it<=P->maxit;it++){
    double dmax=0.0,xmax=0.0;
    assemble(P,P->X,P->T,par,M,ld,R);
    for(i=0;i<N;i++)for(j=0;j<N;j++)S[(size_t)i*N+j]=M[(size_t)i*ld+j];
    for(i=0;i<N;i++)R[i]=-R[i];
    if(lu_solve(N,S,N,R))break;
    for(i=0;i<N0;i++){
      P->X[i]+=R[i];
      if(fabs(R[i])>dmax)dmax=fabs(R[i]);
      if(fabs(P->X[i])>xmax)xmax=fabs(P->X[i]);
    }
    P->T+=R[N0];
    if(P->verbose)printf("  correct it %d  |dX|=%.3e  dT=%.3e  T=%.8f\n",it,dmax,R[N0],P->T);
    if(P->T<=0.0)break;
    if(dmax<P->tol*(1.0+xmax)&&fabs(R[N0])<P->tol*P->T){result=it;break;}
  }
  free(M);free(S);free(R);
  return result;
}

/* ---------------- mesh refinement ---------------- */
/* Interpolate the node values X (and optionally a second array Y of the same
   layout) onto a mesh with twice as many intervals.  Returns 0 on success. */
static int refine_mesh(DDEPO *P,double *Y,double **Ynew)
{
  int n=P->n,m=P->ncol,nt2=2*P->ntst,nm2=nt2*m,g,c;
  double *nX=(double *)malloc(sizeof(double)*(size_t)nm2*n);
  double *nXd=(double *)malloc(sizeof(double)*(size_t)nm2*n);
  double *nY=NULL;
  if(Y)nY=(double *)malloc(sizeof(double)*(size_t)nm2*n);
  if(!nX||!nXd||(Y&&!nY)){free(nX);free(nXd);free(nY);return 1;}
  for(g=0;g<nm2;g++){
    double s=(double)g/nm2;
    for(c=0;c<n;c++){
      nX[g*n+c]=interp(P,P->X,c,s);
      if(Y)nY[g*n+c]=interp(P,Y,c,s);
    }
  }
  free(P->X);free(P->Xd);
  P->X=nX;P->Xd=nXd;P->ntst=nt2;P->nmesh=nm2;
  P->nrefine++;
  if(Ynew)*Ynew=nY;else free(nY);
  update_ref(P);
  return 0;
}

static int can_refine(const DDEPO *P)
{
  return P->dtol>0.0&&2*P->ntst*P->ncol*P->n<=P->maxunk;
}

/* ---------------- pseudo-arclength continuation ---------------- */
/* rms deviation of the orbit from its mean, over all components */
static double rms_dev(const DDEPO *P)
{
  int c,g,n=P->n;
  double s=0.0;
  for(c=0;c<n;c++){
    double mean=0.0;
    for(g=0;g<P->nmesh;g++)mean+=P->X[g*n+c];
    mean/=P->nmesh;
    for(g=0;g<P->nmesh;g++){double d=P->X[g*n+c]-mean;s+=d*d;}
  }
  return sqrt(s/((double)P->nmesh*n));
}

static double wdot(int N0,double T,const double *a,const double *b)
{
  double s=0.0;
  int i;
  for(i=0;i<N0;i++)s+=a[i]*b[i];
  return s/N0+a[N0]*b[N0]/(T*T)+a[N0+1]*b[N0+1];
}

typedef struct { double *M,*S,*R,*w; int Nt; } Work;

static void work_free(Work *W)
{
  free(W->M);free(W->S);free(W->R);free(W->w);
  W->M=W->S=W->R=W->w=NULL;W->Nt=0;
}

static int work_alloc(Work *W,int Nt)
{
  work_free(W);
  W->M=(double *)malloc(sizeof(double)*(size_t)Nt*Nt);
  W->S=(double *)malloc(sizeof(double)*(size_t)Nt*Nt);
  W->R=(double *)malloc(sizeof(double)*Nt);
  W->w=(double *)malloc(sizeof(double)*Nt);
  W->Nt=Nt;
  if(!W->M||!W->S||!W->R||!W->w){work_free(W);return 0;}
  return 1;
}

/* Newton for the full system (collocation, phase, arclength <v,u-up>=0; when
   `last`, the parameter is fixed at pbound instead).  Starts from P->X,P->T,*par. */
static int arc_newton(DDEPO *P,Work *W,const double *v,const double *up,int last,
		      double pbound,double *par,int *iters)
{
  int N0=P->nmesh*P->n,Nt=N0+2,i,j,it;
  double *M=W->M,*S=W->S,*R=W->R,*w=W->w;
  for(it=1;it<=P->maxit;it++){
    double dmax=0.0,xmax=0.0,r=0.0;
    P->setpar(P->ctx,P->ip,*par);
    assemble(P,P->X,P->T,*par,M,Nt,R);
    for(i=0;i<N0;i++)M[(size_t)(Nt-1)*Nt+i]=v[i]/N0;
    M[(size_t)(Nt-1)*Nt+N0]=v[N0]/(P->T*P->T);
    M[(size_t)(Nt-1)*Nt+N0+1]=v[N0+1];
    for(i=0;i<N0;i++)r+=v[i]*(P->X[i]-up[i])/N0;
    r+=v[N0]*(P->T-up[N0])/(P->T*P->T)+v[N0+1]*(*par-up[N0+1]);
    R[Nt-1]=r;
    if(last){
      for(i=0;i<=N0;i++)M[(size_t)(Nt-1)*Nt+i]=0.0;
      M[(size_t)(Nt-1)*Nt+N0+1]=1.0;
      R[Nt-1]=*par-pbound;
    }
    for(i=0;i<Nt;i++)for(j=0;j<Nt;j++)S[(size_t)i*Nt+j]=M[(size_t)i*Nt+j];
    for(i=0;i<Nt;i++)w[i]=-R[i];
    if(lu_solve(Nt,S,Nt,w))return 0;
    for(i=0;i<N0;i++){
      P->X[i]+=w[i];
      if(fabs(w[i])>dmax)dmax=fabs(w[i]);
      if(fabs(P->X[i])>xmax)xmax=fabs(P->X[i]);
    }
    P->T+=w[N0];*par+=w[N0+1];
    if(P->T<=0.0)return 0;
    *iters=it;
    if(dmax<P->tol*(1.0+xmax)&&fabs(w[N0])<P->tol*P->T&&fabs(w[N0+1])<P->tol*(1.0+fabs(*par)))return 1;
  }
  return 0;
}

/* tangent of the solution curve at the current point, oriented like vold (or, if
   vold==NULL, with parameter component of sign dir) */
static int tangent_at(DDEPO *P,Work *W,double par,const double *vold,int dir,double *v)
{
  int N0=P->nmesh*P->n,Nt=N0+2,i,j;
  double *M=W->M,*S=W->S,*R=W->R,*w=W->w;
  P->setpar(P->ctx,P->ip,par);
  assemble(P,P->X,P->T,par,M,Nt,R);
  if(vold==NULL){
    int N=N0+1;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++)S[(size_t)i*N+j]=M[(size_t)i*Nt+j];
      w[i]=-M[(size_t)i*Nt+N0+1]*dir;
    }
    if(lu_solve(N,S,N,w))return 1;
    for(i=0;i<N;i++)v[i]=w[i];
    v[N0+1]=dir;
  }
  else{
    double dp;
    double *vo=(double *)malloc(sizeof(double)*Nt);   /* vold may be the same array as v */
    if(!vo)return 1;
    memcpy(vo,vold,sizeof(double)*Nt);
    vold=vo;
    for(i=0;i<N0;i++)M[(size_t)(Nt-1)*Nt+i]=vold[i]/N0;
    M[(size_t)(Nt-1)*Nt+N0]=vold[N0]/(P->T*P->T);
    M[(size_t)(Nt-1)*Nt+N0+1]=vold[N0+1];
    for(i=0;i<Nt;i++)for(j=0;j<Nt;j++)S[(size_t)i*Nt+j]=M[(size_t)i*Nt+j];
    for(i=0;i<Nt;i++)w[i]=0.0;
    w[Nt-1]=1.0;
    if(lu_solve(Nt,S,Nt,w)){free(vo);return 1;}
    for(i=0;i<Nt;i++)v[i]=w[i];
    dp=wdot(N0,P->T,v,vold);
    if(dp<0.0)for(i=0;i<Nt;i++)v[i]=-v[i];
    free(vo);
  }
  {
    double nv=sqrt(wdot(N0,P->T,v,v));
    for(i=0;i<Nt;i++)v[i]/=nv;
  }
  return 0;
}

/* ---------------- special points: test functions from the multipliers ---------------- */
typedef struct { int nmu,ok; double r[16],i[16],triv; } MUCACHE;

static void cache_save(const DDEPO *P,MUCACHE *c)
{
  c->nmu=P->nmu;c->ok=P->mu_ok;c->triv=P->mu_triv;
  memcpy(c->r,P->mu_r,sizeof(c->r));memcpy(c->i,P->mu_i,sizeof(c->i));
}
static void cache_restore(DDEPO *P,const MUCACHE *c)
{
  P->nmu=c->nmu;P->mu_ok=c->ok;P->mu_triv=c->triv;
  memcpy(P->mu_r,c->r,sizeof(c->r));memcpy(P->mu_i,c->i,sizeof(c->i));
}

/* Computes the multipliers of the current point (kept in P->mu_*) and two test functions:
     tf[0] = Re prod (mu_k + 1)           changes sign when a real multiplier crosses -1 (period doubling),
     tf[1] = Re prod_{i<j} (mu_i mu_j-1)  changes sign when a pair crosses the unit circle (torus;
                                          also for a real pair with product 1, a neutral saddle).
   Returns 2, or 0 if the multipliers could not be computed. */
static int event_tests(DDEPO *P,double *tf)
{
  double mr[16],mi[16],te=0.0,cr,ci,t;
  int cnt,i,j,nm;
  P->mu_ok=0;
  cnt=ddepo_floquet(P,mr,mi,12,&te);
  if(cnt<2)return 0;
  for(i=0;i<cnt;i++){P->mu_r[i]=mr[i];P->mu_i[i]=mi[i];}
  P->nmu=cnt;P->mu_triv=te;P->mu_ok=1;
  nm=cnt<8?cnt:8;
  /* do not cut a complex conjugate pair in two */
  if(nm<cnt&&fabs(mi[nm-1])>1e-12&&fabs(mr[nm]-mr[nm-1])<=1e-9*(1.0+fabs(mr[nm]))&&fabs(mi[nm]+mi[nm-1])<=1e-9*(1.0+fabs(mi[nm])))nm++;
  cr=1.0;ci=0.0;
  for(i=0;i<nm;i++){t=cr*(mr[i]+1.0)-ci*mi[i];ci=cr*mi[i]+ci*(mr[i]+1.0);cr=t;}
  tf[0]=cr;
  cr=1.0;ci=0.0;
  for(i=0;i<nm;i++)
    for(j=i+1;j<nm;j++){
      double pr=mr[i]*mr[j]-mi[i]*mi[j]-1.0,pi=mr[i]*mi[j]+mi[i]*mr[j];
      t=cr*pr-ci*pi;ci=cr*pi+ci*pr;cr=t;
    }
  tf[1]=cr;
  return 2;
}

/* Finds the zero of test function idx on the arc between the previous point (uold) and the
   current one, by regula falsi (Illinois) on the arclength fraction, solving the corrector at
   each trial.  On success P holds the special point and *parsp its parameter value. */
static int locate_event(DDEPO *P,Work *W,const double *vprev,const double *uold,double dsused,
			int idx,double tf0,double tf1,double *parsp)
{
  int N0=P->nmesh*P->n,Nt=N0+2,it,i,side=0,found=0;
  double lo=0.0,hi=1.0,flo=tf0,fhi=tf1,th,f,par=0.0,tf[2];
  double *up=(double *)malloc(sizeof(double)*Nt);
  if(!up)return 0;
  for(it=0;it<40;it++){
    int iters=0;
    th=(flo*hi-fhi*lo)/(flo-fhi);
    if(!(th>lo+1e-3*(hi-lo)&&th<hi-1e-3*(hi-lo)))th=0.5*(lo+hi);
    for(i=0;i<Nt;i++)up[i]=uold[i]+th*dsused*vprev[i];
    memcpy(P->X,up,sizeof(double)*N0);P->T=up[N0];par=up[N0+1];
    if(!arc_newton(P,W,vprev,up,0,0.0,&par,&iters))break;
    P->setpar(P->ctx,P->ip,par);
    if(event_tests(P,tf)!=2)break;
    f=tf[idx];
    if(f*fhi>0.0){hi=th;fhi=f;if(side==-1)flo*=0.5;side=-1;}
    else{lo=th;flo=f;if(side==1)fhi*=0.5;side=1;}
    if(fabs(f)<1e-9*(fabs(tf0)+fabs(tf1))||hi-lo<1e-10){found=1;break;}
  }
  free(up);
  if(found)*parsp=par;
  return found;
}

/* Fold of cycles: the parameter component of the tangent changes sign between the previous point
   (uold, tangent vprev) and the current one.  Regula falsi (Illinois) on the arclength fraction, with
   the corrector and the tangent computed at each trial.  On success P holds the fold and *parsp its parameter. */
static int locate_lp(DDEPO *P,Work *W,const double *vprev,const double *uold,double dsused,
		     double f0,double f1,double *parsp)
{
  int N0=P->nmesh*P->n,Nt=N0+2,it,i,side=0,found=0;
  double lo=0.0,hi=1.0,flo=f0,fhi=f1,th,f,par=0.0;
  double *up=(double *)malloc(sizeof(double)*Nt),*vt=(double *)malloc(sizeof(double)*Nt);
  if(!up||!vt){free(up);free(vt);return 0;}
  ASM_CENTRAL=1;      /* the tangent's parameter component must be accurate near the fold */
  for(it=0;it<40;it++){
    int iters=0;
    th=(flo*hi-fhi*lo)/(flo-fhi);
    if(!(th>lo+1e-3*(hi-lo)&&th<hi-1e-3*(hi-lo)))th=0.5*(lo+hi);
    for(i=0;i<Nt;i++)up[i]=uold[i]+th*dsused*vprev[i];
    memcpy(P->X,up,sizeof(double)*N0);P->T=up[N0];par=up[N0+1];
    if(!arc_newton(P,W,vprev,up,0,0.0,&par,&iters))break;
    P->setpar(P->ctx,P->ip,par);
    if(tangent_at(P,W,par,vprev,0,vt))break;
    f=vt[N0+1];
    if(f*fhi>0.0){hi=th;fhi=f;if(side==-1)flo*=0.5;side=-1;}
    else{lo=th;flo=f;if(side==1)fhi*=0.5;side=1;}
    if(fabs(f)<1e-10||hi-lo<1e-11){found=1;break;}
  }
  ASM_CENTRAL=0;
  free(up);free(vt);
  if(found)*parsp=par;
  return found;
}

int ddepo_continue(DDEPO *P,double pmin,double pmax,int dir,double ds,
		   double dsmin,double dsmax,int nmax,DDEPO_REPORT cb,void *rctx)
{
  int N0,Nt,i,step,npts=0,status=0;
  double par,vold_par,dev,devmax;
  double *u=NULL,*v=NULL,*up=NULL,*uold=NULL;
  Work W;
  memset(&W,0,sizeof(W));

  int hopf=P->hopf;
  double *vh=NULL,*vprev=NULL;
  double tfprev[2]={0.0,0.0},tfc[2]={0.0,0.0},dsused=0.0;
  int have_tf=0;
  P->stop=0;P->nunres=0;P->maxdefect=0.0;P->event=0;
  if(P->ntst*P->ncol*P->n>P->maxunk){npts=-1;goto done;}
  if(hopf){
    /* first-order orbit from a Hopf point: correct it with the amplitude along the mode shape
       fixed, the period and the parameter free (the system is not singular for eps>0) */
    int itn=0;
    N0=P->nmesh*P->n;Nt=N0+2;
    par=P->getpar(P->ctx,P->ip);
    u=(double *)malloc(sizeof(double)*Nt);v=(double *)malloc(sizeof(double)*Nt);
    up=(double *)malloc(sizeof(double)*Nt);uold=(double *)malloc(sizeof(double)*Nt);
    vh=(double *)malloc(sizeof(double)*Nt);
    if(!u||!v||!up||!uold||!vh||!P->phi||!work_alloc(&W,Nt)){npts=-1;goto done;}
    memcpy(vh,P->phi,sizeof(double)*N0);vh[N0]=0.0;vh[N0+1]=0.0;
    memcpy(up,P->X,sizeof(double)*N0);up[N0]=P->T;up[N0+1]=par;
    if(!arc_newton(P,&W,vh,up,0,0.0,&par,&itn)){npts=-2;goto done;}
    P->setpar(P->ctx,P->ip,par);
  }
  else{
    if(ddepo_correct(P)<0){npts=-2;goto done;}
    /* refine the starting orbit until the mesh resolves it */
    while(P->dtol>0.0&&ddepo_defect(P)>P->dtol&&can_refine(P)){
      if(refine_mesh(P,NULL,NULL)){npts=-1;goto done;}
      if(ddepo_correct(P)<0){npts=-2;goto done;}
    }
    par=P->getpar(P->ctx,P->ip);
    N0=P->nmesh*P->n;Nt=N0+2;
    u=(double *)malloc(sizeof(double)*Nt);v=(double *)malloc(sizeof(double)*Nt);
    up=(double *)malloc(sizeof(double)*Nt);uold=(double *)malloc(sizeof(double)*Nt);
    if(!u||!v||!up||!uold||!work_alloc(&W,Nt)){npts=-1;goto done;}
  }
  if(par<pmin||par>pmax){npts=-5;goto done;}      /* the start must lie inside the range */
  if(P->locate){
    int n2=event_tests(P,tfprev);
    have_tf=(n2==2);
    /* a start that lies on a special point (an Extend from a PD point) must not find it again */
    if(have_tf){
      if(fabs(tfprev[0])<1e-7)tfprev[0]=0.0;
      if(fabs(tfprev[1])<1e-7)tfprev[1]=0.0;
    }
  }
  if(cb&&cb(rctx,0,par,P->T,P->X,P->nrefine>0?8:0)){npts=1;goto done;}
  npts=1;

  update_ref(P);
  if(hopf){
    /* orient the tangent so that the amplitude along the mode shape grows */
    if(tangent_at(P,&W,par,vh,0,v)){npts=-3;goto done;}
    P->hopf=0;
  }
  else if(tangent_at(P,&W,par,NULL,dir,v)){npts=-3;goto done;}
  vold_par=v[N0+1];
  dev=devmax=rms_dev(P);

  for(step=1;step<=nmax;step++){
    int ok=0,last=0,iters=0,refined=0,landed=0,tfc_ok=0,ev_stop=0,user_stop=0;
    double pbound=0.0,dcur=0.0;
    /* never step by more than half the current amplitude: do not jump through a Hopf point */
    if(ds>0.5*dev)ds=(0.5*dev>dsmin)?0.5*dev:dsmin;
    memcpy(uold,P->X,sizeof(double)*N0);uold[N0]=P->T;uold[N0+1]=par;
    free(vprev);vprev=(double *)malloc(sizeof(double)*Nt);
    if(vprev)memcpy(vprev,v,sizeof(double)*Nt);
    while(!ok){
      double dsl=ds,ppred;
      if(P->poll&&P->poll(P->pollctx)){user_stop=1;break;}
      last=0;
      ppred=uold[N0+1]+dsl*v[N0+1];
      if(v[N0+1]>0.0&&ppred>=pmax){dsl*=(pmax-uold[N0+1])/(ppred-uold[N0+1]);last=1;pbound=pmax;}
      if(v[N0+1]<0.0&&ppred<=pmin){dsl*=(pmin-uold[N0+1])/(ppred-uold[N0+1]);last=1;pbound=pmin;}
      for(i=0;i<Nt;i++)up[i]=uold[i]+dsl*v[i];
      memcpy(P->X,up,sizeof(double)*N0);P->T=up[N0];par=up[N0+1];
      ok=arc_newton(P,&W,v,up,last,pbound,&par,&iters);
      dsused=dsl;
      if(ok){
	/* a converged point much farther from the predictor than the step is a different
	   solution (typically the equilibrium): reject it */
	double dd=0.0,dT=P->T-up[N0],dpp=par-up[N0+1];
	for(i=0;i<N0;i++){double e=P->X[i]-up[i];dd+=e*e;}
	dd=dd/N0+dT*dT/(P->T*P->T)+dpp*dpp;
	if(sqrt(dd)>1.0*dsl+1e-12)ok=0;
      }
      if(!ok){
	ds*=0.5;
	memcpy(P->X,uold,sizeof(double)*N0);P->T=uold[N0];par=uold[N0+1];
	P->setpar(P->ctx,P->ip,par);
	if(ds<dsmin){P->stop=3;npts=npts>0?npts:-4;goto done;}
      }
    }
    /* accepted: new tangent */
    if(tangent_at(P,&W,par,v,0,v)){/* keep the old tangent */}
    /* special points (period doubling, torus) between the previous point and this one */
    if(P->locate){
      int nt2=event_tests(P,tfc);
      tfc_ok=(nt2==2);
      if(tfc_ok&&have_tf&&vprev){
	int ie;
	for(ie=0;ie<2&&!ev_stop;ie++){
	  if(tfprev[ie]*tfc[ie]<0.0){
	    double *Xc=(double *)malloc(sizeof(double)*N0),*vc=(double *)malloc(sizeof(double)*Nt);
	    double Tc=P->T,parc=par,parsp=0.0;
	    MUCACHE cc;
	    if(Xc&&vc){
	      memcpy(Xc,P->X,sizeof(double)*N0);memcpy(vc,v,sizeof(double)*Nt);
	      cache_save(P,&cc);
	      if(locate_event(P,&W,vprev,uold,dsused,ie,tfprev[ie],tfc[ie],&parsp)){
		P->event=ie+1;
		P->setpar(P->ctx,P->ip,parsp);
		if(cb&&cb(rctx,step,parsp,P->T,P->X,16))ev_stop=1;
		P->event=0;
	      }
	      memcpy(P->X,Xc,sizeof(double)*N0);P->T=Tc;par=parc;memcpy(v,vc,sizeof(double)*Nt);
	      cache_restore(P,&cc);
	      P->setpar(P->ctx,P->ip,par);
	    }
	    free(Xc);free(vc);
	  }
	}
      }
      /* a fold of cycles: the parameter component of the tangent changes sign */
      if(!ev_stop&&vprev&&vprev[N0+1]*v[N0+1]<0.0){
	double *Xc=(double *)malloc(sizeof(double)*N0),*vc=(double *)malloc(sizeof(double)*Nt);
	double Tc=P->T,parc=par,parsp=0.0,tmp[2];
	MUCACHE cc;
	if(Xc&&vc){
	  memcpy(Xc,P->X,sizeof(double)*N0);memcpy(vc,v,sizeof(double)*Nt);
	  cache_save(P,&cc);
	  if(locate_lp(P,&W,vprev,uold,dsused,vprev[N0+1],v[N0+1],&parsp)){
	    P->event=3;
	    P->setpar(P->ctx,P->ip,parsp);
	    event_tests(P,tmp);                 /* the multipliers of the fold point itself */
	    if(cb&&cb(rctx,step,parsp,P->T,P->X,16))ev_stop=1;
	    P->event=0;
	  }
	  memcpy(P->X,Xc,sizeof(double)*N0);P->T=Tc;par=parc;memcpy(v,vc,sizeof(double)*Nt);
	  cache_restore(P,&cc);
	  P->setpar(P->ctx,P->ip,par);
	}
	free(Xc);free(vc);
      }
      if(ev_stop){P->stop=4;break;}
    }
    /* does the mesh resolve this orbit?  If not, refine it and converge again on the finer mesh */
    dcur=(P->dtol>0.0)?ddepo_defect(P):0.0;
    while(P->dtol>0.0&&dcur>P->dtol&&can_refine(P)){
      double *vX=NULL,*vnew=NULL;
      double *upn=NULL;
      int Nt2,N02;
      /* the tangent's node values live on the same mesh: interpolate them too */
      vX=(double *)malloc(sizeof(double)*N0);
      if(!vX){npts=-1;goto done;}
      memcpy(vX,v,sizeof(double)*N0);
      {
	double vT=v[N0],vP=v[N0+1],Tsave=P->T;
	if(refine_mesh(P,vX,&vnew)){free(vX);npts=-1;goto done;}
	free(vX);
	N02=P->nmesh*P->n;Nt2=N02+2;
	free(u);free(v);free(up);free(uold);
	u=(double *)malloc(sizeof(double)*Nt2);v=(double *)malloc(sizeof(double)*Nt2);
	up=(double *)malloc(sizeof(double)*Nt2);uold=(double *)malloc(sizeof(double)*Nt2);
	upn=(double *)malloc(sizeof(double)*Nt2);
	if(!u||!v||!up||!uold||!upn||!work_alloc(&W,Nt2)){free(vnew);free(upn);npts=-1;goto done;}
	memcpy(v,vnew,sizeof(double)*N02);v[N02]=vT;v[N02+1]=vP;
	free(vnew);
	memcpy(upn,P->X,sizeof(double)*N02);upn[N02]=Tsave;upn[N02+1]=par;
	N0=N02;Nt=Nt2;
	{
	  int it2=0;
	  if(!arc_newton(P,&W,v,upn,last,pbound,&par,&it2)){free(upn);P->stop=3;npts=npts>0?npts:-2;goto done;}
	}
	free(upn);
	if(tangent_at(P,&W,par,v,0,v)){/* keep */}
	refined=1;
	dev=rms_dev(P);
	dcur=ddepo_defect(P);
      }
    }
    if(P->dtol>0.0&&dcur>P->dtol){        /* could not refine further: remember that this point is not resolved */
      P->nunres++;
      if(dcur>P->maxdefect)P->maxdefect=dcur;
    }
    /* a refinement or the corrector can leave the parameter just outside the range: land on the bound */
    if(!last&&(par>pmax||par<pmin)){
      double pb=(par>pmax)?pmax:pmin;
      int it3=0;
      memcpy(up,P->X,sizeof(double)*N0);up[N0]=P->T;up[N0+1]=par;
      if(!arc_newton(P,&W,v,up,1,pb,&par,&it3)){P->stop=3;npts=npts>0?npts:-2;goto done;}
      last=1;landed=1;
      if(tangent_at(P,&W,par,v,0,v)){/* keep */}
    }
    /* the tests of the accepted point, for the next step (and the multipliers for the report) */
    if(P->locate){
      if(refined||landed||!tfc_ok){int n2=event_tests(P,tfprev);have_tf=(n2==2);}
      else{tfprev[0]=tfc[0];tfprev[1]=tfc[1];have_tf=1;}
    }
    status=0;
    if(v[N0+1]*vold_par<0.0)status|=1;
    vold_par=v[N0+1];
    update_ref(P);
    dev=rms_dev(P);
    if(user_stop){
      /* the last accepted point is the end of the run: restore it and report it as the end point */
      memcpy(P->X,uold,sizeof(double)*N0);P->T=uold[N0];par=uold[N0+1];
      P->setpar(P->ctx,P->ip,par);
      P->stop=5;
      if(step>1&&cb)cb(rctx,step-1,par,P->T,P->X,2);
      break;
    }
    if(dev>devmax)devmax=dev;
    if(P->amin>0.0&&dev<P->amin*devmax)status|=(4|2);
    if(last)status|=2;
    if(step==nmax)status|=2;
    if(refined)status|=8;
    npts++;
    if(P->verbose)printf("step %d  par=%.8f  T=%.8f  iters=%d  ds=%.4g  ntst=%d\n",step,par,P->T,iters,ds,P->ntst);
    if(status&4)P->stop=2;
    else if(last)P->stop=0;
    else if(step==nmax)P->stop=1;
    if(cb&&cb(rctx,step,par,P->T,P->X,status)){P->stop=4;break;}
    if(last||(status&4))break;
    /* step size control */
    if(iters<=3){ds*=1.5;if(ds>dsmax)ds=dsmax;}
    else if(iters>=6){ds*=0.7;if(ds<dsmin)ds=dsmin;}
  }
done:
  work_free(&W);
  free(u);free(v);free(up);free(uold);free(vh);free(vprev);
  return npts;
}

/* ===================================================================== */
/*  Floquet multipliers                                                  */
/* ===================================================================== */

/* ---- eigenvalues of a general real matrix: balance, reduce to Hessenberg form, shifted QR
        (after EISPACK balanc/elmhes/hqr) ---- */
#define EA(i,j) a[((i)-1)*n+((j)-1)]
static void eig_balance(int n,double *a)
{
  int last=0,i,j;
  double radix=2.0,sqrdx=4.0;
  while(!last){
    last=1;
    for(i=1;i<=n;i++){
      double r=0.0,c=0.0,g,f,s;
      for(j=1;j<=n;j++)
	if(j!=i){c+=fabs(EA(j,i));r+=fabs(EA(i,j));}
      if(c!=0.0&&r!=0.0){
	g=r/radix;f=1.0;s=c+r;
	while(c<g){f*=radix;c*=sqrdx;}
	g=r*radix;
	while(c>g){f/=radix;c/=sqrdx;}
	if((c+r)/f<0.95*s){
	  last=0;g=1.0/f;
	  for(j=1;j<=n;j++)EA(i,j)*=g;
	  for(j=1;j<=n;j++)EA(j,i)*=f;
	}
      }
    }
  }
}

static void eig_hessenberg(int n,double *a)
{
  int m,j,i;
  double x,y,t;
  for(m=2;m<n;m++){
    x=0.0;i=m;
    for(j=m;j<=n;j++)
      if(fabs(EA(j,m-1))>fabs(x)){x=EA(j,m-1);i=j;}
    if(i!=m){
      for(j=m-1;j<=n;j++){t=EA(i,j);EA(i,j)=EA(m,j);EA(m,j)=t;}
      for(j=1;j<=n;j++){t=EA(j,i);EA(j,i)=EA(j,m);EA(j,m)=t;}
    }
    if(x!=0.0){
      for(i=m+1;i<=n;i++){
	y=EA(i,m-1);
	if(y!=0.0){
	  y/=x;EA(i,m-1)=y;
	  for(j=m;j<=n;j++)EA(i,j)-=y*EA(m,j);
	  for(j=1;j<=n;j++)EA(j,m)+=y*EA(j,i);
	}
      }
    }
  }
  for(i=3;i<=n;i++)for(j=1;j<=i-2;j++)EA(i,j)=0.0;   /* drop the stored multipliers */
}

static double sgn_(double a,double b){return (b>=0.0)?fabs(a):-fabs(a);}

static int eig_hqr(int n,double *a,double *wr,double *wi)
{
  int nn,m,l,k,j,its,i,mmin;
  double z=0,y,x,w,v,u,t,s,r=0,q=0,p=0,anorm=0.0;
  for(i=1;i<=n;i++)
    for(j=(i-1>1?i-1:1);j<=n;j++)anorm+=fabs(EA(i,j));
  nn=n;t=0.0;
  while(nn>=1){
    its=0;
    do{
      for(l=nn;l>=2;l--){
	s=fabs(EA(l-1,l-1))+fabs(EA(l,l));
	if(s==0.0)s=anorm;
	if(fabs(EA(l,l-1))+s==s){EA(l,l-1)=0.0;break;}
      }
      x=EA(nn,nn);
      if(l==nn){
	wr[nn-1]=x+t;wi[nn-1]=0.0;nn--;
      }
      else{
	y=EA(nn-1,nn-1);
	w=EA(nn,nn-1)*EA(nn-1,nn);
	if(l==nn-1){
	  p=0.5*(y-x);q=p*p+w;z=sqrt(fabs(q));x+=t;
	  if(q>=0.0){
	    z=p+sgn_(z,p);
	    wr[nn-2]=wr[nn-1]=x+z;
	    if(z!=0.0)wr[nn-1]=x-w/z;
	    wi[nn-2]=wi[nn-1]=0.0;
	  }
	  else{
	    wr[nn-2]=wr[nn-1]=x+p;
	    wi[nn-2]=-(wi[nn-1]=z);
	  }
	  nn-=2;
	}
	else{
	  if(its==60)return 1;
	  if(its==10||its==20){
	    t+=x;
	    for(i=1;i<=nn;i++)EA(i,i)-=x;
	    s=fabs(EA(nn,nn-1))+fabs(EA(nn-1,nn-2));
	    y=x=0.75*s;
	    w=-0.4375*s*s;
	  }
	  ++its;
	  for(m=nn-2;m>=l;m--){
	    z=EA(m,m);r=x-z;s=y-z;
	    p=(r*s-w)/EA(m+1,m)+EA(m,m+1);
	    q=EA(m+1,m+1)-z-r-s;
	    r=EA(m+2,m+1);
	    s=fabs(p)+fabs(q)+fabs(r);
	    p/=s;q/=s;r/=s;
	    if(m==l)break;
	    u=fabs(EA(m,m-1))*(fabs(q)+fabs(r));
	    v=fabs(p)*(fabs(EA(m-1,m-1))+fabs(z)+fabs(EA(m+1,m+1)));
	    if(u+v==v)break;
	  }
	  for(i=m+2;i<=nn;i++){
	    EA(i,i-2)=0.0;
	    if(i!=m+2)EA(i,i-3)=0.0;
	  }
	  for(k=m;k<=nn-1;k++){
	    if(k!=m){
	      p=EA(k,k-1);q=EA(k+1,k-1);r=0.0;
	      if(k!=nn-1)r=EA(k+2,k-1);
	      if((x=fabs(p)+fabs(q)+fabs(r))!=0.0){p/=x;q/=x;r/=x;}
	    }
	    if((s=sgn_(sqrt(p*p+q*q+r*r),p))!=0.0){
	      if(k==m){
		if(l!=m)EA(k,k-1)=-EA(k,k-1);
	      }
	      else EA(k,k-1)=-s*x;
	      p+=s;x=p/s;y=q/s;z=r/s;q/=p;r/=p;
	      for(j=k;j<=nn;j++){
		p=EA(k,j)+q*EA(k+1,j);
		if(k!=nn-1){p+=r*EA(k+2,j);EA(k+2,j)-=p*z;}
		EA(k+1,j)-=p*y;
		EA(k,j)-=p*x;
	      }
	      mmin=nn<k+3?nn:k+3;
	      for(i=l;i<=mmin;i++){
		p=x*EA(i,k)+y*EA(i,k+1);
		if(k!=nn-1){p+=z*EA(i,k+2);EA(i,k+2)-=p*r;}
		EA(i,k+1)-=p*q;
		EA(i,k)-=p;
	      }
	    }
	  }
	}
      }
    }while(l<nn-1);
  }
  return 0;
}
#undef EA

int ddepo_eigenvalues(int n,double *a,double *wr,double *wi)
{
  if(n<1)return 1;
  if(n==1){wr[0]=a[0];wi[0]=0.0;return 0;}
  eig_balance(n,a);
  eig_hessenberg(n,a);
  return eig_hqr(n,a,wr,wi);
}

/* ---- LU with row pivoting for several right-hand sides ---- */
static int lu_factor(int N,double *A,int *piv)
{
  int i,j,k;
  for(k=0;k<N;k++){
    int p=k;
    double mx=fabs(A[(size_t)k*N+k]),t;
    for(i=k+1;i<N;i++){t=fabs(A[(size_t)i*N+k]);if(t>mx){mx=t;p=i;}}
    if(mx<1e-300)return 1;
    piv[k]=p;
    if(p!=k)for(j=0;j<N;j++){t=A[(size_t)k*N+j];A[(size_t)k*N+j]=A[(size_t)p*N+j];A[(size_t)p*N+j]=t;}
    for(i=k+1;i<N;i++){
      double f=A[(size_t)i*N+k]/A[(size_t)k*N+k];
      A[(size_t)i*N+k]=f;
      for(j=k+1;j<N;j++)A[(size_t)i*N+j]-=f*A[(size_t)k*N+j];
    }
  }
  return 0;
}

/* B is N x nc, row-major; overwritten by A^{-1} B */
static void lu_solve_mat(int N,const double *A,const int *piv,double *B,int nc)
{
  int i,j,k;
  /* the factorization swapped whole rows (multipliers included): apply all the swaps first */
  for(k=0;k<N;k++)
    if(piv[k]!=k)for(j=0;j<nc;j++){double t=B[(size_t)k*nc+j];B[(size_t)k*nc+j]=B[(size_t)piv[k]*nc+j];B[(size_t)piv[k]*nc+j]=t;}
  for(k=0;k<N;k++){
    for(i=k+1;i<N;i++){
      double f=A[(size_t)i*N+k];
      if(f!=0.0)for(j=0;j<nc;j++)B[(size_t)i*nc+j]-=f*B[(size_t)k*nc+j];
    }
  }
  for(i=N-1;i>=0;i--){
    for(k=i+1;k<N;k++){
      double f=A[(size_t)i*N+k];
      if(f!=0.0)for(j=0;j<nc;j++)B[(size_t)i*nc+j]-=f*B[(size_t)k*nc+j];
    }
    {double d=A[(size_t)i*N+i];for(j=0;j<nc;j++)B[(size_t)i*nc+j]/=d;}
  }
}

/* ---- coefficients of the variational equation along the orbit ----
   dy/dt = A(t) y(t) + sum_d B_d(t) y(t - tau_d), at the collocation points (index j*ncol+i) */
static int var_coeffs(DDEPO *P,double *Ac,double *Bc,int *nslot,double *tau)
{
  int j,i,k,a,b,u,m=P->ncol,n=P->n,nm=P->nmesh;
  double xc[64],f[64],f1[64],xp[64];
  *nslot=0;
  for(j=0;j<P->ntst;j++)
    for(i=0;i<m;i++){
      double s=((double)j+P->gz[i])/P->ntst;
      int g=j*m+i,nused,slots;
      int uvar[DD_MAXUSED],uslot[DD_MAXUSED];
      double uval[DD_MAXUSED];
      for(a=0;a<n;a++){
	double v=0.0;
	for(k=0;k<=m;k++)v+=P->L[i*(m+1)+k]*P->X[((j*m+k)%nm)*n+a];
	xc[a]=v;
      }
      H_record=1;H_nused=0;H_pvar=-1;H_pslot=-1;H_ph=0.0;
      call_rhs(P,P->X,s,P->T,xc,f);
      H_record=0;
      nused=H_nused;slots=H_nseen;
      for(u=0;u<nused;u++){uvar[u]=H_uvar[u];uslot[u]=H_uslot[u];uval[u]=H_uval[u];}
      for(u=0;u<slots;u++){
	if(u>=*nslot){tau[u]=H_tau[u];*nslot=u+1;}
      }
      /* central differences: the multipliers (a crossing of -1 or of the unit circle) need accurate
         coefficients, and a one-sided difference of a steep nonlinearity is off by about h f''/2 */
      for(b=0;b<n;b++){
	double h=1e-5*(1.0+fabs(xc[b])),f2[64];
	memcpy(xp,xc,sizeof(double)*n);
	xp[b]+=h;
	call_rhs(P,P->X,s,P->T,xp,f1);
	xp[b]-=2.0*h;
	call_rhs(P,P->X,s,P->T,xp,f2);
	for(a=0;a<n;a++)Ac[((size_t)g*n+a)*n+b]=(f1[a]-f2[a])/(2.0*h);
      }
      for(u=0;u<nused;u++){
	double h=1e-5*(1.0+fabs(uval[u])),f2[64];
	H_pvar=uvar[u];H_pslot=uslot[u];H_ph=h;
	call_rhs(P,P->X,s,P->T,xc,f1);
	H_ph=-h;
	call_rhs(P,P->X,s,P->T,xc,f2);
	H_pvar=-1;H_ph=0.0;
	for(a=0;a<n;a++)
	  Bc[(((size_t)uslot[u]*nm+g)*n+a)*n+uvar[u]]=(f1[a]-f2[a])/(2.0*h);
      }
    }
  return 0;
}

/* ---- monodromy matrix (nh x nh, row-major, malloc'd) ---- */
static int monodromy(DDEPO *P,double **Mh,int *nhout,double **RwOut,int *LOut)
{
  int n=P->n,m=P->ncol,nt=P->ntst,nm=P->nmesh;
  int nslot=0,d,q,i,k,kk,a,b,r,L,nh,G0,Ntot;
  double tau[DD_MAXSLOT],theta,T=P->T,thmax=0.0;
  double *Ac=NULL,*Bc=NULL,*Rw=NULL,*K=NULL,*RHS=NULL,*M=NULL;
  int *piv=NULL;
  size_t rowlen;
  Ac=(double *)calloc((size_t)nm*n*n,sizeof(double));
  Bc=(double *)calloc((size_t)DD_MAXSLOT*nm*n*n,sizeof(double));
  if(!Ac||!Bc){free(Ac);free(Bc);return 1;}
  var_coeffs(P,Ac,Bc,&nslot,tau);
  for(d=0;d<nslot;d++){theta=tau[d]/T;if(theta>thmax)thmax=theta;}
  L=(int)ceil(thmax*nt-1e-9);
  if(L<0)L=0;
  nh=(L*m+1)*n;
  if(nh>1600){free(Ac);free(Bc);return 2;}
  G0=-L*m;
  Ntot=(nt+L)*m+1;
  rowlen=(size_t)nh;
  Rw=(double *)calloc((size_t)Ntot*n*rowlen,sizeof(double));
  K=(double *)malloc(sizeof(double)*(size_t)m*n*m*n);
  RHS=(double *)malloc(sizeof(double)*(size_t)m*n*rowlen);
  piv=(int *)malloc(sizeof(int)*m*n);
  if(!Rw||!K||!RHS||!piv){free(Ac);free(Bc);free(Rw);free(K);free(RHS);free(piv);return 1;}
  for(i=0;i<nh;i++)Rw[(size_t)i*rowlen+i]=1.0;           /* initial state: identity on nodes -L m .. 0 */
  for(q=0;q<nt;q++){
    int j=q%nt,N=m*n;
    memset(K,0,sizeof(double)*(size_t)N*N);
    memset(RHS,0,sizeof(double)*(size_t)N*rowlen);
    for(i=0;i<m;i++){
      int g=j*m+i;
      const double *A=Ac+(size_t)g*n*n;
      /* the derivative and A terms: nodes k=0..m of the current interval; k=0 is known */
      for(k=0;k<=m;k++){
	double lp=nt*P->Lp[i*(m+1)+k],lv=P->L[i*(m+1)+k];
	for(a=0;a<n;a++)for(b=0;b<n;b++){
	  double c=((a==b)?lp:0.0)-T*A[a*n+b]*lv;
	  if(k>=1)K[(size_t)(i*n+a)*N+(k-1)*n+b]+=c;
	  else{
	    int node=q*m+0;             /* the known first node */
	    const double *src=Rw+(size_t)((node-G0)*n+b)*rowlen;
	    double *dst=RHS+(size_t)(i*n+a)*rowlen;
	    for(r=0;r<nh;r++)dst[r]+=c*src[r];
	  }
	}
      }
      /* delayed terms */
      for(d=0;d<nslot;d++){
	double u=q+P->gz[i]-tau[d]/T*nt;
	int qd=(int)floor(u);
	double sg=u-qd;
	const double *B=Bc+(((size_t)d*nm+g)*n)*n;
	for(kk=0;kk<=m;kk++){
	  double ld=lag(m,kk,sg);
	  int node=qd*m+kk;
	  for(a=0;a<n;a++)for(b=0;b<n;b++){
	    double c=-T*B[a*n+b]*ld;
	    if(c==0.0)continue;
	    if(node>=q*m+1){                /* an unknown node of the current interval */
	      K[(size_t)(i*n+a)*N+(node-q*m-1)*n+b]+=c;
	    }
	    else{
	      const double *src=Rw+(size_t)((node-G0)*n+b)*rowlen;
	      double *dst=RHS+(size_t)(i*n+a)*rowlen;
	      for(r=0;r<nh;r++)dst[r]+=c*src[r];
	    }
	  }
	}
      }
    }
    if(lu_factor(N,K,piv)){free(Ac);free(Bc);free(Rw);free(K);free(RHS);free(piv);return 3;}
    for(r=0;r<N*(int)rowlen;r++)RHS[r]=-RHS[r];
    lu_solve_mat(N,K,piv,RHS,(int)rowlen);
    /* the new nodes q m+1 .. q m+m */
    for(k=1;k<=m;k++)
      for(a=0;a<n;a++)
	memcpy(Rw+(size_t)(((q*m+k)-G0)*n+a)*rowlen,RHS+(size_t)((k-1)*n+a)*rowlen,sizeof(double)*rowlen);
  }
  /* the state after one period: nodes (nt-L) m .. nt m */
  M=(double *)malloc(sizeof(double)*(size_t)nh*nh);
  if(!M){free(Ac);free(Bc);free(Rw);free(K);free(RHS);free(piv);return 1;}
  for(i=0;i<nh;i++)
    memcpy(M+(size_t)i*nh,Rw+(size_t)(((nt-L)*m-G0)*n+i)*rowlen,sizeof(double)*rowlen);
  free(Ac);free(Bc);free(K);free(RHS);free(piv);
  if(RwOut){*RwOut=Rw;*LOut=L;}else free(Rw);
  *Mh=M;*nhout=nh;
  return 0;
}

typedef struct { double r,i; } CPLX;
static int cmp_mod(const void *x,const void *y)
{
  const CPLX *a=(const CPLX *)x,*b=(const CPLX *)y;
  double ma=a->r*a->r+a->i*a->i,mb=b->r*b->r+b->i*b->i;
  return (ma<mb)-(ma>mb);
}

int ddepo_multipliers(DDEPO *P,double *mur,double *mui,int nmax)
{
  double *M=NULL,*wr,*wi;
  CPLX *c;
  int nh=0,i,cnt;
  if(monodromy(P,&M,&nh,NULL,NULL))return -1;
  wr=(double *)malloc(sizeof(double)*nh);wi=(double *)malloc(sizeof(double)*nh);
  c=(CPLX *)malloc(sizeof(CPLX)*nh);
  if(!wr||!wi||!c){free(M);free(wr);free(wi);free(c);return -1;}
  if(ddepo_eigenvalues(nh,M,wr,wi)){free(M);free(wr);free(wi);free(c);return -1;}
  for(i=0;i<nh;i++){c[i].r=wr[i];c[i].i=wi[i];}
  qsort(c,nh,sizeof(CPLX),cmp_mod);
  cnt=nh<nmax?nh:nmax;
  for(i=0;i<cnt;i++){mur[i]=c[i].r;mui[i]=c[i].i;}
  free(M);free(wr);free(wi);free(c);
  return cnt;
}

int ddepo_floquet(DDEPO *P,double *mur,double *mui,int nmax,double *triv_err)
{
  double *tr,*ti,best=1e300;
  int cnt,i,bi=-1,k;
  tr=(double *)malloc(sizeof(double)*(nmax+1));ti=(double *)malloc(sizeof(double)*(nmax+1));
  if(!tr||!ti){free(tr);free(ti);return -1;}
  cnt=ddepo_multipliers(P,tr,ti,nmax+1);
  if(cnt<1){free(tr);free(ti);return -1;}
  for(i=0;i<cnt;i++){
    double d=(tr[i]-1.0)*(tr[i]-1.0)+ti[i]*ti[i];
    if(d<best){best=d;bi=i;}
  }
  if(triv_err)*triv_err=sqrt(best);
  k=0;
  for(i=0;i<cnt&&k<nmax;i++){
    if(i==bi)continue;
    mur[k]=tr[i];mui[k]=ti[i];k++;
  }
  free(tr);free(ti);
  return k;
}

/* ===================================================================== */
/*  Branch switching at a period doubling                                */
/* ===================================================================== */

/* The critical Floquet eigenfunction of the current orbit for the real multiplier closest to
   mu_target (-1 for a period doubling).  The eigenvector of the monodromy matrix is found by
   inverse iteration and propagated around the orbit with the stored time-stepping rows.
   v receives nmesh*n values (the nodes g=0..nmesh-1 of one period, scaled to maximum 1);
   *mu is the multiplier and *mismatch measures how well v(1)=mu v(0) holds.  Returns 0 on success. */
int ddepo_floquet_mode(DDEPO *P,double mu_target,double *v,double *mu,double *mismatch)
{
  double *M=NULL,*Rw=NULL,*Mc=NULL,*wr=NULL,*wi=NULL,*A=NULL,*x=NULL;
  int *piv=NULL;
  int nh=0,L=0,i,j,k,g,a,it,n=P->n,m=P->ncol,nt=P->ntst,nm=P->nmesh,G0,ret=1;
  double mu0,best=1e300,mx,vend,vmax=0.0,mm=0.0;
  if(monodromy(P,&M,&nh,&Rw,&L))return 1;
  G0=-L*m;
  Mc=(double *)malloc(sizeof(double)*(size_t)nh*nh);
  wr=(double *)malloc(sizeof(double)*nh);wi=(double *)malloc(sizeof(double)*nh);
  A=(double *)malloc(sizeof(double)*(size_t)nh*nh);
  x=(double *)malloc(sizeof(double)*nh);piv=(int *)malloc(sizeof(int)*nh);
  if(!Mc||!wr||!wi||!A||!x||!piv)goto out;
  memcpy(Mc,M,sizeof(double)*(size_t)nh*nh);
  if(ddepo_eigenvalues(nh,Mc,wr,wi))goto out;
  k=-1;
  for(i=0;i<nh;i++){
    double d=fabs(wr[i]-mu_target)+fabs(wi[i]);
    if(d<best){best=d;k=i;}
  }
  if(k<0)goto out;
  mu0=wr[k];
  /* inverse iteration with a shift very close to the eigenvalue */
  for(i=0;i<nh;i++)for(j=0;j<nh;j++)A[(size_t)i*nh+j]=M[(size_t)i*nh+j];
  for(i=0;i<nh;i++)A[(size_t)i*nh+i]-=mu0*(1.0+1e-9)+1e-13;
  if(lu_factor(nh,A,piv))goto out;
  for(i=0;i<nh;i++)x[i]=1.0+0.3*sin(1.7*i);
  for(it=0;it<5;it++){
    lu_solve_mat(nh,A,piv,x,1);
    mx=0.0;
    for(i=0;i<nh;i++)if(fabs(x[i])>mx)mx=fabs(x[i]);
    if(!(mx>0.0)||mx!=mx)goto out;
    for(i=0;i<nh;i++)x[i]/=mx;
  }
  /* node values over the period: row g of the time stepping applied to the eigenvector */
  for(g=0;g<nm;g++)
    for(a=0;a<n;a++){
      const double *row=Rw+(size_t)(((g-G0)*n+a))*nh;
      double sum=0.0;
      for(i=0;i<nh;i++)sum+=row[i]*x[i];
      v[g*n+a]=sum;
      if(fabs(sum)>vmax)vmax=fabs(sum);
    }
  for(a=0;a<n;a++){
    const double *row=Rw+(size_t)(((nt*m-G0)*n+a))*nh;
    double sum=0.0;
    for(i=0;i<nh;i++)sum+=row[i]*x[i];
    vend=sum;
    if(fabs(vend-mu0*v[a])>mm)mm=fabs(vend-mu0*v[a]);
  }
  if(!(vmax>0.0))goto out;
  for(i=0;i<nm*n;i++)v[i]/=vmax;
  *mu=mu0;*mismatch=mm/vmax;
  ret=0;
out:
  free(M);free(Rw);free(Mc);free(wr);free(wi);free(A);free(x);free(piv);
  return ret;
}

/* Set up P2 (twice as many intervals as P1, same ncol and n) with the first-order orbit of the
   doubled period on the period doubling branch:
        x2(s) = x1(2s mod 1) + eps * sign * v(2s mod 1),   T2 = 2 T1,
   with sign = +1 on the first half and -1 on the second (v is anti-periodic).  The next
   ddepo_continue on P2 corrects it with the amplitude along the mode fixed, as for a Hopf start. */
int ddepo_pd_guess(DDEPO *P2,const DDEPO *P1,const double *v,double eps)
{
  int n=P1->n,nm1=P1->nmesh,g,a;
  if(P2->n!=n||P2->ncol!=P1->ncol||P2->ntst!=2*P1->ntst)return 1;
  free(P2->phi);
  P2->phi=(double *)malloc(sizeof(double)*(size_t)P2->nmesh*n);
  if(!P2->phi)return 1;
  for(g=0;g<P2->nmesh;g++){
    int g1=g%nm1;
    double sg=(g<nm1)?1.0:-1.0;
    for(a=0;a<n;a++){
      double ph=sg*v[g1*n+a];
      P2->phi[g*n+a]=ph;
      P2->X[g*n+a]=P1->X[g1*n+a]+eps*ph;
    }
  }
  P2->T=2.0*P1->T;
  P2->hopf=1;
  update_ref(P2);
  return 0;
}

/* ===================================================================== */
/*  Two-parameter continuation of a fold or a period doubling of cycles  */
/* ===================================================================== */
/* Common unknowns z = (X[N0], T, p1, p2, null vector), the null vector starting at offset yo = N0+3.
   Equations (E1 .. E3), plus one arclength equation:
     E1  F(X,T;p1,p2) = 0                 collocation and phase condition       (N0+1)
     fold             E2  F_u (W,wT) = 0  null vector of the collocation Jacobian  (N0+1)
                      E3  |(W,wT)|^2 = 1                                           (1)
     period doubling  E2  L(X,T;p) Y = 0  antiperiodic variational equation       (N0)
                      E3  |Y|^2 = 1                                                (1)
   F_u (W,wT) is the central difference of F along (W,wT).  The antiperiodic operator is handled on a
   doubled mesh: (X,X) with period 2T and the perturbation (Y,-Y) is periodic there, so the same
   collocation code gives both L Y (a central difference) and the matrix L (the doubled Jacobian, first
   block row, first column block minus second).  The Jacobian of E2 with respect to (X,T,p1,p2) is a
   finite difference of E2. */
typedef struct {
  int N0,M,ip2,kind,yo;               /* kind 0: fold, 1: period doubling */
  double *wt;                         /* weights of the arclength norm */
  double *Xp,*Xm,*F1,*F2,*Mj,*K,*A,*rhs,*Rt,*E2a,*E2b,*zt;
  int *piv;
  int has_p2;                         /* kind 1: the doubled problem */
  DDEPO P2;
  double *X2,*R2,*Mj2,*Rt2;
  double *Yref;                       /* kind 2: reference for the phase gauge of Y (2*N0) */
} FoldW;

static void foldw_free(FoldW *F)
{
  free(F->wt);free(F->Xp);free(F->Xm);free(F->F1);free(F->F2);free(F->Mj);free(F->K);free(F->A);
  free(F->rhs);free(F->Rt);free(F->E2a);free(F->E2b);free(F->zt);free(F->piv);
  free(F->X2);free(F->R2);free(F->Mj2);free(F->Rt2);free(F->Yref);
  if(F->has_p2)ddepo_free(&F->P2);
  memset(F,0,sizeof(*F));
}

static int foldw_alloc(FoldW *F,const DDEPO *P,int N0,int ip2,int kind)
{
  int ne2=(kind==0)?N0+1:((kind==1)?N0:2*N0);
  int yo=N0+3+((kind==2)?1:0);                 /* kind 2 has the rotation angle after p2 */
  int M=yo+ne2;
  memset(F,0,sizeof(*F));
  F->N0=N0;F->M=M;F->ip2=ip2;F->kind=kind;F->yo=yo;
  F->wt=(double *)calloc(M,sizeof(double));
  F->Xp=(double *)malloc(sizeof(double)*N0);F->Xm=(double *)malloc(sizeof(double)*N0);
  F->F1=(double *)malloc(sizeof(double)*(N0+1));F->F2=(double *)malloc(sizeof(double)*(N0+1));
  F->Mj=(double *)malloc(sizeof(double)*(size_t)(N0+1)*(N0+2));
  F->K=(double *)calloc((size_t)(M-1)*M,sizeof(double));
  F->A=(double *)malloc(sizeof(double)*(size_t)M*M);
  F->rhs=(double *)malloc(sizeof(double)*M);F->Rt=(double *)malloc(sizeof(double)*(N0+1));
  F->E2a=(double *)malloc(sizeof(double)*(2*N0+2));F->E2b=(double *)malloc(sizeof(double)*(2*N0+2));
  F->zt=(double *)malloc(sizeof(double)*M);F->piv=(int *)malloc(sizeof(int)*M);
  if(!F->wt||!F->Xp||!F->Xm||!F->F1||!F->F2||!F->Mj||!F->K||!F->A||!F->rhs||!F->Rt||!F->E2a||!F->E2b||!F->zt||!F->piv){
    foldw_free(F);return 0;
  }
  if(kind==2){
    F->Yref=(double *)calloc((size_t)2*N0,sizeof(double));
    if(!F->Yref){foldw_free(F);return 0;}
  }
  if(kind==1){
    F->P2.rhs=P->rhs;F->P2.setpar=P->setpar;F->P2.getpar=P->getpar;F->P2.ctx=P->ctx;F->P2.ip=P->ip;
    if(!ddepo_alloc(&F->P2,P->n,2*P->ntst,P->ncol)){foldw_free(F);return 0;}
    F->has_p2=1;
    F->X2=(double *)malloc(sizeof(double)*2*N0);F->R2=(double *)malloc(sizeof(double)*2*N0);
    F->Mj2=(double *)malloc(sizeof(double)*(size_t)(2*N0+1)*(2*N0+2));F->Rt2=(double *)malloc(sizeof(double)*(2*N0+1));
    if(!F->X2||!F->R2||!F->Mj2||!F->Rt2){foldw_free(F);return 0;}
  }
  return 1;
}

static void full_F(DDEPO *P,const double *X,double T,double *F)
{
  colloc_residual(P,X,T,F);
  F[P->nmesh*P->n]=phase_residual(P,X);
}

static void fold_setpars(DDEPO *P,const FoldW *F,const double *z)
{
  /* setting a parameter can be expensive in the host (derived constants, tables): only if it changed */
  if(P->getpar(P->ctx,P->ip)!=z[F->N0+1])P->setpar(P->ctx,P->ip,z[F->N0+1]);
  if(P->getpar(P->ctx,F->ip2)!=z[F->N0+2])P->setpar(P->ctx,F->ip2,z[F->N0+2]);
}

/* ---- fold ---- */
/* E2 = F_u (W,wT) by a central difference */
static void fold_E2(DDEPO *P,FoldW *F,const double *z,double *out)
{
  int N0=F->N0,i;
  double xm=0.0,eps,Tp,Tm;
  for(i=0;i<N0;i++)if(fabs(z[i])>xm)xm=fabs(z[i]);
  eps=1e-5*(1.0+xm);
  fold_setpars(P,F,z);
  for(i=0;i<N0;i++){F->Xp[i]=z[i]+eps*z[N0+3+i];F->Xm[i]=z[i]-eps*z[N0+3+i];}
  Tp=z[N0]+eps*z[2*N0+3];Tm=z[N0]-eps*z[2*N0+3];
  full_F(P,F->Xp,Tp,F->F1);
  full_F(P,F->Xm,Tm,F->F2);
  for(i=0;i<=N0;i++)out[i]=(F->F1[i]-F->F2[i])/(2.0*eps);
}

/* ---- period doubling ---- */
/* E2 = L Y by the central difference on the doubled mesh: the second half of the residual vector is
   the first half with Y -> -Y, so one evaluation at +eps gives both sides */
static void pd_E2(DDEPO *P,FoldW *F,const double *z,double *out)
{
  int N0=F->N0,i,yo=F->yo;
  double xm=0.0,eps;
  for(i=0;i<N0;i++)if(fabs(z[i])>xm)xm=fabs(z[i]);
  eps=1e-5*(1.0+xm);
  fold_setpars(P,F,z);
  for(i=0;i<N0;i++){F->X2[i]=z[i]+eps*z[yo+i];F->X2[N0+i]=z[i]-eps*z[yo+i];}
  colloc_residual(&F->P2,F->X2,2.0*z[N0],F->R2);
  for(i=0;i<N0;i++)out[i]=(F->R2[i]-F->R2[N0+i])/(2.0*eps);
}

/* ---- torus (Neimark-Sacker) ---- */
/* The twisted collocation operator C(theta) of the variational equation, Y(s+1) = exp(i theta) Y(s): the rows are
   the collocation equations, the columns the nodal values of Y.  The coefficients are the collocation Jacobian's,
   with a phase exp(i theta e) where e counts how often a node or a delayed lookup wraps around the period.
   Fills the complex matrix (Cre,Cim) (N0 x N0, row major, zeroed by the caller) if they are not NULL, and/or
   applies the operator to Y = U + iV, giving (oRe,oIm) (zeroed by the caller).  Central differences. */
static void cplx_colloc(DDEPO *P,const double *X,double T,double theta,double *Cre,double *Cim,
			const double *U,const double *V,double *oRe,double *oIm)
{
  int j,i,k,a,b,u,m=P->ncol,n=P->n,nm=P->nmesh,N0=nm*n;
#define CADD(row,col,cr,ci) do{ \
    if(Cre){Cre[(size_t)(row)*N0+(col)]+=(cr);Cim[(size_t)(row)*N0+(col)]+=(ci);} \
    if(oRe){oRe[row]+=(cr)*U[col]-(ci)*V[col];oIm[row]+=(cr)*V[col]+(ci)*U[col];} }while(0)
  for(j=0;j<P->ntst;j++)
    for(i=0;i<m;i++){
      double s=((double)j+P->gz[i])/P->ntst;
      double xc[64],f[64],f1[64],f2[64],xp[64],A[64*64],Bm[DD_MAXUSED][64],tau[DD_MAXSLOT];
      int nused,slots,uvar[DD_MAXUSED],uslot[DD_MAXUSED],jd[DD_MAXUSED],wd[DD_MAXUSED];
      double uval[DD_MAXUSED],sgd[DD_MAXUSED];
      for(a=0;a<n;a++){
	double v=0.0;
	for(k=0;k<=m;k++)v+=P->L[i*(m+1)+k]*X[((j*m+k)%nm)*n+a];
	xc[a]=v;
      }
      H_record=1;H_nused=0;H_pvar=-1;H_pslot=-1;H_ph=0.0;
      call_rhs(P,X,s,T,xc,f);
      H_record=0;
      nused=H_nused;slots=H_nseen;
      for(u=0;u<nused;u++){uvar[u]=H_uvar[u];uslot[u]=H_uslot[u];uval[u]=H_uval[u];}
      for(u=0;u<slots;u++)tau[u]=H_tau[u];
      for(b=0;b<n;b++){
	double h=1e-5*(1.0+fabs(xc[b]));
	memcpy(xp,xc,sizeof(double)*n);
	xp[b]=xc[b]+h;call_rhs(P,X,s,T,xp,f1);
	xp[b]=xc[b]-h;call_rhs(P,X,s,T,xp,f2);
	for(a=0;a<n;a++)A[a*n+b]=(f1[a]-f2[a])/(2.0*h);
      }
      for(u=0;u<nused;u++){
	double h=1e-5*(1.0+fabs(uval[u])),sdr,fl,sg;
	int jj;
	H_pvar=uvar[u];H_pslot=uslot[u];H_ph=h;
	call_rhs(P,X,s,T,xc,f1);
	H_ph=-h;
	call_rhs(P,X,s,T,xc,f2);
	H_pvar=-1;H_ph=0.0;
	for(a=0;a<n;a++)Bm[u][a]=(f1[a]-f2[a])/(2.0*h);
	/* where the delayed value sits: s - tau/T = sd - w with sd in [0,1) and w whole periods back */
	sdr=s-tau[uslot[u]]/T;
	fl=floor(sdr);
	wd[u]=-(int)fl;
	sg=(sdr-fl)*P->ntst;
	jj=(int)floor(sg);
	if(jj>=P->ntst)jj=P->ntst-1;
	if(jj<0)jj=0;
	jd[u]=jj;sgd[u]=sg-jj;
      }
      for(a=0;a<n;a++){
	int row=(j*m+i)*n+a;
	for(k=0;k<=m;k++){
	  int g=(j*m+k)%nm,e=(j*m+k)/nm;                       /* e=1 for the last node of the last interval */
	  double lp=P->Lp[i*(m+1)+k]*P->ntst,cr=cos(theta*e),ci=sin(theta*e);
	  CADD(row,g*n+a,lp*cr,lp*ci);
	  for(b=0;b<n;b++){
	    double c=-T*A[a*n+b]*P->L[i*(m+1)+k];
	    CADD(row,g*n+b,c*cr,c*ci);
	  }
	}
	for(u=0;u<nused;u++){
	  int kd;
	  for(kd=0;kd<=m;kd++){
	    int g=(jd[u]*m+kd)%nm,e=(jd[u]*m+kd)/nm-wd[u];
	    double c=-T*Bm[u][a]*lag(m,kd,sgd[u]),cr=cos(theta*e),ci=sin(theta*e);
	    CADD(row,g*n+uvar[u],c*cr,c*ci);
	  }
	}
      }
    }
#undef CADD
}

/* E2 = C(theta) (U+iV): real part in out[0..N0), imaginary part in out[N0..2N0) */
static void tr_E2(DDEPO *P,FoldW *F,const double *z,double *out)
{
  int N0=F->N0,yo=F->yo;
  fold_setpars(P,F,z);
  memset(out,0,sizeof(double)*2*N0);
  cplx_colloc(P,z,z[N0],z[N0+3],NULL,NULL,z+yo,z+yo+N0,out,out+N0);
}

static void curve_E2(DDEPO *P,FoldW *F,const double *z,double *out)
{
  if(F->kind==0)fold_E2(P,F,z,out);else if(F->kind==1)pd_E2(P,F,z,out);else tr_E2(P,F,z,out);
}

static void curve_residual(DDEPO *P,FoldW *F,const double *z,double *R)
{
  int N0=F->N0,i,yo=F->yo;
  int ne2=(F->kind==0)?N0+1:((F->kind==1)?N0:2*N0);
  double s=0.0;
  fold_setpars(P,F,z);
  full_F(P,z,z[N0],R);
  curve_E2(P,F,z,R+N0+1);
  for(i=0;i<ne2;i++)s+=z[yo+i]*z[yo+i];
  R[N0+1+ne2]=s-1.0;
  if(F->kind==2){                      /* phase gauge: Im <Yref,Y> = 0 */
    s=0.0;
    for(i=0;i<N0;i++)s+=F->Yref[i]*z[yo+N0+i]-F->Yref[N0+i]*z[yo+i];
    R[N0+2+ne2]=s;
  }
}

/* K is (M-1) x M, row major */
static void curve_jacobian(DDEPO *P,FoldW *F,const double *z)
{
  int N0=F->N0,M=F->M,ld=N0+2,r,c,yo=F->yo;
  int ne2=(F->kind==0)?N0+1:((F->kind==1)?N0:2*N0);
  int nfd=N0+2+((F->kind==2)?1:0);          /* last column index of (X,T,p1,p2[,theta]) */
  double dp;
  double *K=F->K;
  memset(K,0,sizeof(double)*(size_t)(M-1)*M);
  fold_setpars(P,F,z);
  assemble(P,z,z[N0],z[N0+1],F->Mj,ld,F->Rt);
  for(r=0;r<=N0;r++){
    for(c=0;c<=N0;c++)K[(size_t)r*M+c]=F->Mj[(size_t)r*ld+c];               /* E1: X,T */
    K[(size_t)r*M+N0+1]=F->Mj[(size_t)r*ld+N0+1];                           /* E1: p1  */
  }
  if(F->kind==0){
    for(r=0;r<=N0;r++)for(c=0;c<=N0;c++)K[(size_t)(N0+1+r)*M+yo+c]=F->Mj[(size_t)r*ld+c];   /* E2: W,wT */
  }
  else if(F->kind==1){
    /* E2 against Y: L = (doubled Jacobian)(first block row), first column block minus the second */
    int ld2=2*N0+2,i;
    for(i=0;i<N0;i++){F->X2[i]=z[i];F->X2[N0+i]=z[i];}
    assemble(&F->P2,F->X2,2.0*z[N0],z[N0+1],F->Mj2,ld2,F->Rt2);
    for(r=0;r<N0;r++)for(c=0;c<N0;c++)
      K[(size_t)(N0+1+r)*M+yo+c]=F->Mj2[(size_t)r*ld2+c]-F->Mj2[(size_t)r*ld2+N0+c];
  }
  else{
    /* E2 against (U,V): the real form [[Cre,-Cim],[Cim,Cre]] of the twisted operator */
    double *Cre=(double *)calloc((size_t)N0*N0,sizeof(double)),*Cim=(double *)calloc((size_t)N0*N0,sizeof(double));
    if(Cre&&Cim){
      cplx_colloc(P,z,z[N0],z[N0+3],Cre,Cim,NULL,NULL,NULL,NULL);
      for(r=0;r<N0;r++)for(c=0;c<N0;c++){
	K[(size_t)(N0+1+r)*M+yo+c]=Cre[(size_t)r*N0+c];
	K[(size_t)(N0+1+r)*M+yo+N0+c]=-Cim[(size_t)r*N0+c];
	K[(size_t)(N0+1+N0+r)*M+yo+c]=Cim[(size_t)r*N0+c];
	K[(size_t)(N0+1+N0+r)*M+yo+N0+c]=Cre[(size_t)r*N0+c];
      }
    }
    free(Cre);free(Cim);
  }
  /* E1: p2 */
  dp=1e-6*(1.0+fabs(z[N0+2]));
  full_F(P,z,z[N0],F->F1);
  P->setpar(P->ctx,F->ip2,z[N0+2]+dp);
  full_F(P,z,z[N0],F->F2);
  P->setpar(P->ctx,F->ip2,z[N0+2]);
  for(r=0;r<=N0;r++)K[(size_t)r*M+N0+2]=(F->F2[r]-F->F1[r])/dp;
  /* E2 against X,T,p1,p2 (and theta) by finite differences of E2 */
  curve_E2(P,F,z,F->E2a);
  memcpy(F->zt,z,sizeof(double)*M);
  for(c=0;c<=nfd;c++){
    double h=1e-5*(1.0+fabs(z[c]));
    F->zt[c]=z[c]+h;
    curve_E2(P,F,F->zt,F->E2b);
    F->zt[c]=z[c];
    for(r=0;r<ne2;r++)K[(size_t)(N0+1+r)*M+c]=(F->E2b[r]-F->E2a[r])/h;
  }
  fold_setpars(P,F,z);
  /* E3: normalisation of the null vector (and the gauge of its phase for a torus) */
  for(c=0;c<ne2;c++)K[(size_t)(N0+1+ne2)*M+yo+c]=2.0*z[yo+c];
  if(F->kind==2)
    for(c=0;c<N0;c++){
      K[(size_t)(N0+2+ne2)*M+yo+c]=-F->Yref[N0+c];
      K[(size_t)(N0+2+ne2)*M+yo+N0+c]=F->Yref[c];
    }
}

/* solve [K; a^T] x = rhs (the last row is a) */
static int fold_solve(FoldW *F,const double *a,double *x)
{
  int M=F->M;
  memcpy(F->A,F->K,sizeof(double)*(size_t)(M-1)*M);
  memcpy(F->A+(size_t)(M-1)*M,a,sizeof(double)*M);
  if(lu_factor(M,F->A,F->piv))return 1;
  lu_solve_mat(M,F->A,F->piv,x,1);
  return 0;
}

/* Newton on E1,E2,E3 and the extra equation a.(z-zref)=0 */
static int fold_newton(DDEPO *P,FoldW *F,double *z,const double *a,const double *zref,int maxit,int *iters)
{
  int M=F->M,it,i;
  double *R=(double *)malloc(sizeof(double)*(M-1));
  double *x=(double *)malloc(sizeof(double)*M);
  int ok=0;
  if(!R||!x){free(R);free(x);return 0;}
  for(it=0;it<=maxit;it++){
    double nr=0.0,ex=0.0;
    curve_residual(P,F,z,R);
    for(i=0;i<M-1;i++){if(!(fabs(R[i])<1e30)){nr=1e30;break;}if(fabs(R[i])>nr)nr=fabs(R[i]);}
    for(i=0;i<M;i++)ex+=a[i]*(z[i]-zref[i]);
    if(fabs(ex)>nr)nr=fabs(ex);
    if(nr<5e-8){ok=1;break;}
    if(it==maxit||nr>1e20)break;
    curve_jacobian(P,F,z);
    for(i=0;i<M-1;i++)x[i]=-R[i];
    x[M-1]=-ex;
    if(fold_solve(F,a,x))break;
    {
      double mx=0.0;
      for(i=0;i<M;i++){if(!(fabs(x[i])<1e30)){mx=1e30;break;}if(fabs(x[i])>mx)mx=fabs(x[i]);}
      if(mx>1e5)break;
    }
    for(i=0;i<M;i++)z[i]+=x[i];
  }
  *iters=it;
  free(R);free(x);
  if(ok){fold_setpars(P,F,z);}
  return ok;
}

/* tangent of the curve at z: [K; c^T] t = e_M, normalised in the weighted norm */
static int fold_tangent(DDEPO *P,FoldW *F,const double *z,const double *c,double *t)
{
  int M=F->M,i;
  double s=0.0;
  curve_jacobian(P,F,z);
  for(i=0;i<M-1;i++)t[i]=0.0;
  t[M-1]=1.0;
  if(fold_solve(F,c,t))return 1;
  for(i=0;i<M;i++)s+=F->wt[i]*t[i]*t[i];
  if(!(s>0.0)||s!=s)return 1;
  s=sqrt(s);
  for(i=0;i<M;i++)t[i]/=s;
  return 0;
}

/* inverse iteration for the null vector of the leading nu x nu block of A (row major, ld), starting from w */
static int null_by_inverse_iteration(const double *A0,int ld,int nu,double *w)
{
  int i,j,it,ret=1;
  double *A=(double *)malloc(sizeof(double)*(size_t)nu*nu);
  int *piv=(int *)malloc(sizeof(int)*nu);
  double amax=0.0,mx;
  if(!A||!piv){free(A);free(piv);return 1;}
  for(i=0;i<nu;i++)for(j=0;j<nu;j++){
    A[(size_t)i*nu+j]=A0[(size_t)i*ld+j];
    if(fabs(A[(size_t)i*nu+j])>amax)amax=fabs(A[(size_t)i*nu+j]);
  }
  for(i=0;i<nu;i++)A[(size_t)i*nu+i]-=1e-9*amax;
  if(lu_factor(nu,A,piv))goto out;
  for(it=0;it<5;it++){
    double s=0.0;
    lu_solve_mat(nu,A,piv,w,1);
    for(i=0;i<nu;i++)s+=w[i]*w[i];
    if(!(s>0.0)||s!=s)goto out;
    s=sqrt(s);
    for(i=0;i<nu;i++)w[i]/=s;
  }
  mx=0.0;
  for(i=0;i<nu;i++)if(fabs(w[i])>mx)mx=fabs(w[i]);
  ret=(mx>0.0)?0:1;
out:
  free(A);free(piv);
  return ret;
}

/* the right null vector of the collocation Jacobian (X,T columns): the fold's (W,wT) */
static int fold_null_vector(DDEPO *P,FoldW *F,const double *X,double T,double p1,double p2,double *w)
{
  int N0=F->N0,Nu=N0+1,ld=N0+2,i;
  P->setpar(P->ctx,P->ip,p1);P->setpar(P->ctx,F->ip2,p2);
  assemble(P,X,T,p1,F->Mj,ld,F->Rt);
  for(i=0;i<Nu;i++)w[i]=1.0+0.2*sin(1.3*i);
  return null_by_inverse_iteration(F->Mj,ld,Nu,w);
}

/* the antiperiodic null vector Y of L at a period doubling; starts from the Floquet eigenfunction */
static int pd_null_vector(DDEPO *P,FoldW *F,const double *X,double T,double p1,double p2,double *Y)
{
  int N0=F->N0,ld2=2*N0+2,i,r,c;
  double *L=(double *)malloc(sizeof(double)*(size_t)N0*N0),mu=0.0,mm=0.0,s=0.0;
  int ret=1;
  if(!L)return 1;
  P->setpar(P->ctx,P->ip,p1);P->setpar(P->ctx,F->ip2,p2);
  for(i=0;i<N0;i++){F->X2[i]=X[i];F->X2[N0+i]=X[i];}
  assemble(&F->P2,F->X2,2.0*T,p1,F->Mj2,ld2,F->Rt2);
  for(r=0;r<N0;r++)for(c=0;c<N0;c++)L[(size_t)r*N0+c]=F->Mj2[(size_t)r*ld2+c]-F->Mj2[(size_t)r*ld2+N0+c];
  /* the start: the Floquet eigenfunction of the monodromy matrix (multiplier near -1) */
  {
    DDEPO Pt=*P;
    int ok;
    Pt.X=(double *)malloc(sizeof(double)*N0);
    if(!Pt.X){free(L);return 1;}
    memcpy(Pt.X,X,sizeof(double)*N0);Pt.T=T;
    ok=(ddepo_floquet_mode(&Pt,-1.0,Y,&mu,&mm)==0&&fabs(mu+1.0)<0.2);
    free(Pt.X);
    if(!ok){free(L);return 2;}
  }
  for(i=0;i<N0;i++)s+=Y[i]*Y[i];
  s=sqrt(s);
  for(i=0;i<N0;i++)Y[i]/=s;
  ret=null_by_inverse_iteration(L,N0,N0,Y);
  free(L);
  return ret;
}

/* start of a torus curve: the rotation angle from the Floquet multipliers, then the complex null vector of
   the twisted operator (inverse iteration on its real form) */
static int tr_start(DDEPO *P,FoldW *F,double p1,double p2,double *theta,double *Y)
{
  int N0=F->N0,N2=2*N0,i,r,c,k,nm,ret=1;
  double mr[70],mi[70],te=0.0,best=1e9,th=0.0;
  double *Cre,*Cim,*R2;
  P->setpar(P->ctx,P->ip,p1);P->setpar(P->ctx,F->ip2,p2);
  nm=ddepo_floquet(P,mr,mi,12,&te);
  for(k=0;k<nm;k++)
    if(fabs(mi[k])>1e-6){
      double d=fabs(hypot(mr[k],mi[k])-1.0);
      if(d<best){best=d;th=atan2(fabs(mi[k]),mr[k]);}
    }
  if(best>0.1)return 2;                                  /* no complex multiplier near the unit circle */
  Cre=(double *)calloc((size_t)N0*N0,sizeof(double));Cim=(double *)calloc((size_t)N0*N0,sizeof(double));
  R2=(double *)malloc(sizeof(double)*(size_t)N2*N2);
  if(!Cre||!Cim||!R2){free(Cre);free(Cim);free(R2);return 1;}
  cplx_colloc(P,P->X,P->T,th,Cre,Cim,NULL,NULL,NULL,NULL);
  for(r=0;r<N0;r++)for(c=0;c<N0;c++){
    R2[(size_t)r*N2+c]=Cre[(size_t)r*N0+c];
    R2[(size_t)r*N2+N0+c]=-Cim[(size_t)r*N0+c];
    R2[(size_t)(N0+r)*N2+c]=Cim[(size_t)r*N0+c];
    R2[(size_t)(N0+r)*N2+N0+c]=Cre[(size_t)r*N0+c];
  }
  for(i=0;i<N2;i++)Y[i]=1.0+0.2*sin(1.1*i);
  ret=null_by_inverse_iteration(R2,N2,N2,Y);
  free(Cre);free(Cim);free(R2);
  if(ret)return 1;
  memcpy(F->Yref,Y,sizeof(double)*N2);
  *theta=th;
  return 0;
}

/* Follow the curve from z (a point on it with a null vector), see ddepo_fold_curve */
static int curve_follow(DDEPO *P,FoldW *F,double *z,double pmin,double pmax,double qmin,double qmax,int dir,
			double ds,double dsmin,double dsmax,int nmax,ddepo_cb2 cb,void *rctx)
{
  int N0=F->N0,M=F->M,i,step,npts=0,iters=0;
  double *t=NULL,*zp=NULL,*zn=NULL,*a=NULL,*c=NULL,*zold=NULL;
  t=(double *)malloc(sizeof(double)*M);zp=(double *)malloc(sizeof(double)*M);zn=(double *)malloc(sizeof(double)*M);
  a=(double *)malloc(sizeof(double)*M);c=(double *)malloc(sizeof(double)*M);zold=(double *)malloc(sizeof(double)*M);
  if(!t||!zp||!zn||!a||!c||!zold){npts=-1;goto done;}
  /* a curve drawn in a finite box should not be a few long chords: limit the step to a fraction of its size */
  if(fabs(pmax-pmin)<1e20&&fabs(qmax-qmin)<1e20){
    double dsbox=0.04*hypot(pmax-pmin,qmax-qmin);
    if(dsmax>dsbox)dsmax=dsbox;
    if(ds>dsmax)ds=dsmax;
  }
  /* weights of the arclength norm: rms of X, T relative, the two parameters; none for the null vector */
  for(i=0;i<M;i++)F->wt[i]=0.0;
  for(i=0;i<N0;i++)F->wt[i]=1.0/N0;
  F->wt[N0]=1.0/(z[N0]*z[N0]);F->wt[N0+1]=1.0;F->wt[N0+2]=1.0;
  if(F->kind==2)F->wt[N0+3]=1.0;               /* the rotation angle */
  update_ref(P);
  /* correct the start with p2 held fixed: the exact point of the curve at this p2 */
  for(i=0;i<M;i++){a[i]=0.0;zold[i]=z[i];}
  a[N0+2]=1.0;
  if(!fold_newton(P,F,z,a,zold,12,&iters)){npts=-2;goto done;}
  memcpy(P->X,z,sizeof(double)*N0);P->T=z[N0];update_ref(P);
  if(F->kind==2){memcpy(F->Yref,z+F->yo,sizeof(double)*2*N0);P->theta=z[N0+3];}
  /* first tangent: a fixed pseudo-random border, then oriented by dir along p1 */
  for(i=0;i<M;i++)c[i]=sin(0.37*i+0.5)+0.1;
  if(fold_tangent(P,F,z,c,t)){npts=-3;goto done;}
  if(dir*t[N0+1]<0.0)for(i=0;i<M;i++)t[i]=-t[i];
  npts=1;
  if(cb&&cb(rctx,0,z[N0+1],z[N0+2],z[N0],z,0)){P->stop=4;goto done;}
  for(step=1;step<=nmax;step++){
    int ok=0,last=0,nland=0;
    double dsl=ds;
    memcpy(zold,z,sizeof(double)*M);
    if(P->poll&&P->poll(P->pollctx)){
      P->stop=5;
      if(cb)cb(rctx,step-1,z[N0+1],z[N0+2],z[N0],z,2);
      break;
    }
    while(!ok){
      int it2=0;
      for(i=0;i<M;i++){zp[i]=zold[i]+dsl*t[i];zn[i]=zp[i];a[i]=F->wt[i]*t[i];}
      if(!fold_newton(P,F,zn,a,zp,8,&it2)){
	dsl*=0.5;
	if(dsl<dsmin){P->stop=3;if(npts>0)goto finish;npts=-2;goto done;}
	last=0;nland=0;
	continue;
      }
      iters=it2;
      {
	double q1=zn[N0+1],q2=zn[N0+2],fr=1.0;
	int out=0;
	if(q1<pmin){out=1;fr=(pmin-zold[N0+1])/(q1-zold[N0+1]);}
	if(q1>pmax){out=1;fr=(pmax-zold[N0+1])/(q1-zold[N0+1]);}
	if(q2<qmin){double f2=(qmin-zold[N0+2])/(q2-zold[N0+2]);out=1;if(f2<fr)fr=f2;}
	if(q2>qmax){double f2=(qmax-zold[N0+2])/(q2-zold[N0+2]);out=1;if(f2<fr)fr=f2;}
	if(out){
	  /* leave the box: land on its boundary by shrinking the step */
	  double tolb=1e-6*(1.0+fabs(pmax-pmin)+fabs(qmax-qmin));
	  int near_b=(fabs(q1-pmin)<tolb||fabs(q1-pmax)<tolb||fabs(q2-qmin)<tolb||fabs(q2-qmax)<tolb);
	  if(near_b&&nland>0){last=1;ok=1;break;}
	  if(nland>=8||!(fr>0.0&&fr<1.0)){
	    P->stop=0;goto finish;       /* accept the last point inside: the previous one, finish */
	  }
	  dsl*=fr;nland++;
	  continue;
	}
      }
      ok=1;
    }
    memcpy(z,zn,sizeof(double)*M);
    memcpy(P->X,z,sizeof(double)*N0);P->T=z[N0];
    fold_setpars(P,F,z);
    update_ref(P);
    npts++;
    if(F->kind==2){
      memcpy(F->Yref,z+F->yo,sizeof(double)*2*N0);     /* the phase gauge follows the null vector */
      P->theta=z[N0+3];
    }
    if(step==nmax){last=1;P->stop=1;}
    /* new tangent, oriented along the old one */
    for(i=0;i<M;i++)c[i]=F->wt[i]*t[i];
    if(!last){
      double *tn=(double *)malloc(sizeof(double)*M);
      if(!tn||fold_tangent(P,F,z,c,tn)){free(tn);P->stop=3;if(cb)cb(rctx,step,z[N0+1],z[N0+2],z[N0],z,2);break;}
      memcpy(t,tn,sizeof(double)*M);free(tn);
    }
    if(cb&&cb(rctx,step,z[N0+1],z[N0+2],z[N0],z,last?2:0)){P->stop=4;break;}
    if(last)break;
    if(iters<=3){ds=dsl*1.5;if(ds>dsmax)ds=dsmax;}
    else if(iters>=6){ds=dsl*0.7;if(ds<dsmin)ds=dsmin;}
    else ds=dsl;
  }
  goto done;
finish:
  if(cb)cb(rctx,step-1,zold[N0+1],zold[N0+2],zold[N0],zold,2);
done:
  free(t);free(zp);free(zn);free(a);free(c);free(zold);
  return npts;
}

static int curve_start(DDEPO *P,int ip2,int kind,double pmin,double pmax,double qmin,double qmax,int dir,
		       double ds,double dsmin,double dsmax,int nmax,ddepo_cb2 cb,void *rctx)
{
  FoldW F;
  int N0=P->nmesh*P->n,npts,rc;
  double *z=NULL,p1,p2;
  P->stop=0;
  if(P->ntst*P->ncol*P->n>P->maxunk)return -1;
  if(!foldw_alloc(&F,P,N0,ip2,kind))return -1;
  z=(double *)malloc(sizeof(double)*F.M);
  if(!z){foldw_free(&F);return -1;}
  p1=P->getpar(P->ctx,P->ip);p2=P->getpar(P->ctx,ip2);
  memcpy(z,P->X,sizeof(double)*N0);z[N0]=P->T;z[N0+1]=p1;z[N0+2]=p2;
  update_ref(P);
  if(kind==0)rc=fold_null_vector(P,&F,P->X,P->T,p1,p2,z+F.yo);
  else if(kind==1)rc=pd_null_vector(P,&F,P->X,P->T,p1,p2,z+F.yo);
  else rc=tr_start(P,&F,p1,p2,z+N0+3,z+F.yo);
  if(rc){free(z);foldw_free(&F);return -4;}
  npts=curve_follow(P,&F,z,pmin,pmax,qmin,qmax,dir,ds,dsmin,dsmax,nmax,cb,rctx);
  free(z);
  foldw_free(&F);
  return npts;
}

int ddepo_fold_curve(DDEPO *P,int ip2,double pmin,double pmax,double qmin,double qmax,int dir,
		     double ds,double dsmin,double dsmax,int nmax,ddepo_cb2 cb,void *rctx)
{
  return curve_start(P,ip2,0,pmin,pmax,qmin,qmax,dir,ds,dsmin,dsmax,nmax,cb,rctx);
}

int ddepo_pd_curve(DDEPO *P,int ip2,double pmin,double pmax,double qmin,double qmax,int dir,
		   double ds,double dsmin,double dsmax,int nmax,ddepo_cb2 cb,void *rctx)
{
  return curve_start(P,ip2,1,pmin,pmax,qmin,qmax,dir,ds,dsmin,dsmax,nmax,cb,rctx);
}

int ddepo_tr_curve(DDEPO *P,int ip2,double pmin,double pmax,double qmin,double qmax,int dir,
		   double ds,double dsmin,double dsmax,int nmax,ddepo_cb2 cb,void *rctx)
{
  return curve_start(P,ip2,2,pmin,pmax,qmin,qmax,dir,ds,dsmin,dsmax,nmax,cb,rctx);
}
