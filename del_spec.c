/* Characteristic roots of a linear(ized) delay equation

        y'(t) = A_0 y(t) + sum_{k=1}^m A_k y(t - tau_k)

   Characteristic equation:  det Delta(lambda)=0,
        Delta(lambda) = lambda I - A_0 - sum_k A_k exp(-lambda tau_k).

   Method (Breda-Maset-Vermiglio):
     1. Discretize the infinitesimal generator of the solution semigroup by
        Chebyshev collocation on [-tau_max,0].  The eigenvalues of the
        resulting (N+1)n x (N+1)n matrix approximate the rightmost roots
        with spectral accuracy (they are poor for roots far to the left,
        which we do not care about).
     2. Each of the rightmost approximations is then polished by Newton's
        method on the bordered system
              Delta(lambda) v = 0 ,   w^H v = 1
        so the final roots are accurate to machine precision regardless of
        N.
   Coefficient layout is that of do_delay_sing(): coef[i+n*j+n*n*k] is the
   (i,j) entry of A_k stored so that A_k[j][i]=coef[j*n+i+n*n*k]  (row=
   equation j, column=variable i).                                       */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "del_stab.h"

void eigen(int n,double *a,double *ev,double *work,int *ierr);

typedef double complex cplx;

/* A_k[r][c] */
#define AK(k,r,c) (coef[(c)+(r)*n+(size_t)n*n*(k)])

/* ----- small dense complex linear solver (Gauss, partial pivoting) ----- */
static int csolve(cplx *a,cplx *b,int n)   /* a is n x n row-major, b rhs */
{
  int i,j,k,p;
  cplx t;
  for(k=0;k<n;k++){
    p=k;
    for(i=k+1;i<n;i++)if(cabs(a[i*n+k])>cabs(a[p*n+k]))p=i;
    if(cabs(a[p*n+k])<1e-300)return 1;
    if(p!=k){
      for(j=0;j<n;j++){t=a[k*n+j];a[k*n+j]=a[p*n+j];a[p*n+j]=t;}
      t=b[k];b[k]=b[p];b[p]=t;
    }
    for(i=k+1;i<n;i++){
      t=a[i*n+k]/a[k*n+k];
      if(t==0)continue;
      for(j=k;j<n;j++)a[i*n+j]-=t*a[k*n+j];
      b[i]-=t*b[k];
    }
  }
  for(i=n-1;i>=0;i--){
    t=b[i];
    for(j=i+1;j<n;j++)t-=a[i*n+j]*b[j];
    b[i]=t/a[i*n+i];
  }
  return 0;
}

static void build_delta(const double *coef,const double *delay,int n,int m,
                        cplx lam,cplx *D,cplx *Dp)
{
  int i,j,k;
  cplx e;
  for(i=0;i<n;i++)for(j=0;j<n;j++){
    D[i*n+j]=(i==j?lam:0.0)-AK(0,i,j);
    Dp[i*n+j]=(i==j?1.0:0.0);
  }
  for(k=1;k<=m;k++){
    e=cexp(-lam*delay[k-1]);
    for(i=0;i<n;i++)for(j=0;j<n;j++){
      D[i*n+j]-=AK(k,i,j)*e;
      Dp[i*n+j]+=delay[k-1]*AK(k,i,j)*e;
    }
  }
}

/* Newton polish of a root; returns 1 on convergence */
static int refine_root(const double *coef,const double *delay,int n,int m,
                       cplx *lam)
{
  int i,j,it,nn=n+1;
  cplx *D=malloc(sizeof(cplx)*n*n),*Dp=malloc(sizeof(cplx)*n*n);
  cplx *J=malloc(sizeof(cplx)*nn*nn),*r=malloc(sizeof(cplx)*nn);
  cplx *v=malloc(sizeof(cplx)*n),*w=malloc(sizeof(cplx)*n);
  cplx l=*lam,s;
  double nv;
  int ok=0;
  if(!D||!Dp||!J||!r||!v||!w)goto done;

  /* starting null vector by inverse iteration at a slightly shifted lambda */
  build_delta(coef,delay,n,m,l+1e-7*(1.0+cabs(l)),D,Dp);
  for(i=0;i<n;i++)v[i]=1.0/(i+1.0);
  for(it=0;it<3;it++){
    memcpy(J,D,sizeof(cplx)*n*n);
    for(i=0;i<n;i++)r[i]=v[i];
    if(csolve(J,r,n))break;
    nv=0;for(i=0;i<n;i++)nv+=creal(r[i]*conj(r[i]));
    nv=sqrt(nv);
    if(nv==0)break;
    for(i=0;i<n;i++)v[i]=r[i]/nv;
  }
  for(i=0;i<n;i++)w[i]=conj(v[i]);       /* normalization w^H v=1 */
  s=0;for(i=0;i<n;i++)s+=w[i]*v[i];
  for(i=0;i<n;i++)v[i]/=s;

  for(it=0;it<40;it++){
    build_delta(coef,delay,n,m,l,D,Dp);
    for(i=0;i<n;i++){
      s=0;for(j=0;j<n;j++)s+=D[i*n+j]*v[j];
      r[i]=-s;
    }
    s=0;for(j=0;j<n;j++)s+=w[j]*v[j];
    r[n]=1.0-s;
    for(i=0;i<n;i++){
      for(j=0;j<n;j++)J[i*nn+j]=D[i*n+j];
      s=0;for(j=0;j<n;j++)s+=Dp[i*n+j]*v[j];
      J[i*nn+n]=s;
    }
    for(j=0;j<n;j++)J[n*nn+j]=w[j];
    J[n*nn+n]=0.0;
    if(csolve(J,r,nn))goto done;
    for(i=0;i<n;i++)v[i]+=r[i];
    l+=r[n];
    if(cabs(r[n])<1e-13*(1.0+cabs(l))){ok=1;break;}
    if(!isfinite(creal(l))||!isfinite(cimag(l)))goto done;
  }
 done:
  if(ok)*lam=l;
  free(D);free(Dp);free(J);free(r);free(v);free(w);
  return ok;
}

/* order by decreasing real part */
static int cmp_re(const void *a,const void *b)
{
  double x=((const double *)a)[0],y=((const double *)b)[0];
  double xi=((const double *)a)[1],yi=((const double *)b)[1];
  if(x>y)return -1;
  if(x<y)return 1;
  if(xi>yi)return -1;
  if(xi<yi)return 1;
  return 0;
}

/* Returns the number of roots written to ev (<=nev), as pairs (re,im),
   sorted by decreasing real part.  *nunstab is the number of roots with
   Re>tol among all roots found.  *nrefined is the number that passed
   Newton polishing (the rest are collocation approximations).
   Returns -1 on failure.                                                */
int delay_spectrum(double *coef,double *delay,int n,int m,int nev,
                   double *ev,int *nunstab,int *nrefined,int nmax)
{
  int N,dim,i,j,k,l,ierr,ncand,nout=0,nref=0,nun=0;
  double tmax=0.0,*M=NULL,*work=NULL,*cand=NULL,*x=NULL,*bw=NULL,*D=NULL;
  double *lag=NULL;
  const double tol=1e-9;

  for(k=0;k<m;k++)if(delay[k]>tmax)tmax=delay[k];
  if(m==0||tmax<=0.0){
    /* ordinary ODE: eigenvalues of A_0 (eigen wants column-major) */
    N=0;
  } else {
    N=nmax/n-1;
    if(N>40)N=40;
    if(N<10)N=10;
  }
  dim=n*(N+1);
  M=calloc((size_t)dim*dim,sizeof(double));
  work=malloc(sizeof(double)*(dim+1));
  cand=malloc(sizeof(double)*2*(dim+1));
  if(!M||!work||!cand){free(M);free(work);free(cand);return -1;}

  if(N==0){
    for(i=0;i<n;i++)for(j=0;j<n;j++)M[i+j*n]=AK(0,i,j);
  } else {
    /* Chebyshev points x_j=cos(j pi/N) on [-1,1]; theta_j=tmax/2 (x_j-1) */
    x=malloc(sizeof(double)*(N+1));
    bw=malloc(sizeof(double)*(N+1));
    D=malloc(sizeof(double)*(N+1)*(N+1));
    lag=malloc(sizeof(double)*(N+1));
    for(j=0;j<=N;j++){
      x[j]=cos(j*3.14159265358979323846/N);
      bw[j]=((j%2)?-1.0:1.0)*((j==0||j==N)?0.5:1.0);
    }
    /* differentiation matrix, D[i][j] */
    for(i=0;i<=N;i++){
      double rowsum=0;
      for(j=0;j<=N;j++){
        if(i!=j){
          D[i*(N+1)+j]=(bw[j]/bw[i])/(x[i]-x[j]);
          rowsum+=D[i*(N+1)+j];
        }
      }
      D[i*(N+1)+i]=-rowsum;
    }
    for(i=0;i<=N;i++)for(j=0;j<=N;j++)D[i*(N+1)+j]*=2.0/tmax;
    /* generator matrix, stored column-major (M[r+c*dim]) :
       block row 0: A0 at node 0 + interpolated delayed terms
       block row i>=1: derivative rows                                    */
    for(i=1;i<=N;i++)
      for(j=0;j<=N;j++)
        for(l=0;l<n;l++)
          M[(i*n+l)+(size_t)(j*n+l)*dim]=D[i*(N+1)+j];
    for(l=0;l<n;l++)for(j=0;j<n;j++)
      M[l+(size_t)j*dim]+=AK(0,l,j);
    for(k=1;k<=m;k++){
      /* Lagrange basis values at theta=-tau_k  <->  x*=1-2 tau_k/tmax */
      double xs=1.0-2.0*delay[k-1]/tmax,den=0.0;
      int hit=-1;
      for(j=0;j<=N;j++){
        if(fabs(xs-x[j])<1e-14){hit=j;break;}
        lag[j]=bw[j]/(xs-x[j]);
        den+=lag[j];
      }
      if(hit>=0){for(j=0;j<=N;j++)lag[j]=(j==hit)?1.0:0.0;}
      else for(j=0;j<=N;j++)lag[j]/=den;
      for(j=0;j<=N;j++)
        for(l=0;l<n;l++)for(i=0;i<n;i++)
          M[l+(size_t)(j*n+i)*dim]+=AK(k,l,i)*lag[j];
    }
  }
  eigen(dim,M,cand,work,&ierr);
  if(ierr!=0){
    free(M);free(work);free(cand);free(x);free(bw);free(D);free(lag);
    return -1;
  }
  qsort(cand,dim,2*sizeof(double),cmp_re);

  /* polish the rightmost ones.  Take everything with Re>-margin so the
     unstable count is not truncated, and at least nev+6 candidates.     */
  ncand=0;
  for(i=0;i<dim;i++)if(cand[2*i]>tol)ncand++;
  ncand=ncand+nev+6;
  if(ncand>dim)ncand=dim;
  {
    double *out=malloc(sizeof(double)*2*ncand);
    int nout2=0;
    if(!out){free(M);free(work);free(cand);free(x);free(bw);free(D);free(lag);return -1;}
    for(i=0;i<ncand;i++){
      cplx l0=cand[2*i]+I*cand[2*i+1];
      int dup=0;
      if(N>0){
        if(refine_root(coef,delay,n,m,&l0))nref++;
        else if(cand[2*i]<-0.5*(1.0+cabs(l0))&&i>=nev) continue; /* spurious tail */
      }
      for(j=0;j<nout2;j++)
        if(fabs(out[2*j]-creal(l0))<1e-7*(1+cabs(l0))&&
           fabs(out[2*j+1]-cimag(l0))<1e-7*(1+cabs(l0)))dup=1;
      if(dup)continue;
      out[2*nout2]=creal(l0);out[2*nout2+1]=cimag(l0);nout2++;
    }
    qsort(out,nout2,2*sizeof(double),cmp_re);
    for(i=0;i<nout2;i++)if(out[2*i]>tol)nun++;
    nout=nout2<nev?nout2:nev;
    for(i=0;i<nout;i++){ev[2*i]=out[2*i];ev[2*i+1]=out[2*i+1];}
    free(out);
  }
  *nunstab=nun;
  *nrefined=nref;
  free(M);free(work);free(cand);free(x);free(bw);free(D);free(lag);
  return nout;
}
