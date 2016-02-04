/* equation free modeling 
example is due to Carlo Laing */
#include <math.h>
#define MAXN 500

double x[MAXN],s[MAXN],i0[MAXN];


double BoxMuller;
int BoxMullerFlag;
double drand48();
#define Tpi 6.283185307

double norm()
{

  double fac,r,v1,v2;
  if(BoxMullerFlag==0){ 
    do {
      v1=2.0*drand48()-1.0;
      v2=2.0*drand48()-1.0;
      r=v1*v1+v2*v2;
    } while(r>=1.0);
    fac=sqrt(-2.0*log(r)/r);
    BoxMuller=v1*fac;
    BoxMullerFlag=1;
    return(v2*fac);
  }
  else {
    BoxMullerFlag=0;
    return(BoxMuller);
  }
}

double onerun(double dt, int n,int ntrials, int niter, int nstart,double ibar,double isig,double S,double tau)
{
  int ni,nt,i;
  double cs[MAXN];
  double ssum=0.0;
  double sbar;
  double sstart,f=(double)(niter-nstart)*dt;
  int n1=nstart+1;
  for(nt=0;nt<ntrials;nt++){
    for(i=0;i<n;i++){ /* initialize */
      x[i]=Tpi*drand48();
      i0[i]=ibar+isig*norm();
      s[i]=S;
    }

    for(ni=0;ni<niter;ni++){ /* integrate x,s for niter times */
      sbar=0;
      for(i=0;i<n;i++) /* get synaptic drive */
	sbar+=s[i];
      sbar=sbar/(double) n;
      if(ni==nstart)sstart=sbar;
      for(i=0;i<n;i++){
	cs[i]=cos(x[i]);
	s[i]+=dt*(.01*pow((1-cs[i]),10.0)-s[i])/tau;
	x[i]+=dt*(1-cs[i]+(1+cs[i])*(i0[i]+sbar));
      } /* all have been updated */
    }
    /* approx slope */
    for(i=0;i<n;i++) 
	sbar+=s[i];
    sbar=sbar/(double) n;
    ssum+=((sbar-sstart)/f); 
  }
  return (ssum/(double)ntrials);
}
	
main()
{ 
  double dt=.05,ibar=-.5,isig=.3333,S,tau=10.0;
  int ntrials=100,niter=200,nstart=50,n=100;
  int i;
  double sout;
  for(i=0;i<51;i++){
    S=-1.0+.08*i;
    sout=onerun(dt, n,ntrials, niter, nstart,ibar,isig, S,tau);
    printf("%g %g \n",S,sout);
  }
}
