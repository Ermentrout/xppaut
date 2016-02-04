/* equation free modeling 
example is due to Carlo Laing */
#include <math.h>
#include <sys/time.h>
#define MAXN 500

double x[MAXN],s[MAXN],i0[MAXN];

/* since the normal distribution isnt part of C... have to make a norma
   numbee generator */
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

/* this runs the simulation for a while and averages the slope over a number
   of trials 
*/
double onerun(double dt, int n,int ntrials, int niter, int nstart,double ibar,double isig,double nsig,double S,double tau)
{
  int ni,nt,i;
  double cs[MAXN];
  double ssum=0.0;
  double sbar;
  double sdt=1./sqrt(dt);
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
	x[i]+=dt*(1-cs[i]+(1+cs[i])*(i0[i]+sbar+sdt*norm()*nsig));
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
	
efmx(double *in,double *out,int nin,int nout,double *var,double *con)
{ 
  double S=in[0],ibar=in[1],isig=in[2],nsig=in[3],tau=in[4];
  double dt=.05;
  int ntrials=(int)in[6],niter=(int)in[5],nstart=(int)(.25*niter),n=(int)in[7];
  unsigned long start;
  unsigned long stop;
  
  struct timeval timer;
  gettimeofday(&timer,NULL);
  start=timer.tv_sec * 1000 + timer.tv_usec / 1000 ;
  /*  printf("nin=%d nout=%d\n",nin,nout);
      printf(" ntrial=%d  niter=%d nstart=%d ibar=%g isig=%g",ntrials,niter,nstart, ibar,isig); */

  out[0]=onerun(dt, n,ntrials, niter,nstart,ibar,isig,nsig,S,tau);
    gettimeofday(&timer,NULL);
   stop=timer.tv_sec * 1000 + timer.tv_usec / 1000 ;
   out[1]=(double)(stop-start);
   
}
