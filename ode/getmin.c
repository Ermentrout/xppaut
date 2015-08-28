



#include <math.h>
fun(double *in,double *out,int nin,int nout,double *var,double *con)
{
  int nmin=(int)in[0],i;
 int nmax=(int)in[1];
 double d=var[nmin+1],dnew;
 for(i=nmin;i<=nmax;i++){
   dnew=fabs(var[i+1]);
   if(dnew<d){
     d=dnew;
   }
 }
 out[0]=d;
}
