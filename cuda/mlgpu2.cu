/*

Linux:
nvcc -lcublas -o MLGPU.SO --shared -Xcompiler -fpic mlgpu.cu 
  export PATH=/usr/local/cuda-7.0/bin:$PATH
 export LD_LIBRARY_PATH=/usr/local/cuda-7.0/lib64:$LD_LIBRARY_PATH

Mac:
nvcc -lcublas -o MLGPU.SO --shared -Xcompiler -fPIC  -Xcompiler -O3 mlgpu.cu

export PATH=/Developer/NVIDIA/CUDA-7.0/bin:$PATH
export DYLD_LIBRARY_PATH=/Developer/NVIDIA/CUDA-7.0/lib:$DYLD_LIBRARY_PATH


*/

#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include <curand.h> 
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#define MAXPAR 50
#define NTHREAD 128

#define real double

__constant__ real devp[MAXPAR];



real *devsum,*devy,*devyp,*devwgt;
/* needed for double to float */
real *hosty,*hostyp,*hostwgt;

int allocflag=0;
int recopywgtflag=0;

cublasHandle_t handle;
/* ML parameters   */ 

#define iapp devp[0]
#define phi devp[1]
#define va devp[2]
#define vb devp[3]
#define vc devp[4]
#define vd devp[5]
#define gca devp[6]
#define vk devp[7]
#define vl devp[8]
#define gk devp[9]
#define gl devp[10]
#define tsyn devp[11]
#define gsyn devp[12]
#define vsyn devp[13]
#define vth devp[14]
#define vshp devp[15]

/*  ML Functions */

__device__ real minf(real v)
{
  return .5*(1+tanh((v-va)/vb));
}
__device__ real ninf(real v)
{

  return .5*(1+tanh((v-vc)/vd));
}
__device__ real lamn(real  v)
{
  return phi*cosh((v-vc)/(2*vd));
}
__device__ real s_inf(real v)
{
  
  return 1.0/(1.0+exp(-(v-vth)/vshp)); 

}

/*  this is the right-hand side evaluation */



__global__ void update_rhs(real *y,real *yp,real *stot, int n)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;
    real v,w,s;
 /*   if(i<10){
    printf(" %d %g %g %g \n",i,y[i],y[i+n],y[i+2*n]);
    }
  */
    if(i<n){
     v=y[i];
     w=y[i+n];
     s=y[i+2*n];
     yp[i]=iapp-gl*(v-vl)-gca*(v-1.0)*minf(v)-gk*w*(v-vk)-gsyn*stot[i]*(v-vsyn);
     yp[i+n]=lamn(v)*(ninf(v)-w);
     yp[i+2*n]=(sinf(v)-s)/tsyn;
    }
}


/*  I will assume for now that the weight do not change throughout the simulation
    so they will only be loaded into the device once
*/
void allocate_ram(int n, double *w)
{
  int n3=3*n,nn=n*n;
  int i;
  if(allocflag==1)return;  /* already allocated */
  hostwgt=(real *)malloc(nn*sizeof(real));
  cudaMalloc((void**)&devwgt,nn*sizeof(real));
  /*  copy one time only for now */
  for(i=0;i<nn;i++)
      hostwgt[i]=(real)w[i];
  cudaMemcpy(devwgt,hostwgt,nn*sizeof(real),cudaMemcpyHostToDevice); 
  hosty=(real *)malloc(n3*sizeof(real));
  cudaMalloc((void**)&devy,n3*sizeof(real));
  hostyp=(real *)malloc(n3*sizeof(real));
  cudaMalloc((void**)&devyp,n3*sizeof(real));
    
  cudaMalloc((void**)&devsum,n*sizeof(real));
 cublasCreate(&handle);
  allocflag=1;  

}

/* do this every time */
void copy_par_to_dev(double *par)
{
  real dp[MAXPAR];
  int i;
  for(i=0;i<MAXPAR;i++)
    dp[i]=(real)par[i];
  cudaMemcpyToSymbol(devp,dp,sizeof(dp));
}

/* do this every time */

void copy_to_dev(int n,double  *w,double *y)
{
  int i;
  int n3=3*n,nn=n*n;

  for(i=0;i<n3;i++)
    hosty[i]=y[i];
  cudaMemcpy(devy,hosty,n3*sizeof(real),cudaMemcpyHostToDevice);
  if(recopywgtflag){

    for(i=0;i<nn;i++)
      hostwgt[i]=w[i];
  cudaMemcpy(devwgt,hostwgt,nn*sizeof(real),cudaMemcpyHostToDevice);
  }
}  

void copy_from_dev(double *ydot, int nn)
{
  int i;
  cudaMemcpy(hostyp,devyp,nn*sizeof(real),cudaMemcpyDeviceToHost);
  for(i=0;i<nn;i++)
    ydot[i]=(double)hostyp[i];
}






void update_sums(int n)
{
  real alpha=1.0,beta=0.0;


  /* single  */
/*  cublasSgemv(handle,CUBLAS_OP_T,n,n,&alpha,devwgt,n,&devy[2*n],1,&beta,devsum,1); */
  /* double */
 cublasDgemv(handle,CUBLAS_OP_T,n,n,&alpha,devwgt,n,&devy[2*n],1,&beta,devsum,1); 

 
}


extern "C"

void MLGPU(int nn,int ivar, double *par,double *var,double *z[50],double *ydot)
{
  int n=nn/3;
  double *y;
  int nblock=(int)((n+NTHREAD-1)/NTHREAD);
  y=var+ivar;
  allocate_ram( n, z[0]);
  copy_par_to_dev(par);
  copy_to_dev(n,z[0],y);
  update_sums(n);
  update_rhs<<<nblock,128>>>(devy,devyp,devsum,n);
  copy_from_dev(ydot,nn);
 
}


