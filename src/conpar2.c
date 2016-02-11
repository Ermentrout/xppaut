#include "auto_f2c.h"
#include "auto_c.h"
#include "auto_types.h"

#ifdef TIME
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
static double time_start (void) {
  struct rusage time;
  double seconds,microseconds;
  getrusage(RUSAGE_SELF,&time);
  seconds = (double)time.ru_utime.tv_sec;
  microseconds = (double)time.ru_utime.tv_usec;
  return seconds + microseconds/1e6;
}
static double time_end(double start) {
  struct rusage time;
  double seconds,microseconds;
  getrusage(RUSAGE_SELF,&time);
  seconds = (double)time.ru_utime.tv_sec;
  microseconds = (double)time.ru_utime.tv_usec;
  return (seconds + microseconds/1e6)-start;
}
#endif


/*This is the process function.  It is meant to be called either
  on a SMP using shared memory, or wrapped inside another
  routine for message passing*/

void *conpar_process(void * arg)
{
  integer icf_dim1, irf_dim1, d_dim1;
  integer a_dim1, a_dim2, b_dim1, b_dim2, c_dim1, c_dim2;
  
  /* Local variables */
  integer ipiv, jpiv, itmp;
  doublereal tpiv;
  integer i, l, k1, k2, m1, m2, ic, ir;
  doublereal rm;
  integer ir1, irp;
  doublereal piv;
  integer icp1;

  integer *nov, *nra, *nca;
  doublereal *a;
  integer *ncb;
  doublereal *b;
  integer *nbc, *nrc;
  doublereal *c, *d;
  integer *irf, *icf;
  integer loop_start,loop_end;



  nov = ((conpar_parallel_arglist *)arg)->nov;
  nra = ((conpar_parallel_arglist *)arg)->nra;
  nca = ((conpar_parallel_arglist *)arg)->nca;
  a = ((conpar_parallel_arglist *)arg)->a;
  ncb = ((conpar_parallel_arglist *)arg)->ncb;
  b = ((conpar_parallel_arglist *)arg)->b;
  nbc = ((conpar_parallel_arglist *)arg)->nbc;
  nrc = ((conpar_parallel_arglist *)arg)->nrc;
  c = ((conpar_parallel_arglist *)arg)->c;
  d = ((conpar_parallel_arglist *)arg)->d;
  irf = ((conpar_parallel_arglist *)arg)->irf;
  icf = ((conpar_parallel_arglist *)arg)->icf;
  loop_start = ((conpar_parallel_arglist *)arg)->loop_start;
  loop_end = ((conpar_parallel_arglist *)arg)->loop_end;

  /* In the default case we don't need to do anything special */
  if(global_conpar_type == CONPAR_DEFAULT) {
    ;
  }
  /* In the message passing case we set d to be
     0.0, do a sum here, and then do the final
     sum (with the true copy of d) in the
     master */
  else if (global_conpar_type == CONPAR_MPI) {
    for(i=0;i<(*ncb)*(*nrc);i++)
      d[i]=0.0;
  }
  /* In the shared memory case we create a local
     variable for doing this threads part of the
     sum, then we do a final sum into shared memory
     at the end */
  else if (global_conpar_type == CONPAR_PTHREADS) {

    ;

  }

  /* Note that the summation of the adjacent overlapped part of C */
  /* is delayed until REDUCE, in order to merge it with other communications.*/
  /* NA is the local NTST. */
  
  irf_dim1 = *nra;
  icf_dim1 = *nca;
  d_dim1 = *ncb;

  a_dim1 = *nca;
  a_dim2 = *nra;
  b_dim1 = *ncb;
  b_dim2 = *nra;
  c_dim1 = *nca;
  c_dim2 = *nrc;
  
  /* Condensation of parameters (Elimination of local variables). */
  m1 = *nov + 1;
  m2 = *nca - *nov;

  for (i = loop_start;i < loop_end; i++) {
    for (ic = m1; ic <= m2; ++ic) {
      ir1 = ic - *nov + 1;
      irp = ir1 - 1;
      icp1 = ic + 1;
      /*	     **Search for pivot (Complete pivoting) */
      piv = 0.0;
      ipiv = irp;
      jpiv = ic;
      for (k1 = irp; k1 <= *nra; ++k1) {
	int irf_k1_i = irf[-1 + k1 + i*irf_dim1];
	for (k2 = ic; k2 <= m2; ++k2) {
	  int icf_k2_i = icf[-1 + k2 + i*icf_dim1];
	  tpiv = a[-1 + icf_k2_i + a_dim1*(-1 + irf_k1_i + a_dim2*i)];
	  if (tpiv < 0.0) {
	    tpiv = -tpiv;
	  }
	  if (piv < tpiv) {
	    piv = tpiv;
	    ipiv = k1;
	    jpiv = k2;
	  }
	}
      }
      /*	     **Move indices */
      itmp = icf[-1 + ic + i*icf_dim1];
      icf[-1 + ic + i*icf_dim1] = icf[-1 + jpiv + i*icf_dim1];
      icf[-1 + jpiv + i*icf_dim1] = itmp;
      itmp = irf[-1 + irp + i*irf_dim1];
      irf[-1 + irp + i*irf_dim1] = irf[-1 + ipiv + i*irf_dim1];
      irf[-1 + ipiv + i*irf_dim1] = itmp;
      {
	int icf_ic_i = icf[-1 + ic + i*icf_dim1];
	int irf_irp_i = irf[-1 + irp + i*irf_dim1];
	int a_offset2 = a_dim1*(-1 + irf_irp_i + a_dim2*i);
	int b_offset2 = b_dim1*(-1 + irf_irp_i + b_dim2*i);
	/*	     **End of pivoting; elimination starts here */
	for (ir = ir1; ir <= *nra; ++ir) {
	  int irf_ir_i = irf[-1 + ir + i*irf_dim1];
	  int a_offset1 = a_dim1*(-1 + irf_ir_i + a_dim2*i);
	  int b_offset1 = b_dim1*(-1 + irf_ir_i + b_dim2*i);
	  rm = a[-1 + icf_ic_i + a_dim1*(-1 + irf_ir_i + a_dim2*i)]/a[-1 + icf_ic_i + a_dim1*(-1 + irf_irp_i + a_dim2*i)];
	  a[-1 + icf_ic_i + a_dim1*(-1 + irf_ir_i + a_dim2*i)] = rm;
	  if (rm != (double)0.) {
	    for (l = 0; l < *nov; ++l) {
	      a[l + a_offset1] -= rm * a[l + a_offset2];
	    }
	    for (l = icp1 -1; l < *nca; ++l) {
	      int icf_l_i = icf[l + i*icf_dim1];
	      a[-1 + icf_l_i + a_offset1] -= rm * a[-1 + icf_l_i + a_offset2];
	    }
	    for (l = 0; l < *ncb; ++l) {
	      b[l + b_offset1] -= rm * b[l + b_offset2];
	    }
	  }
	}
	for (ir = *nbc + 1; ir <= *nrc; ++ir) {
	  int c_offset1 = c_dim1*(-1 + ir + c_dim2*i);
	  int d_offset1 = (-1 + ir)*d_dim1;
	  rm = c[-1 + icf_ic_i + c_dim1*(-1 + ir + c_dim2*i)]/a[-1 + icf_ic_i + a_dim1*(-1 + irf_irp_i + a_dim2*i)];
	  c[-1 + icf_ic_i + c_dim1*(-1 + ir + c_dim2*i)]=rm;
	  if (rm != (double)0.) {
	    for (l = 0; l < *nov; ++l) {
	      c[l + c_offset1] -= rm * a[l + a_offset2];
	    }
	    for (l = icp1 -1 ; l < *nca; ++l) {
	      int icf_l_i = icf[l + i*icf_dim1];
	      c[-1 + icf_l_i + c_offset1] -= rm * a[-1 + icf_l_i + a_offset2];
	    }
	    for (l = 0; l < *ncb; ++l) {
	      /* 
		 A little explanation of what is going on here
		 is in order I believe.  This array is
		 created by a summation across all workers,
		 hence it needs a mutex to avoid concurrent
		 writes (in the shared memory case) or a summation
		 in the master (in the message passing case).
		 Since mutex's can be somewhat slow, we will do the
		 summation into a local variable, and then do a
		 final summation back into global memory when the
		 main loop is done.
	      */
	      /* Nothing special for the default case */
	      if(global_conpar_type == CONPAR_DEFAULT) {
		d[l + d_offset1] -= rm * b[l + b_offset2];
	      }
	      /* In the message passing case we sum into d,
		 which is a local variable initialized to 0.0.
		 We then sum our part with the masters part
		 in the master. */
	      else if (global_conpar_type == CONPAR_MPI) {
		d[l + d_offset1] -= rm * b[l + b_offset2];
	      }
	      /* In the shared memory case we sum into a local
		 variable our contribution, and then sum
		 into shared memory at the end (inside a mutex */
	      else if (global_conpar_type == CONPAR_PTHREADS) {

		;

	      }
	    }
	  }
	}
      }
    }
  }
  return NULL;
}



int 
conpar_default_wrapper(integer *nov, integer *na, integer *nra, integer *nca, doublereal *a, integer *ncb, doublereal *b, integer *nbc, integer *nrc, doublereal *c, doublereal *d, integer *irf, integer *icf)

{
    conpar_parallel_arglist data;
    data.nov = nov;
    data.nra = nra;
    data.nca = nca;
    data.a = a;
    data.ncb = ncb;
    data.b = b;
    data.nbc = nbc;
    data.nrc = nrc;
    data.c = c;
    data.d = d;
    data.irf = irf;
    data.icf = icf;
    data.loop_start = 0;
    data.loop_end = *na;
    conpar_process(&data);
    return 0;
}


int 
conpar(integer *nov, integer *na, integer *nra, integer *nca, doublereal *a, integer *ncb, doublereal *b, integer *nbc, integer *nrc, doublereal *c, doublereal *d, integer *irf, integer *icf)
{
  /* Aliases for the dimensions of the arrays */
  integer icf_dim1, irf_dim1;
  
  /* Local variables */
  integer i,j;
  integer nex;

  /*<	      NEX=NCA-2*NOV >*/
  irf_dim1 = *nra;
  icf_dim1 = *nca;

  /* Function Body */
  nex = *nca - (*nov << 1);
  if (nex == 0) {
    return 0;
  }
  
  /*     Initialization */
  for (i = 0; i <*na; ++i) {
    for (j = 0; j < *nra; ++j) {
      irf[j + i * irf_dim1] = j+1;
    }
    for (j = 0; j < *nca; ++j) {
      icf[j + i * icf_dim1] = j+1;
    }
  }

  switch(global_conpar_type) {
  default:
    conpar_default_wrapper(nov, na, nra, nca, a, 
			    ncb, b, nbc, nrc, c, d,irf, icf);
    break;
  }
  return 0;
} 








