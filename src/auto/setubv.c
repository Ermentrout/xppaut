#include "auto_f2c.h"
#include "auto_c.h"
#include "auto_types.h"

#ifdef PTHREADS
pthread_mutex_t mutex_for_funi = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef TIME
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
static double time_start (void) {
  struct timeval time;
  double seconds,microseconds;
  gettimeofday(&time,NULL);
  seconds = (double)time.tv_sec;
  microseconds = (double)time.tv_usec;
  return seconds + microseconds/1e6;
}
static double time_end(double start) {
  struct timeval time;
  double seconds,microseconds;
  gettimeofday(&time,NULL);
  seconds = (double)time.tv_sec;
  microseconds = (double)time.tv_usec;
  return (seconds + microseconds/1e6)-start;
}
#endif

void *setubv_make_aa_bb_cc(void * arg)
{  
  /* System generated locals */
  integer aa_dim1, aa_dim2, bb_dim1, bb_dim2, cc_dim1,
    cc_dim2, ups_dim1, 
    uoldps_dim1, udotps_dim1, upoldp_dim1,
    dbc_dim1, dicd_dim1, wploc_dim1, dfdu_dim1, dfdp_dim1, wp_dim1, wt_dim1;
  
  /* Local variables */
  integer i, j, k, l, m;
  integer k1, l1;
  integer i1,j1;

  integer ib, ic, jj;
  doublereal dt;  
  integer ib1, ic1;
  integer jp1;
  doublereal ddt;

  setubv_parallel_arglist *larg =  (setubv_parallel_arglist *)arg;

  doublereal *dicd, *ficd, *dfdp, *dfdu, *uold;
  doublereal *f;
  doublereal *u, *wploc;
  doublereal *dbc, *fbc, *uic, *uio, *prm, *uid, *uip, *ubc0, *ubc1;
  
  doublereal *ups = larg->ups;
  doublereal *upoldp = larg->upoldp;
  doublereal *udotps = larg->udotps;
  doublereal *uoldps = larg->uoldps;

  doublereal *aa = larg->aa;
  doublereal *bb = larg->bb;
  doublereal *cc = larg->cc;

  doublereal *wp = larg->wp;
  doublereal *wt = larg->wt;


  dicd = (doublereal *)malloc(sizeof(doublereal)*(larg->nint)*(larg->ndim + NPARX));
  ficd = (doublereal *)malloc(sizeof(doublereal)*(larg->nint));
  dfdp = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim)*NPARX);
  dfdu = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim)*(larg->ndim));
  uold = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));
  f    = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));
  u    = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));
  wploc= (doublereal *)malloc(sizeof(doublereal)*(larg->ncol)*(larg->ncol+1));
  dbc  = (doublereal *)malloc(sizeof(doublereal)*(larg->nbc)*(2*larg->ndim + NPARX));
  fbc  = (doublereal *)malloc(sizeof(doublereal)*(larg->nbc));
  uic  = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));
  uio  = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));
  prm  = (doublereal *)malloc(sizeof(doublereal)*NPARX);
  uid  = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));
  uip  = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));
  ubc0 = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));
  ubc1 = (doublereal *)malloc(sizeof(doublereal)*(larg->ndim));

  upoldp_dim1 = larg->ndxloc;
  udotps_dim1 = larg->ndxloc;
  uoldps_dim1 = larg->ndxloc;
  ups_dim1 = larg->ndxloc;
  dicd_dim1 = larg->nint;
  dbc_dim1 = larg->nbc;
  dfdu_dim1 = larg->ndim;
  dfdp_dim1 = larg->ndim;

  bb_dim1 = larg->ncb;
  bb_dim2 = larg->nra;
  
  cc_dim1 = larg->nca;
  cc_dim2 = larg->nrc;
  
  aa_dim1 = larg->nca;
  aa_dim2 = larg->nra;

  wploc_dim1 = larg->ncol + 1;
  wp_dim1 = larg->ncol + 1;
  wt_dim1 = larg->ncol + 1;

  /* Generate AA and BB: */
  
  /*      Partition the mesh intervals */
  /*jj will be replaced with loop_start and loop_end*/
  for (jj = larg->loop_start; jj < larg->loop_end; ++jj) {
    j = jj;
    jp1 = j + 1;
    dt = larg->dtm[j];
    ddt = 1. / dt;
    for (ic = 0; ic < larg->ncol; ++ic) {
      for (ib = 0; ib < larg->ncol + 1; ++ib) {
	ARRAY2D(wploc, ib, ic) = ddt * ARRAY2D(wp,ib, ic);
      }
    }
    /*this loop uses the loop_offset variable since up and uoldps
      and sent by the MPI version in their entirety, but
      loop_start and loop_end have been shifted.  The loop_offset
      variable contains the original value of loop_start and removes
      the shift*/
    for (ic = 0; ic < larg->ncol; ++ic) {
      for (k = 0; k < larg->ndim; ++k) {
	u[k] = ARRAY2D(wt, larg->ncol, ic) * ARRAY2D(ups, jp1 + larg->loop_offset, k);
	uold[k] = ARRAY2D(wt, larg->ncol, ic) * ARRAY2D(uoldps, jp1 + larg->loop_offset, k);
	for (l = 0; l < larg->ncol; ++l) {
	  l1 = l * larg->ndim + k;
	  u[k] += ARRAY2D(wt, l, ic) * ARRAY2D(ups, j + larg->loop_offset, l1);
	  uold[k] += ARRAY2D(wt, l, ic) * ARRAY2D(uoldps, j + larg->loop_offset, l1);
	}
      }

      for (i = 0; i < NPARX; ++i) {
	prm[i] = larg->par[i];
      }
      /*  
	  Ok this is a little wierd, so hold tight.  This function
	  is actually a pointer to a wrapper function, which eventually
	  calls the user defined func_.  Which wrapper is used
	  depends on what kind of problem it is.  The need for
	  the mutex is because some of these wrappers use a common
	  block for temporary storage 
	  NOTE!!!:  The icni and bcni wrappers do the same thing,
	  so if they ever get parallelized they need to be
	  checked as well.
      */
#ifdef PTHREADS
      pthread_mutex_lock(&mutex_for_funi);
#endif
      (*(larg->funi))(larg->iap, larg->rap, larg->ndim, u, uold, larg->icp, prm, 2, f, dfdu, dfdp);
#ifdef PTHREADS
      pthread_mutex_unlock(&mutex_for_funi);
#endif

      ic1 = ic * (larg->ndim);
      for (ib = 0; ib < larg->ncol + 1; ++ib) {
	double wt_tmp=ARRAY2D(wt, ib, ic);
	double wploc_tmp=ARRAY2D(wploc, ib, ic);
	ib1 = ib * larg->ndim;
	for (i = 0; i < larg->ndim; ++i) {
	  ARRAY3D(aa, ib1 + i, ic1 + i, jj) = wploc_tmp;
	  for (k = 0; k < larg->ndim; ++k) {
	    ARRAY3D(aa, ib1 + k, ic1 + i, jj) -= wt_tmp * ARRAY2D(dfdu, i, k);
	  }
	}
      }
      for (i = 0; i < larg->ndim; ++i) {
	for (k = 0; k < larg->ncb; ++k) {
	  ARRAY3D(bb, k, ic1 + i,  jj) = -ARRAY2D(dfdp, i, larg->icp[k]);
	}
      }
    }
  
  }

  /*     Generate CC : */
  
  /*     Boundary conditions : */
  if (larg->nbc > 0) {
    for (i = 0; i < larg->ndim; ++i) {
      ubc0[i] = ARRAY2D(ups, 0, i);
      ubc1[i] = ARRAY2D(ups, larg->na, i);
    }
    
#ifdef PTHREADS
    pthread_mutex_lock(&mutex_for_funi);
#endif
    (*(larg->bcni))(larg->iap, larg->rap, larg->ndim, larg->par, 
	    larg->icp, larg->nbc, ubc0, ubc1, fbc, 2, dbc);
#ifdef PTHREADS
    pthread_mutex_unlock(&mutex_for_funi);
#endif
    for (i = 0; i < larg->nbc; ++i) {
      for (k = 0; k < larg->ndim; ++k) {
	/*NOTE!!
	  This needs to split up.  Only the first processor does the first part
	  and only the last processors does the last part.*/
	if(larg->loop_offset + larg->loop_start == 0) {
	  ARRAY3D(cc, k, i, 0) = ARRAY2D(dbc, i, k);
	}
	if(larg->loop_offset + larg->loop_end == larg->na) {
	  ARRAY3D(cc, larg->nra + k, i, larg->na-1 - larg->loop_offset) = 
	    ARRAY2D(dbc ,i , larg->ndim + k);
	}
      }
    }
  }
  
  /*     Integral constraints : */
  if (larg->nint > 0) {
    for (jj = larg->loop_start; jj < larg->loop_end; ++jj) {
      j = jj;
      jp1 = j + 1;
      for (k = 0; k < (larg->ncol + 1); ++k) {
	for (i = 0; i < larg->ndim; ++i) {
	  i1 = k * larg->ndim + i;
	  j1 = j;
	  if (k+1 == (larg->ncol + 1)) {
	    i1 = i;
	  }
	  if (k+1 == (larg->ncol + 1)) {
	    j1 = jp1;
	  }
	  uic[i] = ARRAY2D(ups, j1 + larg->loop_offset, i1);
	  uio[i] = ARRAY2D(uoldps, j1 + larg->loop_offset, i1);
	  uid[i] = ARRAY2D(udotps, j1 + larg->loop_offset, i1);
	  uip[i] = ARRAY2D(upoldp, j1 + larg->loop_offset, i1);
	}
	
#ifdef PTHREADS
	pthread_mutex_lock(&mutex_for_funi);
#endif
	(*(larg->icni))(larg->iap, larg->rap, larg->ndim, larg->par, 
		larg->icp, larg->nint, 
		uic, uio, uid, uip, ficd, 2, dicd);
#ifdef PTHREADS
	pthread_mutex_unlock(&mutex_for_funi);
#endif
	
	for (m = 0; m < larg->nint; ++m) {
	  for (i = 0; i < larg->ndim; ++i) {
	    k1 = k * larg->ndim + i;
	    ARRAY3D(cc, k1, larg->nbc + m, jj) = 
	      larg->dtm[j] * larg->wi[k ] * ARRAY2D(dicd, m, i);
	  }
	}
      }
    }
  }
  /*     Pseudo-arclength equation : */
  for (jj = larg->loop_start; jj < larg->loop_end; ++jj) {
    for (i = 0; i < larg->ndim; ++i) {
      for (k = 0; k < larg->ncol; ++k) {
	k1 = k * larg->ndim + i;
	ARRAY3D(cc, k1 , larg->nrc - 1, jj) = 
	  larg->dtm[jj] * larg->thu[i] * larg->wi[k] * 
	  ARRAY2D(udotps, jj + larg->loop_offset, k1);
      }
      ARRAY3D(cc, larg->nra + i, larg->nrc -1, jj) = 
	larg->dtm[jj] * larg->thu[i] * larg->wi[larg->ncol] * 
	ARRAY2D(udotps, jj + 1 + larg->loop_offset, i);
    }
  }

  free(dicd );
  free(ficd );
  free(dfdp );
  free(dfdu );
  free(uold );
  free(f    );
  free(u    );
  free(wploc);
  free(dbc  );
  free(fbc  );
  free(uic  );
  free(uio  );
  free(prm  );
  free(uid  );
  free(uip  );
  free(ubc0 );
  free(ubc1 );

  return NULL;

}

#ifdef PTHREADS
int 
setubv_threads_wrapper(setubv_parallel_arglist data)
{
  setubv_parallel_arglist *send_data;
  int i;
  pthread_t *th;
  void * retval;
  pthread_attr_t attr;
  int retcode;
  
  th = (pthread_t *)malloc(sizeof(pthread_t)*global_num_procs);
  send_data = (setubv_parallel_arglist *)malloc(sizeof(setubv_parallel_arglist)*global_num_procs);
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);

  for(i=0;i<global_num_procs;i++) {
    setubv_parallel_arglist_copy(&send_data[i],data);
    send_data[i].loop_start = (i*(data.na))/global_num_procs;
    send_data[i].loop_end = ((i+1)*(data.na))/global_num_procs;
    send_data[i].loop_offset = 0;
    retcode = pthread_create(&th[i], &attr, setubv_make_aa_bb_cc, (void *) &send_data[i]);
    if (retcode != 0) fprintf(stderr, "create %d failed %d\n", i, retcode);
  }
  for(i=0;i<global_num_procs;i++) {
    retcode = pthread_join(th[i], &retval);
    if (retcode != 0) fprintf(stderr, "join %d failed %d\n", i, retcode);
  }  
  free(send_data);
  free(th);
  return 0;
}
#endif

#ifdef MPI
int 
setubv_mpi_wrapper(setubv_parallel_arglist data)
{
  integer loop_start,loop_end;
  integer loop_start_tmp,loop_end_tmp;
  integer loop_offset;
  int i,comm_size;
  int *aa_counts,*aa_displacements;
  int *bb_counts,*bb_displacements;
  int *cc_counts,*cc_displacements;
  int *dtm_counts,*dtm_displacements;

  MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
  aa_counts=(int *)malloc(sizeof(int)*comm_size);
  aa_displacements=(int *)malloc(sizeof(int)*comm_size);
  bb_counts=(int *)malloc(sizeof(int)*comm_size);
  bb_displacements=(int *)malloc(sizeof(int)*comm_size);
  cc_counts=(int *)malloc(sizeof(int)*comm_size);
  cc_displacements=(int *)malloc(sizeof(int)*comm_size);
  dtm_counts=(int *)malloc(sizeof(int)*comm_size);
  dtm_displacements=(int *)malloc(sizeof(int)*comm_size);
  aa_counts[0] = 0;
  aa_displacements[0] = 0;
  bb_counts[0] = 0;
  bb_displacements[0] = 0;
  cc_counts[0] = 0;
  cc_displacements[0] = 0;
  dtm_counts[0] = 0;
  dtm_displacements[0] = 0;

  
  for(i=1;i<comm_size;i++){
    
    /*Send message to get worker into setubv mode*/
    {
      int message=AUTO_MPI_SETUBV_MESSAGE;
      MPI_Send(&message,1,MPI_INT,i,0,MPI_COMM_WORLD);
    }
    loop_start = ((i-1)*(data.na))/(comm_size - 1);
    loop_end = ((i)*(data.na))/(comm_size - 1);
    aa_counts[i] = (data.nca)*(data.nra)*(loop_end-loop_start);
    aa_displacements[i] = (data.nca)*(data.nra)*loop_start;
    bb_counts[i] = (data.ncb)*(data.nra)*(loop_end-loop_start);
    bb_displacements[i] = (data.ncb)*(data.nra)*loop_start;
    cc_counts[i] = (data.nca)*(data.nrc)*(loop_end-loop_start);
    cc_displacements[i] = (data.nca)*(data.nrc)*loop_start;
    dtm_counts[i] = (loop_end-loop_start);
    dtm_displacements[i] = (loop_start);

    loop_start_tmp = 0;
    loop_end_tmp = loop_end-loop_start;
    MPI_Send(&loop_start_tmp ,1,MPI_LONG,i,0,MPI_COMM_WORLD);
    MPI_Send(&loop_end_tmp   ,1,MPI_LONG,i,0,MPI_COMM_WORLD);
    loop_offset = loop_start;
    MPI_Send(&loop_offset    ,1,MPI_LONG,i,0,MPI_COMM_WORLD);
  }

  {
    integer params[11];
    params[0]=data.na;
    params[1]=data.ndim;
    params[2]=data.ips;
    params[3]=data.ncol;
    params[4]=data.nbc;
    params[5]=data.nint;
    params[6]=data.ncb;
    params[7]=data.nrc;
    params[8]=data.nra;
    params[9]=data.nca;
    params[10]=data.ndxloc;
    MPI_Bcast(params     ,11,MPI_LONG,0,MPI_COMM_WORLD);
  }    

  {
    int position=0;
    void *buffer;
    int bufsize;
    int size_int,size_double;
    int niap,nrap;
    /* Here we compute the number of elements in the iap and rap structures.
       Since each of the structures is homogeneous we just divide the total
       size by the size of the individual elements.*/
    niap = sizeof(iap_type)/sizeof(integer);
    nrap = sizeof(rap_type)/sizeof(doublereal);
    MPI_Pack_size(niap+NPARX2,MPI_LONG,MPI_COMM_WORLD,&size_int);
    MPI_Pack_size(nrap+NPARX2+
		  (data.ndxloc)*(data.ndim)*(data.ncol)+
		  (data.ndxloc)*(data.ndim)*(data.ncol)+
		  (data.ncol + 1)*(data.ncol)+
		  (data.ncol + 1)*(data.ncol)+
		  (data.ncol + 1)+
		  (data.ndxloc)*(data.ndim)*(data.ncol)+
		  (data.ndxloc)*(data.ndim)*(data.ncol)+
		  (data.ndim)*8+
		  NPARX+
		  NPARX,
		  MPI_DOUBLE,MPI_COMM_WORLD,&size_double);
    bufsize = size_int + size_double;
    buffer=malloc((unsigned)bufsize);

    MPI_Pack(data.iap    ,niap,MPI_LONG,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.rap    ,nrap,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    /**********************************************/
    MPI_Pack(data.par    ,NPARX2,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.icp    ,NPARX2,MPI_LONG,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.ups    ,(data.ndxloc)*(data.ndim)*(data.ncol),MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.uoldps ,(data.ndxloc)*(data.ndim)*(data.ncol),MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.wp     ,(data.ncol + 1)*(data.ncol),MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.wt     ,(data.ncol + 1)*(data.ncol),MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.wi     ,(data.ncol + 1),MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.udotps ,(data.ndxloc)*(data.ndim)*(data.ncol),MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.upoldp ,(data.ndxloc)*(data.ndim)*(data.ncol),MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);

    MPI_Pack(data.thu    ,(data.ndim)*8,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.thl    ,NPARX,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(data.rldot  ,NPARX,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    
    MPI_Bcast(buffer     ,position,MPI_PACKED,0,MPI_COMM_WORLD);
  }

  MPI_Scatterv(data.dtm        ,dtm_counts,dtm_displacements,MPI_DOUBLE,
	       NULL,0,MPI_DOUBLE,
	       0,MPI_COMM_WORLD);

  /* Worker runs here */
  return 0;
}
#endif

int 
setubv_default_wrapper(setubv_parallel_arglist data)
{
  setubv_make_aa_bb_cc((void *)&data);
  return 0;
}

int 
setubv(integer ndim, integer ips, integer na, integer ncol, integer nbc, integer nint, integer ncb, integer nrc, integer nra, integer nca, 
       FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), integer ndxloc, iap_type *iap, rap_type *rap, doublereal *par, integer *icp, 
       doublereal rds, doublereal *aa, doublereal *bb, doublereal *cc, doublereal *dd, doublereal *fa, doublereal *fc, doublereal *rlcur, 
       doublereal *rlold, doublereal *rldot, doublereal *ups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *dups, 
       doublereal *dtm, doublereal *thl, doublereal *thu, doublereal *p0, doublereal *p1)
{
  /* System generated locals */
  integer aa_dim1, aa_dim2, bb_dim1, bb_dim2, cc_dim1,
    cc_dim2, dd_dim1;

  /* Local variables */
  integer i, j, k;

  doublereal *wi, *wp, *wt;
  
  wi   = (doublereal *)malloc(sizeof(doublereal)*(ncol+1) );
  wp   = (doublereal *)malloc(sizeof(doublereal)*(ncol)*(ncol+1) );
  wt   = (doublereal *)malloc(sizeof(doublereal)*(ncol)*(ncol+1) );

  dd_dim1 = ncb;
  
  bb_dim1 = ncb;
  bb_dim2 = nra;
  
  cc_dim1 = nca;
  cc_dim2 = nrc;
  
  aa_dim1 = nca;
  aa_dim2 = nra;

  wint(ncol + 1, wi);
  genwts(ncol, ncol + 1, wt, wp);
  
  /* Initialize to zero. */
  for (i = 0; i < nrc; ++i) {
    fc[i] = 0.;
    for (k = 0; k < ncb; ++k) {
      ARRAY2D(dd, k, i) = 0.;
    }
  }

  /* Set constants. */
  for (i = 0; i < ncb; ++i) {
    par[icp[i]] = rlcur[i];
  }
  
  /*  NA is the local node's mesh interval number. */
  
  for (i = 0; i < na; ++i) {
    for (j = 0; j < nra; ++j) {
      for (k = 0; k < nca; ++k) {
	ARRAY3D(aa, k, j, i) = 0.;
      }
    }
    for (j = 0; j < nra; ++j) {
      for (k = 0; k < ncb; ++k) {
	ARRAY3D(bb, k, j, i) = 0.;
      }
    }
    for (j = 0; j < nca; ++j) {
      for (k = 0; k < nrc; ++k) {
	ARRAY3D(cc, j, k, i) = 0.;
      }
    }
  }

  /*     ** Time evolution computations (parabolic systems) */
  if (ips == 14 || ips == 16) {
    rap->tivp = rlold[0];
  } 
 
  {
    setubv_parallel_arglist arglist;
    setubv_parallel_arglist_constructor(ndim, ips, na, ncol, nbc, nint, ncb, 
					nrc, nra, nca, funi, icni, ndxloc, iap, rap, 
					par, icp, aa, bb, cc, dd, fa, fc, ups, 
					uoldps, udotps, upoldp, dtm, wp, wt, wi, 
					thu, thl, rldot, bcni, &arglist);
  
    switch(global_setubv_type) {

#ifdef PTHREADS
    case SETUBV_PTHREADS:
      setubv_threads_wrapper(arglist);
      break;
#endif

#ifdef MPI
    case SETUBV_MPI:
      if(global_verbose_flag)
	printf("Setubv MPI start\n");
      setubv_mpi_wrapper(arglist);
      if(global_verbose_flag)
	printf("Setubv MPI end\n");
      break;
#endif

    default:
      setubv_default_wrapper(arglist);
      break;
    }
    setubv_make_fa(arglist);
    setubv_make_fc_dd(arglist,dups,rlcur,rlold,rds);
  }

  free(wi   );
  free(wp   );
  free(wt   );
  return 0;
}

void setubv_make_fa(setubv_parallel_arglist larg) {
  integer i,j,k,l;
  integer ic,k1,ib;
  integer jj,jp1,l1,ic1;
  doublereal dt,ddt;

  doublereal *ups = larg.ups;
  integer ups_dim1 = larg.ndxloc;

  doublereal *uoldps = larg.uoldps;
  integer uoldps_dim1 = larg.ndxloc;

  doublereal *wp = larg.wp;
  integer wp_dim1 = larg.ncol + 1;

  doublereal *wt = larg.wt;
  integer wt_dim1 = larg.ncol + 1;
  
  doublereal *fa = larg.fa;
  integer fa_dim1 = larg.nra;
  
  doublereal *wploc= (doublereal *)malloc(sizeof(doublereal)*(larg.ncol)*(larg.ncol+1));
  integer wploc_dim1 = larg.ncol + 1;
  
  doublereal *dfdp = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim)*NPARX);
  doublereal *dfdu = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim)*(larg.ndim));
  doublereal *u    = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));
  doublereal *uold = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));
  doublereal *f    = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));
  doublereal *prm  = (doublereal *)malloc(sizeof(doublereal)*NPARX);

  for (jj = 0; jj < larg.na; ++jj) {
    j = jj;
    jp1 = j + 1;
    dt = larg.dtm[j];
    ddt = 1. / dt;
    for (ic = 0; ic < larg.ncol; ++ic) {
      for (ib = 0; ib < larg.ncol + 1; ++ib) {
	ARRAY2D(wploc, ib, ic) = ddt * ARRAY2D(wp,ib, ic);
      }
    }
    for (ic = 0; ic < larg.ncol; ++ic) {
      for (k = 0; k < larg.ndim; ++k) {
	u[k] = ARRAY2D(wt, larg.ncol, ic) * ARRAY2D(ups, jp1, k);
	uold[k] = ARRAY2D(wt, larg.ncol, ic) * ARRAY2D(uoldps, jp1, k);
	for (l = 0; l < larg.ncol; ++l) {
	  l1 = l * larg.ndim + k;
	  u[k] += ARRAY2D(wt, l, ic) * ARRAY2D(ups, j + larg.loop_offset, l1);
	  uold[k] += ARRAY2D(wt, l, ic) * ARRAY2D(uoldps, j + larg.loop_offset, l1);
	}
      }

      for (i = 0; i < NPARX; ++i) {
	prm[i] = larg.par[i];
      }
      (*(larg.funi))(larg.iap, larg.rap, larg.ndim, u, uold, larg.icp, prm, 2, f, dfdu, dfdp);

      ic1 = ic * (larg.ndim);
      for (i = 0; i < larg.ndim; ++i) {
	ARRAY2D(fa,ic1 + i, jj) = f[i] - ARRAY2D(wploc, larg.ncol, ic) * ARRAY2D(ups, jp1 + larg.loop_offset, i);
	for (k = 0; k < larg.ncol; ++k) {
	  k1 = k * larg.ndim + i;
	  ARRAY2D(fa, ic1 + i, jj) -= ARRAY2D(wploc, k, ic) * ARRAY2D(ups, j + larg.loop_offset, k1);
	}
      }
    }
  
  }
  free(wploc);
  free(dfdp);
  free(dfdu);
  free(u);
  free(uold);
  free(f);
  free(prm);
  
}


void setubv_make_fc_dd(setubv_parallel_arglist larg, doublereal *dups, doublereal *rlcur, 
	     doublereal *rlold, doublereal rds) {
  integer i,j,jj,jp1,k,i1,m,j1;
  doublereal rlsum;

  integer dups_dim1 = larg.ndxloc;
  
  doublereal *dd = larg.dd;
  integer dd_dim1 = larg.ncb;

  doublereal *ups = larg.ups;
  integer ups_dim1 = larg.ndxloc;

  doublereal *uoldps = larg.uoldps;
  integer uoldps_dim1 = larg.ndxloc;
  
  doublereal *udotps = larg.udotps;
  integer udotps_dim1 = larg.ndxloc;
  
  doublereal *upoldp = larg.upoldp;
  integer upoldp_dim1 = larg.ndxloc;
  
  integer dbc_dim1 = larg.nbc;
  doublereal *dbc  = (doublereal *)malloc(sizeof(doublereal)*(larg.nbc)*(2*larg.ndim + NPARX));
  doublereal *fbc  = (doublereal *)malloc(sizeof(doublereal)*(larg.nbc));
  doublereal *ubc0 = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));
  doublereal *ubc1 = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));
  integer dicd_dim1 = larg.nint;
  doublereal *dicd = (doublereal *)malloc(sizeof(doublereal)*(larg.nint)*(larg.ndim + NPARX));
  doublereal *ficd = (doublereal *)malloc(sizeof(doublereal)*(larg.nint));
  doublereal *uic  = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));
  doublereal *uio  = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));
  doublereal *uid  = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));
  doublereal *uip  = (doublereal *)malloc(sizeof(doublereal)*(larg.ndim));

  /* Boundary condition part of FC */
  if (larg.nbc > 0) {
    for (i = 0; i < larg.ndim; ++i) {
      ubc0[i] = ARRAY2D(ups, 0, i);
      ubc1[i] = ARRAY2D(ups, larg.na, i);
    }
    
    (*(larg.bcni))(larg.iap, larg.rap, larg.ndim, larg.par, 
		   larg.icp, larg.nbc, ubc0, ubc1, fbc, 2, dbc);
    for (i = 0; i < larg.nbc; ++i) {
      larg.fc[i] = -fbc[i];
      for (k = 0; k < larg.ncb; ++k) {
	ARRAY2D(dd, k, i) = 
	  ARRAY2D(dbc, i, (larg.ndim *2) + larg.icp[k]);
      }
    }
    /*       Save difference : */
    for (j = 0; j < larg.na + 1; ++j) {
      for (i = 0; i < larg.nra; ++i) {
	ARRAY2D(dups, j, i) = ARRAY2D(ups, j, i) - ARRAY2D(uoldps, j, i);
      }
    }
  }

  /* Integral constraint part of FC */
  if (larg.nint > 0) {
    for (jj = larg.loop_start; jj < larg.loop_end; ++jj) {
      j = jj;
      jp1 = j + 1;
      for (k = 0; k < (larg.ncol + 1); ++k) {
	for (i = 0; i < larg.ndim; ++i) {
	  i1 = k * larg.ndim + i;
	  j1 = j;
	  if (k+1 == (larg.ncol + 1)) {
	    i1 = i;
	  }
	  if (k+1 == (larg.ncol + 1)) {
	    j1 = jp1;
	  }
	  uic[i] = ARRAY2D(ups, j1, i1);
	  uio[i] = ARRAY2D(uoldps, j1, i1);
	  uid[i] = ARRAY2D(udotps, j1, i1);
	  uip[i] = ARRAY2D(upoldp, j1, i1);
	}
	
	(*(larg.icni))(larg.iap, larg.rap, larg.ndim, larg.par, 
		larg.icp, larg.nint, 
		uic, uio, uid, uip, ficd, 2, dicd);
	
	for (m = 0; m < larg.nint; ++m) {
	  larg.fc[larg.nbc + m] -= larg.dtm[j] * larg.wi[k] * ficd[m];
	  for (i = 0; i < larg.ncb; ++i) {
	    ARRAY2D(dd, i, larg.nbc + m) += 
	      larg.dtm[j] * larg.wi[k] * ARRAY2D(dicd, m, larg.ndim + larg.icp[i]);
	  }
	}
      }
    }
  }

  for (i = 0; i < larg.ncb; ++i) {
    ARRAY2D(dd, i, (larg.nrc-1)) = larg.thl[larg.icp[i]] * larg.rldot[i];
  }

  rlsum = 0.;
  for (i = 0; i < larg.ncb; ++i) {
    rlsum += larg.thl[larg.icp[i]] * (rlcur[i] - rlold[i]) * larg.rldot[i];
  }

  larg.fc[larg.nrc-1] = rds - rinpr(larg.iap, &(larg.ndim), &(larg.ndxloc), larg.udotps, dups, larg.dtm, larg.thu) - rlsum;

  free(dbc);
  free(fbc);
  free(ubc0);
  free(ubc1);
  free(dicd);
  free(ficd);
  free(uic);
  free(uio);
  free(uid);
  free(uip);

}

/* Copy a setubv_parallel_arglist */
void setubv_parallel_arglist_copy(setubv_parallel_arglist *output,
				  const setubv_parallel_arglist input) {
  memcpy(output,&input,sizeof(setubv_parallel_arglist));
}


/* Fill in a setubv_parallel_arglist for the individual variables */
void setubv_parallel_arglist_constructor(integer ndim, integer ips, integer na, integer ncol, 
					 integer nbc, integer nint, integer ncb, integer nrc, integer nra, integer nca, 
					 FUNI_TYPE((*funi)), ICNI_TYPE((*icni)), integer ndxloc, iap_type *iap, rap_type *rap, doublereal *par, 
					 integer *icp, doublereal *aa, doublereal *bb, 
					 doublereal *cc, doublereal *dd, doublereal *fa, doublereal *fc, doublereal *ups, 
					 doublereal *uoldps, doublereal *udotps, 
					 doublereal *upoldp, doublereal *dtm, 
					 doublereal *wp, doublereal *wt, doublereal *wi,
					 doublereal *thu, doublereal *thl,
					 doublereal *rldot, BCNI_TYPE((*bcni)), setubv_parallel_arglist *data) {
  data->ndim   = ndim;
  data->ips    = ips;
  data->ncol   = ncol;
  data->nbc    = nbc;
  data->nint   = nint;
  data->ncb    = ncb;
  data->nrc    = nrc;
  data->nra    = nra;
  data->nca    = nca;
  data->na     = na;
  data->funi   = funi;
  data->icni   = icni;
  data->ndxloc = ndxloc;
  data->iap    = iap;
  data->rap    = rap;
  data->par    = par;
  data->icp    = icp;
  data->aa     = aa;
  data->bb     = bb;
  data->cc     = cc;
  data->dd     = dd;
  data->fa     = fa;
  data->fc     = fc;
  data->ups    = ups;
  data->uoldps = uoldps;
  data->udotps = udotps;
  data->upoldp = upoldp;
  data->dtm    = dtm;
  data->loop_start = 0;
  data->loop_end   = na;
  data->loop_offset = 0;
  data->wp     = wp;
  data->wt     = wt;
  data->wi     = wi;
  data->thu    = thu;
  data->thl    = thl;
  data->rldot  = rldot;
  data->bcni   = bcni;
}  









