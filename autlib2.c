/* Autlib2.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "auto_f2c.h"
#include "auto_c.h"
/* #include "malloc.h" */

/*#define ACCES_TEST*/
#ifdef ACCES_TEST
struct {
  double *a,*b,*c;
} test ;
#endif

typedef struct {
  doublereal *a;
  doublereal *b;
  doublereal *c;
  doublereal *d;
  doublereal *a1;
  doublereal *a2;
  doublereal *s1;
  doublereal *s2;
  doublereal *bb;
  doublereal *cc;
  doublereal *faa;
  doublereal *ca1;

  integer *icf;
  integer *irf;
  integer *ipr;
  integer *icf11;
  integer *icf1;
  integer *icf2;
  integer *np;
} main_auto_storage_type;

 main_auto_storage_type main_auto_storage={NULL,NULL,NULL,NULL,
						   NULL,NULL,NULL,NULL,
						   NULL,NULL,NULL,NULL,
						   NULL,NULL,NULL,NULL,
						   NULL,NULL,NULL};

void print_jacobian(iap_type iap,main_auto_storage_type data) {
  int i,j,k,l;
  int num_rows_A = iap.ndim * iap.ncol;
  int num_columns_A = iap.ndim * (iap.ncol + 1);
  int num_columns_B = iap.nfpr;
  int num_rows_C = iap.nbc + iap.nint + 1;
  int numblocks = iap.ntst;
  FILE *fp;
  static int num_calls=0;
  char filename[80];

  sprintf(filename,"jacobian%03d",num_calls);
  fp=fopen(filename,"w");
  num_calls++;

  for(i=0;i<numblocks;i++){
    for(j=0;j<num_rows_A;j++){
      /* Print zeros in front first */
      for(k=0;k<i*(num_columns_A-iap.ndim);k++)
	fprintf(fp,"%18.10e ",0.0);
      /* Now print line from block */
      for(k=0;k<num_columns_A;k++)
	fprintf(fp,"%18.10e ",data.a[k + j*num_columns_A + i*num_rows_A*num_columns_A]);
      /* Now put zeros at end of line */
      for(k=i*(num_columns_A-iap.ndim)+num_columns_A;k<(num_columns_A-iap.ndim)*numblocks+iap.ndim;k++)
	fprintf(fp,"%18.10e ",0.0);
      /* Put in B */
      for(k=0;k<num_columns_B;k++)
	fprintf(fp,"%18.10e ",data.b[k + j*num_columns_B + i*num_rows_A*num_columns_B]);
      fprintf(fp,"\n");
    }
  }

  /*For printing out C there needs to be a summation of the edge guys*/
  for(j=0;j<num_rows_C;j++) {
    /*The first num_rows_A columns are ok as the are*/
    for(k=0;k<(num_columns_A-iap.ndim);k++)
      fprintf(fp,"%18.10e ",data.c[k + j*num_columns_A + 0*num_rows_C*num_columns_A]);
    /* Now print out the rest of the blocks, doing a summation at the beginning of each */
    for(i=1;i<numblocks;i++) {
      for(k=0;k<iap.ndim;k++)
	fprintf(fp,"%18.10e ",data.c[k+ num_columns_A-iap.ndim + j*num_columns_A + (i-1)*num_rows_C*num_columns_A] +
		data.c[k + j*num_columns_A + i*num_rows_C*num_columns_A]);
      for(k=iap.ndim;k<num_columns_A-iap.ndim;k++)
	fprintf(fp,"%18.10e ",data.c[k + j*num_columns_A + i*num_rows_C*num_columns_A]);
    }
    /*Now print out last column*/
    for(k=num_columns_A-iap.ndim;k<num_columns_A;k++)
      fprintf(fp,"%18.10e ",data.c[k + j*num_columns_A + (numblocks-1)*num_rows_C*num_columns_A]);
    for(l=0;l<num_columns_B;l++)
      fprintf(fp,"%18.10e ",data.d[l + j*num_columns_B]);
    fprintf(fp,"\n");
  }


  fclose(fp);

}

void print_ups_rlcur(iap_type iap,doublereal *ups,doublereal *rlcur) {
  FILE *fp;
  static int num_calls=0;
  char filename[80];
  int i;
  
  sprintf(filename,"ups_rlcur%03d",num_calls);
  fp=fopen(filename,"w");
  num_calls++;
  for(i=0;i<(iap.ndim)*(iap.ncol)*(iap.ntst) + iap.ndim;i++)
    fprintf(fp,"%18.10e\n",ups[i]);
  for(i=0;i<iap.nfpr;i++)
    fprintf(fp,"%18.10e\n",rlcur[i]);

  fclose(fp);

}

void print_fa_fc(iap_type iap,doublereal *fa,doublereal *fc,char *filename) {
  FILE *fp;
  int i,j;
  int num_rows_A = iap.ndim * iap.ncol;
  int numblocks = iap.ntst;

  fp=fopen(filename,"w");

  for(i=0;i<numblocks;i++)
    for(j=0;j<num_rows_A;j++)
      fprintf(fp,"%18.10e\n",fa[j+i*num_rows_A]);
  for(i=0;i<iap.nfpr+iap.ndim;i++)
    fprintf(fp,"%10.10e\n",fc[i]);

  fclose(fp);

}

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*           Setting up of the Jacobian and right hand side */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ------ */
/* Subroutine */ int 
solvbv(integer *ifst, iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), doublereal *rds, integer *nllv, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *dups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *dtm, doublereal *fa, doublereal *fc, doublereal *p0, doublereal *p1, doublereal *thl, doublereal *thu)
{
  
  
  /* Local variables */
  
  integer ndim;
  logical ipar;
  integer ncol, nclm, nfpr, nint, nrow, ntst, ntst0;
  
  doublereal *ff, *ft;
  
  integer nbc, iid, iam;
  doublereal det;
  integer ips, nrc;
  
  integer kwt;
  
 

  /*     N AX is the local N TSTX, which is smaller than the global N TSTX. */
  /*     NODES is the total number of nodes. */
  
  
  /* Sets up and solves the linear equations for one Newton/Chord iteration 
   */
  
  
  /* Most of the required memory is allocated below */
  /* This is an interesting section of code.  The main point
     is that setubv and conpar only get called when ifst
     is 1.  This is a optimization since you can solve
     the system using the previously factored jacobian.
     One thing to watch out for is that two seperate calls
     of solvbv_ talk to each other through these arrays,
     so it is only safe to get rid of them when ifst is
     1 (since their entries will then be recreated in conpar
     and setubv).
  */

  ff = (double *)malloc(sizeof(double)*((iap->ndim * iap->ncol ) )* iap->ntst + 1);
  ft = (double *)malloc(sizeof(double)*((iap->ndim * iap->ncol ) )* (iap->ntst + 1));



   if (*ifst==1){
     /* printf("I am freeing and allocating stuff \n");  */
    /* The formulas used for the allocation are somewhat complex, but they
       are based on following macros (the space after the first letter is 
       for the scripts which detect these things automatically, the original
       name does not have the space:
       
       M 1AAR =  (((iap->ndim * iap->ncol ) + iap->ndim ) )      
       M 2AA  =	((iap->ndim * iap->ncol ) )                     
       N AX   =	(iap->ntst /NODES+1)                            
       M 1BB  =	(NPARX)                                         
       M 2BB  =	((iap->ndim * iap->ncol ) )                     
       M 1CC  =	((((iap->ndim * iap->ncol ) + iap->ndim ) ) )   
       M 2CC  =	(((iap->ndim +3) +NINTX+1) )                    
       M 1DD  =	(((iap->ndim +3) +NINTX+1) )                    
       M 2DD  =	(NPARX)                                         
       N RCX  =	((iap->ndim +3) +NINTX+1)                       
       N CLMX =	((iap->ndim * iap->ncol ) + iap->ndim )         
       N ROWX =	(iap->ndim * iap->ncol )                        
    */
    
    /* Free floating point arrays */
    free(main_auto_storage.a);
    free(main_auto_storage.b);
    free(main_auto_storage.c);
    free(main_auto_storage.d);
    free(main_auto_storage.a1);
    free(main_auto_storage.a2);
    free(main_auto_storage.s1);
    free(main_auto_storage.s2);
    free(main_auto_storage.bb);
    free(main_auto_storage.cc);
    free(main_auto_storage.faa);
    free(main_auto_storage.ca1);
    
    /* Free integer arrays */
    free(main_auto_storage.icf);
    free(main_auto_storage.irf);
    free(main_auto_storage.ipr);
    free(main_auto_storage.icf11);
    free(main_auto_storage.icf1);
    free(main_auto_storage.icf2);
    free(main_auto_storage.np);

    /*(M 1AAR*M 2AA*N AX) */
    main_auto_storage.a=(doublereal *)malloc(sizeof(doublereal)*((((iap->ndim * iap->ncol ) + iap->ndim ) ) * 
								 ((iap->ndim * iap->ncol ) ) * 
								 (iap->ntst +1) )); 
    /*(M 1BB*M 2BB*N AX)*/ 
    main_auto_storage.b=(doublereal *)malloc(sizeof(doublereal)*((NPARX) * ((iap->ndim * iap->ncol ) ) * (iap->ntst +1) ) );
    /*(M 1CC*M 2CC*N AX)*/ 
    main_auto_storage.c=(doublereal *)malloc(sizeof(doublereal)*(((((iap->ndim * iap->ncol ) + iap->ndim ) ) ) * 
								   ((iap->nbc +iap->nint+1) ) * (iap->ntst +1) ));
    /*(M 1DD*M 2DD)*/ 
    main_auto_storage.d=(doublereal *)malloc(sizeof(doublereal)*(((iap->nbc +iap->nint+1) ) * (NPARX) ) );
    /*(iap->ndim * iap->ndim *N AX)*/ 
    main_auto_storage.a1=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim * iap->ndim * (iap->ntst +1) ) );
    /*(iap->ndim * iap->ndim *N AX)*/ 
    main_auto_storage.a2=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim * iap->ndim * (iap->ntst +1) )); 
    /*(iap->ndim * iap->ndim *N AX)*/ 
    main_auto_storage.s1=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim * iap->ndim * (iap->ntst +1) )); 
    /*(iap->ndim * iap->ndim *N AX)*/ 
    main_auto_storage.s2=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim * iap->ndim * (iap->ntst +1) )); 
    /*(iap->ndim *N PARX*N AX)*/ 
    main_auto_storage.bb=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim *NPARX* (iap->ntst +1) ) );
    /*(N RCX* iap->ndim *N AX+1)*/ 
    main_auto_storage.cc=(doublereal *)malloc(sizeof(doublereal)*((iap->nbc + iap->nint + 1) * iap->ndim * (iap->ntst  + 1) + 1));

    /*(iap->ndim *N AX)*/ 
    main_auto_storage.faa=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim * (iap->ntst +1) ) );

    /*(iap->ndim * iap->ndim *K REDO)*/ 
    main_auto_storage.ca1=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim * iap->ndim *KREDO) );
    
    /*(N CLMX*N AX)*/ 
    main_auto_storage.icf=(integer *)malloc(sizeof(integer)*(((iap->ndim * iap->ncol ) + iap->ndim ) * (iap->ntst +1) ) );
    /*(N ROWX*N AX)*/ 
    main_auto_storage.irf=(integer *)malloc(sizeof(integer)*((iap->ndim * iap->ncol ) * (iap->ntst +1) ) );
    /*(iap->ndim *N AX)*/ 
    main_auto_storage.ipr=(integer *)malloc(sizeof(integer)*(iap->ndim * (iap->ntst +1) ) );
    /*(iap->ndim *K REDO)*/ 
    main_auto_storage.icf11=(integer *)malloc(sizeof(integer)*(iap->ndim *KREDO) );
    /*(iap->ndim *N AX)*/ 
    main_auto_storage.icf1=(integer *)malloc(sizeof(integer)*(iap->ndim * (iap->ntst +1) ));
    /*(iap->ndim *N AX)*/ 
    main_auto_storage.icf2=(integer *)malloc(sizeof(integer)*(iap->ndim * (iap->ntst +1) )); 
    /*(2)*/ 
    main_auto_storage.np=(integer *)malloc(sizeof(integer)*(2) );
     } 
  
  
  iam = iap->mynode;
  kwt = iap->numnodes;
  if (kwt > 1) {
    ipar = TRUE_;
  } else {
    ipar = FALSE_;
  }
  
  ndim = iap->ndim;
  ips = iap->ips;
  ntst = iap->ntst;
  ncol = iap->ncol;
  nbc = iap->nbc;
  nint = iap->nint;
  iid = iap->iid;
  nfpr = iap->nfpr;
  nrc = nbc + nint + 1;
  nrow = ndim * ncol;
  nclm = nrow + ndim;
  
  if (kwt > ntst) {
    printf("NTST is less than the number of nodes\n");
    exit(0);
  } else {
    partition(&ntst, &kwt, main_auto_storage.np);
  }
  
  /*     NTST0 is the global one, NTST is the local one. */
  /*     The value of NTST may be different in different nodes. */
  ntst0 = ntst;
  ntst = main_auto_storage.np[iam];
  
  if (*ifst == 1) {
    setubv(ndim, ips, ntst, ncol, nbc, nint, nfpr, nrc, nrow, nclm,
	   funi, bcni, icni, *ndxloc, iap, rap, par, icp, 
	   *rds, main_auto_storage.a, main_auto_storage.b, main_auto_storage.c, main_auto_storage.d, ft, fc, rlcur, 
	   rlold, rldot, ups, uoldps, udotps, upoldp, dups, 
	   dtm, thl, thu, p0, p1);
#ifdef ACCES_TEST
    test.a=main_auto_storage.a;
    test.b=main_auto_storage.b;
    test.c=main_auto_storage.c;
    main_auto_storage.a = NULL;
    main_auto_storage.b = NULL;
    main_auto_storage.c = NULL;
#endif
  } else {
    setrhs(&ndim, &ips, &ntst, &ntst0, main_auto_storage.np, &ncol, &nbc, &nint, &
	   nfpr, &nrc, &nrow, &nclm, &iam, &kwt, &ipar, funi, bcni, icni,
	   ndxloc, iap, rap, par, icp, rds, ft, fc, rlcur, 
	   rlold, rldot, ups, uoldps, udotps, upoldp, dups, dtm, thl, 
	   thu, p0, p1);
  }
  /*     The matrix D and FC are set to zero for all nodes except the first.
   */
  if (iam > 0) {
    setfcdd(ifst, main_auto_storage.d, fc, &nfpr, &nrc);
  }

#ifdef MATLAB_OUTPUT
  print_jacobian(*iap,main_auto_storage);
  {
    static num_calls = 0;
    char filename[80];
    sprintf(filename,"before%03d",num_calls);
    num_calls++;
    print_fa_fc(*iap,ft,fc,filename);
  }
#endif
  brbd(main_auto_storage.a, main_auto_storage.b, main_auto_storage.c, main_auto_storage.d, ft, fc, p0, p1, 
       ifst, &iid, nllv, &det, &ndim, &ntst, &nbc, &nrow, &nclm, &nfpr, &
       nrc, &iam, &kwt, &ipar, main_auto_storage.a1, main_auto_storage.a2, main_auto_storage.bb, 
       main_auto_storage.cc, main_auto_storage.faa, main_auto_storage.ca1, main_auto_storage.s1, main_auto_storage.s2, 
       main_auto_storage.icf11, main_auto_storage.ipr, main_auto_storage.icf1, main_auto_storage.icf2, 
       main_auto_storage.irf, main_auto_storage.icf);
#ifdef ACCES_TEST
    main_auto_storage.a = test.a;
    main_auto_storage.b = test.b;
    main_auto_storage.c = test.c;
#endif
  
  /*
    This is some stuff from the parallel version that isn't needed anymore 
    ----------------------------------------------------------------------
    lenft = ntst * nrow << 3;
    lenff = ntst0 * nrow << 3;
    jtmp1 = M 2AA;   I added spaces so these don't get flagged as header file macro dependancies
    jtmp2 = M 3AA;   I added spaces so these don't get flagged as header file macro dependancies
    lenff2 = jtmp1 * (jtmp2 + 1) << 3;
  */
  if (ipar) {
    /*        Global concatenation of the solution from each node. */
    integer tmp;
    gcol();
    tmp = iap->ntst + 1;
    faft(ff, fa, &ntst0, &nrow, ndxloc);
  } else {
    integer tmp;
    tmp = iap->ntst + 1;
    faft(ft, fa, &ntst0, &nrow, ndxloc);
  }
#ifdef MATLAB_OUTPUT
  {
    static num_calls = 0;
    char filename[80];
    sprintf(filename,"after%03d",num_calls);
    num_calls++;
    print_fa_fc(*iap,ft,fc,filename);
  }
#endif  

  rap->det = det;
  free(ff);
  free(ft);
  return 0;
} /* solvbv_ */


/*     ---------- ------- */
/* Subroutine */ int 
setfcdd(integer *ifst, doublereal *dd, doublereal *fc, integer *ncb, integer *nrc)
{
  /* System generated locals */
  integer dd_dim1;

    /* Local variables */
  integer i, j;



  /* Parameter adjustments */
  /*--fc;*/
  dd_dim1 = *ncb;
    
  for (i = 0; i < *nrc; ++i) {
    if (*ifst == 1) {
      for (j = 0; j < *ncb; ++j) {
	ARRAY2D(dd, j, i) = 0.;
      }
    }
    fc[i] = 0.;
  }


  return 0;
} /* setfcdd_ */


/*     ---------- ---- */
/* Subroutine */ int 
faft(doublereal *ff, doublereal *fa, integer *ntst, integer *nrow, integer *ndxloc)
{
  /* System generated locals */
  integer fa_dim1, ff_dim1;

    /* Local variables */
  integer i, j;



  /* Parameter adjustments */
  ff_dim1 = *nrow;
  fa_dim1 = *ndxloc;
    
  for (i = 0; i < *ntst; ++i) {
    for (j = 0; j < *nrow; ++j) {
      ARRAY2D(fa, i, j) = ARRAY2D(ff, j, i);
    }
  }

  return 0;
} /* faft_ */


/*     ---------- --------- */
/* Subroutine */ int 
partition(integer *n, integer *kwt, integer *m)
{
    /* Local variables */
  integer i, s, t;


  /*     Linear distribution of NTST over all nodes */

    /* Parameter adjustments */
    /*--m;*/

    
  t = *n / *kwt;
  s = *n % *kwt;

  for (i = 0; i < *kwt; ++i) {
    m[i] = t;
  }

  for (i = 0; i < s; ++i) {
    ++m[i];
  }

  return 0;
} /* partition_ */


/*     ------- -------- ------ */
integer 
mypart(integer *iam, integer *np)
{
  /* System generated locals */
  integer ret_val;

    /* Local variables */
  integer i, k;




  /*     Partition the mesh */


    /* Parameter adjustments */
    /*--np;*/

    
  k = 0;
  for (i = 0; i < *iam; ++i) {
    k += np[i];
  }
  ret_val = k;

  return ret_val;
} /* mypart_ */


/*     ---------- ------ */
/* Subroutine */ int 
setrhs(integer *ndim, integer *ips, integer *na, integer *ntst, integer *np, integer *ncol, integer *nbc, integer *nint, integer *ncb, integer *nrc, integer *nra, integer *nca, integer *iam, integer *kwt, logical *ipar, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), integer *ndxloc, iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *rds, doublereal *fa, doublereal *fc, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *ups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *dups, doublereal *dtm, doublereal *thl, doublereal *thu, doublereal *p0, doublereal *p1)
{
  /* System generated locals */
  integer ups_dim1, dups_dim1, uoldps_dim1, 
    udotps_dim1, upoldp_dim1, 
    fa_dim1, wt_dim1, wp_dim1, wploc_dim1;

  integer i, j, k, l, m;
  integer mpart, i1, j1, k1, l1;

  doublereal rlsum;
  integer ib, ic, jj;
  integer ic1;

  integer jp1;
  integer ncp1;
  doublereal dt,ddt;

  doublereal *dicd, *ficd, *dfdp, *dfdu, *uold;
  doublereal *f;
  doublereal *u, *wploc;
  doublereal *wi, *wp, *wt;
  doublereal *dbc, *fbc, *uic, *uio, *prm, *uid, *uip, *ubc0, *ubc1;

  dicd = (doublereal *)malloc(sizeof(doublereal)*(iap->nint)*(iap->ndim + NPARX));
  ficd = (doublereal *)malloc(sizeof(doublereal)*(iap->nint));
  dfdp = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*NPARX);
  dfdu = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*(iap->ndim));
  uold = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  f    = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  u    = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  wploc= (doublereal *)malloc(sizeof(doublereal)*(iap->ncol)*(iap->ncol+1));
  wi   = (doublereal *)malloc(sizeof(doublereal)*(iap->ncol+1) );
  wp   = (doublereal *)malloc(sizeof(doublereal)*(iap->ncol)*(iap->ncol+1) );
  wt   = (doublereal *)malloc(sizeof(doublereal)*(iap->ncol)*(iap->ncol+1) );
  dbc  = (doublereal *)malloc(sizeof(doublereal)*(iap->nbc)*(2*iap->ndim + NPARX));
  fbc  = (doublereal *)malloc(sizeof(doublereal)*(iap->nbc));
  uic  = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  uio  = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  prm  = (doublereal *)malloc(sizeof(doublereal)*NPARX);
  uid  = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  uip  = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  ubc0 = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  ubc1 = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));


  /* Parameter adjustments */
  /*--np;*/
  /*--par;*/
  /*--icp;*/
  /*--fc;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  /*--dtm;*/
  /*--thl;*/
  /*--thu;*/
  fa_dim1 = *nra;
  dups_dim1 = *ndxloc;
  upoldp_dim1 = *ndxloc;
  udotps_dim1 = *ndxloc;
  uoldps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
  wt_dim1 = iap->ncol+1;
  wp_dim1 = iap->ncol+1;
  wploc_dim1 = iap->ncol+1;
    
  *iam = iap->mynode;
  *kwt = iap->numnodes;
  if (*kwt > 1) {
    *ipar = TRUE_;
  } else {
    *ipar = FALSE_;
  }

  wint(*ncol + 1, wi);
  genwts(*ncol, iap->ncol + 1, wt, wp);
  /* Initialize to zero. */
  for (i = 0; i < *nrc; ++i) {
    fc[i] = 0.;
  }

  /* Set constants. */
  ncp1 = *ncol + 1;
  for (i = 0; i < *ncb; ++i) {
    par[icp[i]] = rlcur[i];
  }

  /* Generate FA : */

/*      Partition the mesh intervals. */
  mpart = mypart(iam, np);

  for (jj = 0; jj < *na; ++jj) {
    j = jj + mpart;
    jp1 = j + 1;
    dt = dtm[j];
    ddt = 1. / dt;
    for (ic = 0; ic < *ncol; ++ic) {
      for (ib = 0; ib < ncp1; ++ib) {
	ARRAY2D(wploc, ib, ic) = ddt * ARRAY2D(wp, ib, ic);
      }
    }
    for (ic = 0; ic < *ncol; ++ic) {
      for (k = 0; k < *ndim; ++k) {
	u[k] = ARRAY2D(wt, *ncol, ic) * ARRAY2D(ups, jp1, k);
	uold[k] = ARRAY2D(wt, *ncol, ic) * ARRAY2D(uoldps, jp1, k);
	for (l = 0; l < *ncol; ++l) {
	  l1 = l * *ndim + k;
	  u[k] += ARRAY2D(wt, l, ic) * ARRAY2D(ups, j, l1);
	  uold[k] += ARRAY2D(wt, l, ic) * ARRAY2D(uoldps, j, l1);
	}
      }
      /*     ** Time evolution computations (parabolic systems) */
      if (*ips == 14 || *ips == 16) {
	rap->tivp = rlold[0];
      }
      for (i = 0; i < NPARX; ++i) {
	prm[i] = par[i];
      }
      (*funi)(iap, rap, *ndim, u, uold, icp, prm, 2, f, 
	      dfdu, dfdp);
      ic1 = ic * *ndim;
      for (i = 0; i < *ndim; ++i) {
	ARRAY2D(fa, ic1 + i, jj) = f[i] - ARRAY2D(wploc, *ncol, ic) * ARRAY2D(ups, jp1, i);
	for (k = 0; k < *ncol; ++k) {
	  k1 = k * *ndim + i;
	  ARRAY2D(fa, ic1 + i, jj) -= ARRAY2D(wploc, k, ic) * ARRAY2D(ups, j, k1);
	}
      }
      /* L1: */
    }
    /* L2: */
  }

  /*     Generate FC : */

/*     Boundary conditions : */

  if (*nbc > 0) {
    for (i = 0; i < *ndim; ++i) {
      ubc0[i] = ARRAY2D(ups, 0, i);
      ubc1[i] = ARRAY2D(ups, *ntst, i);
    }
    (*bcni)(iap, rap, *ndim, par, icp, *nbc, ubc0, ubc1, 
	    fbc, 2, dbc);
    for (i = 0; i < *nbc; ++i) {
      fc[i] = -fbc[i];
    }
    /*       Save difference : */
    for (j = 0; j < *ntst + 1; ++j) {
      for (i = 0; i < *nra; ++i) {
	ARRAY2D(dups, j, i) = ARRAY2D(ups, j, i) - ARRAY2D(uoldps, j, i);
      }
    }
  }

  /*     Integral constraints : */
  if (*nint > 0) {
    for (jj = 0; jj < *na; ++jj) {
      j = jj + mpart;
      jp1 = j + 1;
      for (k = 0; k < ncp1; ++k) {
	for (i = 0; i < *ndim; ++i) {
	  i1 = k * *ndim + i;
	  j1 = j;
	  if (k + 1 == ncp1) {
	    i1 = i;
	  }
	  if (k + 1 == ncp1) {
	    j1 = jp1;
	  }
	  uic[i] = ARRAY2D(ups, j1, i1);
	  uio[i] = ARRAY2D(uoldps, j1, i1);
	  uid[i] = ARRAY2D(udotps, j1, i1);
	  uip[i] = ARRAY2D(upoldp, j1, i1);
	}
	(*icni)(iap, rap, *ndim, par, icp, *nint, uic, 
		uio, uid, uip, ficd, 2, dicd);
	for (m = 0; m < *nint; ++m) {
	  fc[*nbc + m] -= dtm[j] * wi[k] * ficd[m];
	}
      }
    }
  }

  /*     Pseudo-arclength equation : */
  rlsum = 0.;
  for (i = 0; i < *ncb; ++i) {
    rlsum += thl[icp[i]] * (rlcur[i] - rlold[i]) * rldot[i];
  }

  fc[-1 + *nrc] = *rds - rinpr(iap, ndim, ndxloc, udotps, 
			       dups, dtm, thu) - rlsum;

  free(dicd );
  free(ficd );
  free(dfdp );
  free(dfdu );
  free(uold );
  free(f    );
  free(u    );
  free(wploc);
  free(wi   );
  free(wp   );
  free(wt   );
  free(dbc  );
  free(fbc  );
  free(uic  );
  free(uio  );
  free(prm  );
  free(uid  );
  free(uip  );
  free(ubc0 );
  free(ubc1 );

  return 0;
} /* setrhs_ */


/*     ---------- ---- */
/* Subroutine */ int 
brbd(doublereal *a, doublereal *b, doublereal *c, doublereal *d, doublereal *fa, doublereal *fc, doublereal *p0, doublereal *p1, integer *ifst, integer *idb, integer *nllv, doublereal *det, integer *nov, integer *na, integer *nbc, integer *nra, integer *nca, integer *ncb, integer *nrc, integer *iam, integer *kwt, logical *par, doublereal *a1, doublereal *a2, doublereal *bb, doublereal *cc, doublereal *faa, doublereal *ca1, doublereal *s1, doublereal *s2, integer *icf11, integer *ipr, integer *icf1, integer *icf2, integer *irf, integer *icf)
{
  doublereal *e;
  doublereal *fcc;
  doublereal *sol1,*sol2,*sol3;

  e = (doublereal *)malloc(sizeof(doublereal)*(*nov + *nrc)*(*nov + *nrc));
  fcc = (doublereal *)malloc(sizeof(doublereal)*((*nov + *nrc) + (2*(*nov)*(*nov))+1));

  sol1 = (doublereal *)malloc(sizeof(doublereal)*(*nov)*(*na + 1));
  sol2 = (doublereal *)malloc(sizeof(doublereal)*(*nov)*(*na + 1));
  sol3 = (doublereal *)malloc(sizeof(doublereal)*(*nov)*(*na + 1));
  
  /* Local */

  /* Parameter adjustments */
  /*--icf;*/
  /*--irf;*/
  /*--icf2;*/
  /*--icf1;*/
  /*--ipr;*/
  /*--icf11;*/
  /*--s2;*/
  /*--s1;*/
  /*--ca1;*/
  /*--faa;*/
  /*--cc;*/
  /*--bb;*/
  /*--a2;*/
  /*--a1;*/
  /*--p1;*/
  /*--p0;*/
  /*--fc;*/
  /*--fa;*/
  /*--d;*/
  /*--c;*/
  /*--b;*/
  /*--a;*/

    
  if (*idb > 4 && *iam == 0) {
#ifndef ACCES_TEST
    print1(nov, na, nra, nca, ncb, nrc, a, b, c, d, &
    	   fa[0], fc);
#endif
  }
  if (*ifst == 1) {
#ifdef ACCES_TEST
    a = test.a;
    b = test.b;
    c = test.c;
#endif
    conpar(nov, na, nra, nca, a, ncb, b, nbc, nrc, c, d, irf, icf);
    copycp(iam, kwt, na, nov, nra, nca, a, ncb, b, nrc, c, 
	   a1, a2, bb, cc, irf);
  }

  if (*nllv == 0) {
    conrhs(nov, na, nra, nca, a, nbc, nrc, c, fa, fc, 
	   irf, icf, iam);
    cpyrhs(na, nov, nra, faa, fa, irf);
  } else {
#ifdef RANDY_FIX
    /* The faa array needs to be intialized as well, since it 
       it used in the dimrge_ rountine to print stuff out,
       and in the bcksub_ routine for actual computations! */
    {
      integer k;
      for(k=0;k<((*nov) * (*na + 1));k++)
	faa[k]=0.0;
    }
    setzero(fa, fc, na, nra, nrc);
#else
    setzero(fa, fc, na, nra, nrc);
    cpyrhs(na, nov, nra, faa, fa, irf);
#endif
  }

  if (*ifst == 1) {
    reduce(iam, kwt, par, a1, a2, bb, cc, d, na, 
	   nov, ncb, nrc, s1, s2, ca1, icf1, icf2, 
	   icf11, ipr, nbc);
  }

  if (*nllv == 0) {
    redrhs(iam, kwt, par, a1, a2, cc, faa, fc, na, 
	   nov, ncb, nrc, ca1, icf1, icf2, icf11, ipr,nbc);
  }

  dimrge(iam, kwt, par, e, cc, d, fc, ifst, na, 
	 nrc, nov, ncb, idb, nllv, fcc, p0, p1, det, s1, a2,
	 faa, bb);

  bcksub(iam, kwt, par, s1, s2, a2, bb, faa, fc, 
	 fcc, sol1, sol2, sol3, na, nov, ncb, icf2);

  infpar(iam, par, a, b, fa, sol1, sol2, fc, na, nov, nra, 
	 nca, ncb, irf, icf);

  free(e);
  free(fcc);
  free(sol1);
  free(sol2);
  free(sol3);
  return 0;
} /* brbd_ */


/*     ---------- ------- */
/* Subroutine */ int 
setzero(doublereal *fa, doublereal *fc, integer *na, integer *nra, integer *nrc)
{
  /* System generated locals */
  integer fa_dim1;

    /* Local variables */
  integer i, j;

    /* Parameter adjustments */
    /*--fc;*/
  fa_dim1 = *nra;
    
  for (i = 0; i < *na; ++i) {
    for (j = 0; j < *nra; ++j) {
      ARRAY2D(fa, j, i) = 0.;
    }
  }

  for (i = 0; i < *nrc; ++i) {
    fc[i] = 0.;
  }

  return 0;
} /* setzero_ */


/*     ---------- ------ */
/* Subroutine */ int 
conrhs(integer *nov, integer *na, integer *nra, integer *nca, doublereal *a, integer *nbc, integer *nrc, doublereal *c, doublereal *fa, doublereal *fc, integer *irf, integer *icf, integer *iam)
{
  /* System generated locals */
  integer icf_dim1, irf_dim1, a_dim1, a_dim2, c_dim1, c_dim2, fa_dim1;

    /* Local variables */
  integer nbcp1, i, icfic, irfir, m1, m2, ic, ir, irfirp, ir1, nex,
    irp;


    /* Parameter adjustments */
    /*--fc;*/
  irf_dim1 = *nra;
  fa_dim1 = *nra;
  icf_dim1 = *nca;
  a_dim1 = *nca;
  a_dim2 = *nra;
  c_dim1 = *nca;
  c_dim2 = *nrc;
    
  nex = *nca - (*nov * 2);
  if (nex == 0) {
    return 0;
  }

  /* Condensation of right hand side. */

  nbcp1 = *nbc + 1;
  m1 = *nov + 1;
  m2 = *nov + nex;

  for (i = 0; i < *na; ++i) {
    for (ic = *nov; ic < m2; ++ic) {
      ir1 = ic - *nov + 1;
      irp = ir1 - 1;
      irfirp = ARRAY2D(irf, irp, i);
      icfic = ARRAY2D(icf, ic, i);
      for (ir = ir1; ir < *nra; ++ir) {
	irfir = ARRAY2D(irf, ir, i);
	if (ARRAY3D(a, (icfic - 1), (irfir - 1), i) != (double)0.) {
	  ARRAY2D(fa, (irfir - 1), i) -= ARRAY3D(a, (icfic - 1), (irfir - 1), i) * ARRAY2D(fa, (irfirp - 1), i);
	}
      }
      for (ir = *nbc; ir < *nrc; ++ir) {
	if (ARRAY3D(c, (icfic - 1), ir, i) != (double)0.) {
	  fc[ir] -= ARRAY3D(c, (icfic - 1), ir, i) * ARRAY2D(fa, (irfirp - 1), i);
	}
      }
    }
  }

  return 0;
} /* conrhs_ */


/*     ---------- ------ */
/* Subroutine */ int 
copycp(integer *iam, integer *kwt, integer *na, integer *nov, integer *nra, integer *nca, doublereal *a, integer *ncb, doublereal *b, integer *nrc, doublereal *c, doublereal *a1, doublereal *a2, doublereal *bb, doublereal *cc, integer *irf)
{
  /* System generated locals */
  integer irf_dim1, a_dim1, a_dim2, b_dim1, b_dim2, c_dim1, c_dim2, a1_dim1, a1_dim2,
    a2_dim1, a2_dim2, bb_dim1, bb_dim2, cc_dim1,cc_dim2;

  /* Local variables */
  integer i, irfir, ic, ir, ic1, nap1;


/* Local */

/* Copies the condensed sytem generated by CONPAR into workspace. */

    /* Parameter adjustments */
  a2_dim1 = *nov;
  a2_dim2 = *nov;
  a1_dim1 = *nov;
  a1_dim2 = *nov;
  irf_dim1 = *nra;
  a_dim1 = *nca;
  a_dim2 = *nra;
  bb_dim1 = *nov;
  bb_dim2 = *ncb;
  b_dim1 = *ncb;
  b_dim2 = *nra;
  cc_dim1 = *nov;
  cc_dim2 = *nrc;
  c_dim1 = *nca;
  c_dim2 = *nrc;
    
  nap1 = *na + 1;
  for (i = 0; i < *na; ++i) {
    for (ir = 0; ir < *nov; ++ir) {
      irfir = ARRAY2D(irf, *nra - *nov + ir, i);
      for (ic = 0; ic < *nov; ++ic) {
	ic1 = *nca - *nov + ic;
	ARRAY3D(a1,  ir, ic, i) = ARRAY3D(a, ic, (irfir - 1), i);
	ARRAY3D(a2, ir, ic, i) = ARRAY3D(a, ic1, (irfir - 1), i);
      }
      for (ic = 0; ic < *ncb; ++ic) {
	ARRAY3D(bb, ir, ic, i) = ARRAY3D(b, ic, (irfir - 1), i);
      }
    }
  }

  for (i = 0; i < nap1; ++i) {
    for (ir = 0; ir < *nrc; ++ir) {
      for (ic = 0; ic < *nov; ++ic) {
	if (i + 1 == 1) {
	  ARRAY3D(cc, ic, ir, i) = ARRAY3D(c, ic, ir, i);
	} else if (i + 1 == nap1) {
	  ARRAY3D(cc, ic, ir, i) = ARRAY3D(c, *nra + ic, ir, (i - 1));
	} else {
	  ARRAY3D(cc, ic, ir, i) = ARRAY3D(c, ic, ir, i) + ARRAY3D(c, *nra + ic, ir, (i - 1));
	}
      }
    }
  }

  return 0;
} /* copycp_ */


/*     ---------- ------ */
/* Subroutine */ int 
cpyrhs(integer *na, integer *nov, integer *nra, doublereal *faa, doublereal *fa, integer *irf)
{
  /* System generated locals */
  integer irf_dim1, fa_dim1, faa_dim1;

  /* Local variables */
  integer i, irfir, ir;

/*     **Copy the RHS */
    /* Parameter adjustments */
  faa_dim1 = *nov;
  irf_dim1 = *nra;
  fa_dim1 = *nra;

    
  for (i = 0; i < *na; ++i) {
    for (ir = 0; ir < *nov; ++ir) {
      irfir = ARRAY2D(irf, *nra - *nov + ir, i);
      ARRAY2D(faa, ir, i) = ARRAY2D(fa, (irfir - 1), i);
    }
  }

  return 0;
} /* cpyrhs_ */


/*     ---------- ------ */
/* Subroutine */ int 
reduce(integer *iam, integer *kwt, logical *par, doublereal *a1, doublereal *a2, doublereal *bb, doublereal *cc, doublereal *dd, integer *na, integer *nov, integer *ncb, integer *nrc, doublereal *s1, doublereal *s2, doublereal *ca1, integer *icf1, integer *icf2, integer *icf11, integer *ipr, integer *nbc)
{
  /* System generated locals */
  integer icf1_dim1, icf2_dim1, icf11_dim1, a1_dim1, a1_dim2, a2_dim1, a2_dim2, 
    s1_dim1, s1_dim2, s2_dim1, s2_dim2, 
    bb_dim1, bb_dim2, cc_dim1, cc_dim2, 
    dd_dim1, ca1_dim1, ca1_dim2,
    ipr_dim1;

    /* Local variables */
  logical oddc[KREDO];
  integer niam, ibuf, ismc[KREDO], irmc[KREDO], info, irmm[KREDO], 
    ismm[KREDO], nlev, itmp;
  doublereal zero, tpiv;
  real xkwt;
  integer nbcp1, ibuf1, ipiv1, jpiv1, ipiv2, jpiv2, i, k, l;

  logical evenc[KREDO];

  integer i1, i2, k1, k2, i3, l1, iprow, k3, l2, l3, ic, ir;
  doublereal rm;
  logical master[KREDO];
  integer ib1, ib2, myleft[KREDO];

  logical worker[KREDO];
  integer ir1, iprown, iprown2, ism[KREDO], irm[KREDO], nrcmnbc;
  doublereal tmp;
  integer myleftc[KREDO];
  logical notsend;
  integer nap1, myright[KREDO], nam1, len1, len2, icp1;
  doublereal piv1, piv2;
  doublereal *buf=NULL;


  /* Parameter adjustments */
  ipr_dim1 = *nov;
  icf11_dim1 = *nov;
  icf2_dim1 = *nov;
  icf1_dim1 = *nov;
  ca1_dim1 = *nov;
  ca1_dim2 = *nov;
  s2_dim1 = *nov;
  s2_dim2 = *nov;
  s1_dim1 = *nov;
  s1_dim2 = *nov;
  a2_dim1 = *nov;
  a2_dim2 = *nov;
  a1_dim1 = *nov;
  a1_dim2 = *nov;
  dd_dim1 = *ncb;
  bb_dim1 = *nov;
  bb_dim2 = *ncb;
  cc_dim1 = *nov;
  cc_dim2 = *nrc;
    
  zero = 0.;
  nbcp1 = *nbc + 1;
  nap1 = *na + 1;
  nam1 = *na - 1;
  nrcmnbc = *nrc - *nbc;
  len1 = (*nov * (*nrc - *nbc)) * 8;
  len2 = (*nov + *nrc - *nbc + 1) * 8;
  xkwt = (real) (*kwt);
  {
    real tmp = r_lg10(xkwt) / r_lg10(2.0);
    nlev = i_nint(&tmp);
  }
  notsend = TRUE_;

/*     FOR EACH REURSIVE LEVEL, CALCULATE THE MASTER(HOLDING THE */
/*     PIVOT ROW AFTER ROW SWAPPING) NODE WHICH WILL SEND THE */
/*     PIVOT ROW TO THE CORRESPONDING WORKER NODE WHICH IS DISTANCED */
/*     2**(K-1) FROM THE MASTER WHERE K IS THE RECURSIVE LEVEL NUMBER. */
/*     THE CORRESPONDING MESSAGE TYPE IN EACH RECURSIVE LEVEL IS */
/*     ALSO CALCULATED HERE. */

/* For each level in the recursion, determine the master node */
/* (holding the pivot row after row swapping), which will send the */
/* pivot row to the corresponding worker node at distance 2**(K-1) */
/* from the master. Here K is the level in the recursion. */
/* The message type at each level in the recursion is also determined. */

  if (*par) {

    for (i = 0; i < nlev; ++i) {

      oddc[i] = FALSE_;
      evenc[i] = FALSE_;
      master[i] = FALSE_;
      worker[i] = FALSE_;

      k1 = pow_ii(2, i);
      k2 = k1 * 2;
      niam = *iam / k1;

      if (notsend) {

	if (niam % 2 == 0) {

	  master[i] = TRUE_;
	  notsend = FALSE_;
	  ism[i] = (i + 1) + *iam;
	  irm[i] = ism[i] + k1;
	  myright[i] = *iam + k1;
	  irmm[i] = (i + 1) + *iam + 1 + (*kwt * 2);
	  ismc[i] = (i + 1) + *iam + *kwt;
	  myleftc[i] = *iam - (k1 - 1);

	} else {

	  worker[i] = TRUE_;
	  ism[i] = (i + 1) + *iam;
	  irm[i] = ism[i] - k1;
	  myleft[i] = *iam - k1;

	}

      }

      k = *iam % k2;
      if (k == k1) {
	evenc[i] = TRUE_;
	ismm[i] = (i + 1) + *iam + (*kwt * 2);
      }

      if (*iam % k2 == 0) {
	oddc[i] = TRUE_;
	irmc[i] = (i + 1) + *iam + *kwt + (k1 - 1);
      }

      /* L1: */
    }
  }

  /* Initialization */

  for (i = 0; i < *na; ++i) {
    for (k1 = 0; k1 < *nov; ++k1) {
      ARRAY2D(icf1, k1, i) = k1 + 1;
      ARRAY2D(icf2, k1, i) = k1 + 1;
      ARRAY2D(ipr, k1, i) = k1 + 1;
      for (k2 = 0; k2 < *nov; ++k2) {
	ARRAY3D(s2, k1, k2, i) = 0.;
	ARRAY3D(s1, k1, k2, i) = 0.;
      }
    }
  }

  for (ir = 0; ir < *nov; ++ir) {
    for (ic = 0; ic < *nov; ++ic) {
      ARRAY2D(s1, ir, ic) = ARRAY2D(a1, ir, ic);
    }
  }

  /* The reduction process is done concurrently */
  for (i1 = 0; i1 < nam1; ++i1) {

    i2 = i1 + 1;
    i3 = i2 + 1;

    for (ic = 0; ic < *nov; ++ic) {
      icp1 = ic + 1;

      /* Complete pivoting; rows are swapped physically, columns swap in
	 dices */
      piv1 = zero;
      ipiv1 = ic + 1;
      jpiv1 = ic + 1;
      for (k1 = ic; k1 < *nov; ++k1) {
	for (k2 = ic; k2 < *nov; ++k2) {
	  tpiv = ARRAY3D(a2, k1, ARRAY2D(icf2, k2, i1) - 1, i1);
	  if (tpiv < zero) {
	    tpiv = -tpiv;
	  }
	  if (piv1 < tpiv) {
	    piv1 = tpiv;
	    ipiv1 = k1 + 1;
	    jpiv1 = k2 + 1;
	  }
	}
      }

      piv2 = zero;
      ipiv2 = 1;
      jpiv2 = ic + 1;
      for (k1 = 0; k1 < *nov; ++k1) {
	for (k2 = ic; k2 < *nov; ++k2) {
	  tpiv = ARRAY3D(a1, k1, ARRAY2D(icf1, k2, i2) - 1, i2);
	  if (tpiv < zero) {
	    tpiv = -tpiv;
	  }
	  if (piv2 < tpiv) {
	    piv2 = tpiv;
	    ipiv2 = k1 + 1;
	    jpiv2 = k2 + 1;
	  }
	}
      }

      if (piv1 >= piv2) {
	ARRAY2D(ipr, ic, i1) = ipiv1;
	itmp = ARRAY2D(icf2, ic, i1);
	ARRAY2D(icf2, ic, i1) = ARRAY2D(icf2, (jpiv1 - 1), i1);
	ARRAY2D(icf2, (jpiv1 - 1), i1) = itmp;
	itmp = ARRAY2D(icf1, ic, i2);
	ARRAY2D(icf1, ic, i2) = ARRAY2D(icf1, (jpiv1 - 1), i2);
	ARRAY2D(icf1, (jpiv1 - 1), i2) = itmp;
	/* Swapping */
	for (l = 0; l < *nov; ++l) {
	  tmp = ARRAY3D(s1, ic, l, i1);
	  ARRAY3D(s1, ic, l, i1) = ARRAY3D(s1, (ipiv1 - 1), l, i1);
	  ARRAY3D(s1, (ipiv1 - 1), l, i1) = tmp;
	  if (l >= ic) {
	    tmp = ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
	    ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1) = 
	      ARRAY3D(a2, (ipiv1 - 1), ARRAY2D(icf2, l, i1) - 1, i1);
	    ARRAY3D(a2, (ipiv1 - 1), ARRAY2D(icf2, l, i1) - 1, i1) = tmp;
	  }
	  tmp = ARRAY3D(s2, ic, l, i1);
	  ARRAY3D(s2, ic, l, i1) = ARRAY3D(s2, (ipiv1 - 1), l, i1);
	  ARRAY3D(s2, (ipiv1 - 1), l, i1) = tmp;
	}

	for (l = 0; l < *ncb; ++l) {
	  tmp = ARRAY3D(bb, ic, l, i1);
	  ARRAY3D(bb, ic, l, i1) = ARRAY3D(bb, (ipiv1 - 1), l, i1);
	  ARRAY3D(bb, (ipiv1 - 1), l, i1) = tmp;
	}
      } else {
	ARRAY2D(ipr, ic, i1) = *nov + ipiv2;
	itmp = ARRAY2D(icf2, ic, i1);
	ARRAY2D(icf2, ic, i1) = ARRAY2D(icf2, (jpiv2 - 1), i1);
	ARRAY2D(icf2, (jpiv2 - 1), i1) = itmp;
	itmp = ARRAY2D(icf1, ic, i2);
	ARRAY2D(icf1, ic, i2) = ARRAY2D(icf1, (jpiv2 - 1), i2);
	ARRAY2D(icf1, (jpiv2 - 1), i2) = itmp;
	/* Swapping */
	for (l = 0; l < *nov; ++l) {
	  if (l >= ic) {
	    tmp = ARRAY3D(a2, ic , ARRAY2D(icf2, l, i1) - 1, i1);
	    ARRAY3D(a2, ic , ARRAY2D(icf2, l, i1) - 1, i1) = 
	      ARRAY3D(a1, (ipiv2 - 1), ARRAY2D(icf2, l, i1) - 1, i2);
	    ARRAY3D(a1, (ipiv2 - 1), ARRAY2D(icf2, l, i1) - 1, i2) = tmp;
	  }
	  tmp = ARRAY3D(s2, ic, l, i1);
	  ARRAY3D(s2, ic, l, i1) = ARRAY3D(a2, (ipiv2 - 1), l, i2);
	  ARRAY3D(a2, (ipiv2 - 1), l, i2) = tmp;
	  tmp = ARRAY3D(s1, ic, l, i1);
	  ARRAY3D(s1, ic, l, i1) = ARRAY3D(s1, (ipiv2 - 1), l, i2);
	  ARRAY3D(s1, (ipiv2 - 1), l, i2) = tmp;
	}
	for (l = 0; l < *ncb; ++l) {
	  tmp = ARRAY3D(bb, ic, l, i1);
	  ARRAY3D(bb, ic, l, i1) = ARRAY3D(bb, (ipiv2 - 1), l, i2);
	  ARRAY3D(bb, (ipiv2 - 1), l, i2) = tmp;
	}
      }

      /* End of pivoting; Elimination starts here */

      for (ir = icp1; ir < *nov; ++ir) {
	rm = ARRAY3D(a2, ir, ARRAY2D(icf2, ic, i1) - 1, i1) / 
	  ARRAY3D(a2, ic, ARRAY2D(icf2, ic, i1) - 1, i1);
	ARRAY3D(a2, ir, ARRAY2D(icf2, ic, i1) - 1, i1) = rm;

	if (rm != (double)0.) {
	  for (l = icp1; l < *nov; ++l) {
	    ARRAY3D(a2, ir, ARRAY2D(icf2, l, i1) - 1, i1) -= 
	      rm * ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
	  }

	  for (l = 0; l < *nov; ++l) {
	    ARRAY3D(s1, ir, l, i1) -= rm * ARRAY3D(s1, ic, l, i1);
	    ARRAY3D(s2, ir, l, i1) -= rm * ARRAY3D(s2, ic, l, i1);
	  }

	  for (l = 0; l < *ncb; ++l) {
	    ARRAY3D(bb, ir, l, i1) -= rm * ARRAY3D(bb, ic, l, i1);
	  }
	}
      }

      for (ir = 0; ir < *nov; ++ir) {
	rm = ARRAY3D(a1, ir, ARRAY2D(icf1, ic, i2) - 1, i2) / 
	  ARRAY3D(a2, ic, ARRAY2D(icf2, ic, i1) - 1, i1);
	ARRAY3D(a1, ir, ARRAY2D(icf1, ic, i2) - 1, i2) = rm;

	if (rm != (double)0.) {
	  for (l = icp1; l < *nov; ++l) {
	    ARRAY3D(a1, ir, ARRAY2D(icf1, l, i2) - 1, i2) -= 
	      rm * ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
	  }
	  for (l = 0; l < *nov; ++l) {
	    ARRAY3D(s1, ir, l, i2) -= rm * ARRAY3D(s1, ic, l, i1);
	    ARRAY3D(a2, ir, l, i2) -= rm * ARRAY3D(s2, ic, l, i1);
	  }
	  for (l = 0; l < *ncb; ++l) {
	    ARRAY3D(bb, ir, l, i2) -= rm * ARRAY3D(bb, ic, l, i1);
	  }
	}
      }

      for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
	rm = ARRAY3D(cc, ARRAY2D(icf2, ic, i1) - 1, ir, i2) / 
	  ARRAY3D(a2, ic, ARRAY2D(icf2, ic, i1) - 1, i1);
	ARRAY3D(cc, ARRAY2D(icf2, ic, i1) - 1, ir, i2) = rm;

	if (rm != (double)0.) {
	  for (l = icp1; l < *nov; ++l) {
	    ARRAY3D(cc, ARRAY2D(icf2, l, i1) - 1, ir, i2) -= 
	      rm * ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
	  }
	  for (l = 0; l < *nov; ++l) {
	    ARRAY3D(cc, l, ir, 0) -= rm * ARRAY3D(s1, ic, l, i1);
	    ARRAY3D(cc, l, ir, i3) -= rm * ARRAY3D(s2, ic, l, i1);
	  }
	  for (l = 0; l < *ncb; ++l) {
	    ARRAY2D(dd, l, ir) -= rm * ARRAY3D(bb, ic, l, i1);
	  }
	}
      }

      /* L2: */
    }
    /* L3: */
  }

  /* Initialization */
  for (i = 0; i < *nov; ++i) {
    ARRAY2D(icf2, i, (*na - 1)) = i + 1;
  }

  /*     INTER NODES REDUCE IS DONE VIA COMMUNICATION */
  /*     BETWEEN MASTER NODES AND WORKER NODES. */
  /*     THE SUMMATION OF THE OVERLAPED PART C IN THE */
  /*     NEIGHBOR NODES IN THE CONDENSATION OF PARAMETE */
  /*     ROUTINE IS DELAYED TO HERE TO SUM. */

/* Inter node reduction is done via communication between master node */
/* and worker nodes. The summation over the overlapped part C of */
/*neighboring nodes in the condensation of parameters is delayed until her
e.*/
  if (*par) {

    for (i = 0; i < nlev; ++i) {

      if (master[i]) {
	crecv();
	for (ir = nbcp1; ir <= *nrc; ++ir) {
	  ir1 = ir - *nbc;
	  for (ic = 1; ic <= *nov; ++ic) {
	    l1 = ir1 * *nov + ic;
	    ARRAY3D(cc, ic, ir, (nap1 - 1)) += buf[l1 + 1];
	  }
	}
	for (ir = 0; ir < *nov; ++ir) {
	  for (ic = 0; ic < *nov; ++ic) {
	    ARRAY3D(s2, ir, ic, *na) = 0.;
	  }
	}
      }

      if (evenc[i]) {
	csend();
      }

      if (worker[i]) {
	for (ir = 0; ir < *nov; ++ir) {
	  for (ic = 0; ic < *nov; ++ic) {
	    ARRAY3D(ca1, ir, ic, i) = ARRAY3D(s1, ir, ic, (*na - 1));
	    ARRAY3D(s1, ir, ic, (*na - 1)) = 0.;
	  }
	}

	for (l = 0; l < *nov; ++l) {
	  ARRAY2D(icf11, l, i) = l + 1;
	}
      }

      for (ic = 0; ic < *nov; ++ic) {

	icp1 = ic + 1;
	iprow = *nov - ic + 1;
	iprown = iprow + *nov;
	iprown2 = iprown + *nov;
	ib1 = iprown2 + *ncb + 1;
	ib2 = ib1 + 1;
	ibuf = (ib2 + 1) * 8;
	ibuf1 = (ib2 + *nrc - *nbc) * 8;

	if (master[i]) {

	  /* PIVOTING (COMPLETE PIVOTING) */

	  piv1 = zero;
	  ipiv1 = ic + 1;
	  jpiv1 = ic + 1;
	  for (k1 = ic; k1 < *nov; ++k1) {
	    for (k2 = ic; k2 < *nov; ++k2) {
	      k3 = ARRAY2D(icf2, k2, (*na - 1));
	      tpiv = ARRAY3D(a2, k1, k3, (*na - 1));
	      if (tpiv < zero) {
		tpiv = -tpiv;
	      }
	      if (piv1 < tpiv) {
		piv1 = tpiv;
		ipiv1 = k1 + 1;
		jpiv1 = k2 + 1;
	      }
	    }
	  }

	  crecv();

	  jpiv2 = i_dnnt(&buf[ib1 + 1]);
	  ipiv2 = i_dnnt(&buf[ib2 + 1]);

	  piv2 = buf[1];
	  if (piv2 < 0.) {
	    piv2 = -piv2;
	  }

	  if (piv1 >= piv2) {

	    ARRAY2D(ipr, ic, (*na - 1)) = ipiv1;
	    itmp = ARRAY2D(icf2, ic, (*na - 1));
	    ARRAY2D(icf2, ic, (*na - 1)) = ARRAY2D(icf2, (jpiv1 - 1), (*na - 1));
	    ARRAY2D(icf2, (jpiv1 - 1), (*na - 1)) = itmp;

	    /* Send pivot row to worker */
	    for (l = 0; l < *nov; ++l) {
	      if (l >= ic) {
		l1 = l - ic + 2;
		l2 = ARRAY2D(icf2, l, (*na - 1)) - 1;
		buf[l1 + 1] = ARRAY3D(a2, (ipiv1 - 1), l2, (*na - 1));
	      }
	      l1 = iprow + l;
	      l2 = iprown + l;
	      buf[l1 + 1] = ARRAY3D(s1, (ipiv1 - 1), l, (*na - 1));
	      buf[l2 + 1] = ARRAY3D(s2, (ipiv1 - 1), l, (*na - 1));
	    }

	    for (l = 0; l < *nbc; ++l) {
	      l1 = iprown2 + l;
	      buf[l1 + 1] = ARRAY3D(bb, (ipiv1 - 1), l, (*na - 1));
	    }

	    buf[ib1 + 1] = (doublereal) jpiv1;

	    for (l = nbcp1 - 1; l < *nrc; ++l) {
	      l1 = l - *nbc;
	      l2 = ib1 + l1;
	      l3 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
	      buf[l2 + 1] = ARRAY3D(cc, l3, l, (nap1 - 1));
	    }

	    l1 = ib2 + nrcmnbc;
	    buf[l1 + 1] = 0.;
	    csend();

	    /* Row swapping */
	    for (l = 0; l < *nov; ++l) {
	      tmp = ARRAY3D(s1, ic, l, (*na - 1));
	      ARRAY3D(s1, ic, l, (*na - 1)) = ARRAY3D(s1, (ipiv1 - 1), l, (*na - 1));
	      ARRAY3D(s1, (ipiv1 - 1), l, (*na - 1)) = tmp;
	      if (l >= ic) {
		l1 = ARRAY2D(icf2, l, (*na - 1)) - 1;
		tmp = ARRAY3D(a2, ic, l1, (*na - 1));
		ARRAY3D(a2, ic, l1, (*na - 1)) = ARRAY3D(a2, (ipiv1 - 1), l1, (*na - 1));
		ARRAY3D(a2, (ipiv1 - 1), l1, (*na - 1)) = tmp;
	      }
	      tmp = ARRAY3D(s2, ic, l, (*na - 1));
	      ARRAY3D(s2, ic, l, (*na - 1)) = ARRAY3D(s2, (ipiv1 - 1), l, (*na - 1));
	      ARRAY3D(s2, (ipiv1 - 1), l, (*na - 1)) = tmp;
	    }

	    for (l = 0; l < *ncb; ++l) {
	      tmp = ARRAY3D(bb, ic, l, (*na - 1));
	      ARRAY3D(bb, ic, l, (*na - 1)) = ARRAY3D(bb, (ipiv1 - 1), l, (*na - 1));
	      ARRAY3D(bb, (ipiv1 - 1), l, (*na - 1)) = tmp;
	    }

	  } else {

	    ARRAY2D(ipr, ic, (*na - 1)) = *nov + ipiv2;
	    jpiv1 = jpiv2;
	    itmp = ARRAY2D(icf2, ic, (*na - 1));
	    ARRAY2D(icf2, ic, (*na - 1)) = ARRAY2D(icf2, (jpiv1 - 1), (*na - 1));
	    ARRAY2D(icf2, (jpiv1 - 1), (*na - 1)) = itmp;

	    for (l = 0; l < *nov; ++l) {
	      if (l >= ic) {
		l1 = l - ic + 2;
		l2 = ARRAY2D(icf2, l, (*na - 1)) - 1;
		tmp = buf[l1 + 1];
		buf[l1 + 1] = ARRAY3D(a2, ic, l2, (*na - 1));
		ARRAY3D(a2, ic, l2, (*na - 1)) = tmp;
	      }
	      l1 = iprow + l;
	      l2 = iprown + l;
	      tmp = buf[l1 + 1];
	      buf[l1 + 1] = ARRAY3D(s1, ic, l, (*na - 1));
	      ARRAY3D(s1, ic, l, (*na - 1)) = tmp;
	      tmp = buf[l2 + 1];
	      buf[l2 + 1] = ARRAY3D(s2, ic, l, (*na - 1));
	      ARRAY3D(s2, ic, l, (*na - 1)) = tmp;
	    }

	    for (l = 0; l < *nbc; ++l) {
	      l1 = iprown2 + l;
	      tmp = buf[l1 + 1];
	      buf[l1 + 1] = ARRAY3D(bb, ic, l, (*na - 1));
	      ARRAY3D(bb, ic, l, (*na - 1)) = tmp;
	    }

	    buf[ib1 + 1] = (doublereal) jpiv2;

	    for (l = nbcp1; l <= *nrc; ++l) {
	      l1 = l - *nbc;
	      l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
	      l3 = ib1 + l1  + 1;
	      buf[l3 + 1] = ARRAY3D(cc, l2, l, (nap1 - 1 ));
	    }
	    l1 = ib2 + nrcmnbc;
	    buf[l1 + 1] = 1.;

	    csend();


	  }
	  /* End pivoting in master */

	  /* Send data to worker nodes */
	  for (l = 0; l < *nov; ++l) {
	    buf[l + 1] = ARRAY3D(s1, ic, l, (*na - 1));
	  }

	  for (l = nbcp1 - 1; l <*nrc; ++l) {
	    l1 = l - *nbc;
	    l2 = ARRAY2D(icf2, ic, (*na -1 )) - 1;
	    l3 = *nov + l1;
	    buf[l3 + 1] = ARRAY3D(cc, l2, l, (nap1 - 1));
	  }

	  l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
	  l1 = *nov + nrcmnbc;
	  buf[l1 + 1] = ARRAY3D(a2, ic, l2, (*na - 1));

	  csend();

	  /* Elimination */
	  for (ir = icp1; ir < *nov; ++ir) {
	    l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
	    rm = ARRAY3D(a2, ir, l2, (*na - 1)) / ARRAY3D(a2, ic, l2, (*na - 1));
	    ARRAY3D(a2, ir, l2, (*na - 1)) = rm;
	    if (rm != zero) {
	      for (l = icp1; l < *nov; ++l) {
		l1 = ARRAY2D(icf2, l, (*na - 1)) - 1;
		ARRAY3D(a2, ir, l1, (*na - 1)) -= rm * ARRAY3D(a2, ic, l1, (*na - 1));
	      }
	      for (l = 0; l < *nov; ++l) {
		ARRAY3D(s1, ir, l, (*na - 1)) -= rm *ARRAY3D(s1, ic, l, (*na - 1));
		ARRAY3D(s2, ir, l, (*na - 1)) -= rm *ARRAY3D(s2, ic, l, (*na - 1));
	      }
	      for (l = 0; l < *ncb; ++l) {
		ARRAY3D(bb, ir, l, (*na - 1)) -= rm *ARRAY3D(bb, ic, l, (*na - 1));
	      }
	    }
	  }

	  for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
	    l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
	    rm = ARRAY3D(cc, l2, ir, (nap1 - 1)) / ARRAY3D(a2, ic, l2, (*na - 1));
	    ARRAY3D(cc, l2, ir, (nap1 - 1)) = rm;
	    if (rm != zero) {
	      for (l = icp1; l <= *nov; ++l) {
		l1 = ARRAY2D(icf2, l, (*na - 1)) - 1;
		ARRAY3D(cc, l1, ir, (nap1 - 1)) -= rm * ARRAY3D(a2, ic, l1, (*na - 1));
	      }
	      for (l = 0; l < *nbc; ++l) {
		ARRAY2D(dd, l, ir) -= rm * ARRAY3D(bb, ic, l, (*na - 1));
	      }
	    }
	  }

	}

	if (worker[i]) {

	  /* Pivoting */
	  piv2 = zero;
	  ipiv2 = 1;
	  jpiv2 = ic + 1;
	  for (k1 = 0; k1 < *nov; ++k1) {
	    for (k2 = ic; k2 < *nov; ++k2) {
	      k3 = ARRAY2D(icf11, k2, i) - 1;
	      tpiv = ARRAY3D(ca1, k1, k3, i);
	      if (tpiv < zero) {
		tpiv = -tpiv;
	      }
	      if (piv2 < tpiv) {
		piv2 = tpiv;
		ipiv2 = k1 + 1;
		jpiv2 = k2 + 1;
	      }
	    }
	  }

	  itmp = ARRAY2D(icf11, ic, i);
	  ARRAY2D(icf11, ic, i) = ARRAY2D(icf11, (jpiv2 - 1), i);
	  ARRAY2D(icf11, (jpiv2 - 1), i) = itmp;

	  for (l = 0; l < *nov; ++l) {
	    if (l >= ic) {
	      l1 = l - ic + 2;
	      l2 = ARRAY2D(icf11, l, i) - 1;
	      buf[l1 + 1] = ARRAY3D(ca1, (ipiv2 - 1), l2, (*na - 1));
	    }
	    l1 = iprow + l;
	    l2 = l1 + *nov;
	    buf[l1 + 1] = ARRAY3D(s1, (ipiv2 - 1), l, i);
	    buf[l2 + 1] = ARRAY3D(a2, (ipiv2 - 1), l, (*na - 1));
	  }

	  for (l = 0; l < *ncb; ++l) {
	    l1 = iprown2 + l;
	    buf[l1 + 1] = ARRAY3D(bb, (ipiv2 - 1), l, (*na - 1));
	  }

	  buf[ib1 + 1] = (doublereal) jpiv2;
	  buf[ib2 + 1] = (doublereal) ipiv2;

	  csend();
	  crecv();

	  l1 = ib2 + nrcmnbc;
	  info = i_dnnt(&buf[l1 + 1]);

	  if (info == 1) {
	    /* Send pivot row to master */
	    for (l = 0; l < *nov; ++l) {
	      if (l >= ic) {
		l1 = l - ic + 2;
		l2 = ARRAY2D(icf11, l, i) - 1;
		tmp = ARRAY3D(ca1, (ipiv2 - 1), l2, i);
		ARRAY3D(ca1, (ipiv2 - 1), l2, i) = buf[l1 + 1];
		buf[l1 + 1] = tmp;
	      }
	      l1 = iprow + l;
	      l2 = l1 + *nov;
	      tmp = ARRAY3D(s1, (ipiv2 - 1), l, (*na - 1));
	      ARRAY3D(s1, (ipiv2 - 1), l, (*na - 1)) = buf[l1 + 1];
	      buf[l1 + 1] = tmp;
	      tmp = ARRAY3D(a2, (ipiv2 - 1), l, (*na - 1));
	      ARRAY3D(a2, (ipiv2 - 1), l, (*na - 1)) = buf[l2 + 1];
	      buf[l2 + 1] = tmp;
	    }
	    for (l = 0; l < *nbc; ++l) {
	      l1 = iprown2 + l;
	      tmp = ARRAY3D(bb, (ipiv2 - 1), l, (*na - 1));
	      ARRAY3D(bb, (ipiv2 - 1), l, (*na - 1)) = buf[l1 + 1];
	      buf[l1 + 1] = tmp;
	    }
	  } else {

	    itmp = ARRAY2D(icf11, ic, i);
	    ARRAY2D(icf11, ic, i) = ARRAY2D(icf11, (jpiv2 - 1), i);
	    ARRAY2D(icf11, (jpiv2 - 1), i) = itmp;

	    jpiv2 = i_dnnt(&buf[ib1 + 1]);
	    itmp = ARRAY2D(icf11, ic, i);
	    ARRAY2D(icf11, ic, i) = ARRAY2D(icf11, (jpiv2 - 1), i);
	    ARRAY2D(icf11, (jpiv2 - 1), i) = itmp;
	  }

	  /* Elimination */
	  for (ir = 1; ir <= *nov; ++ir) {
	    l2 = ARRAY2D(icf11, ic, i) - 1;
	    rm = ARRAY3D(ca1, ir, l2, i) / buf[1];
	    ARRAY3D(ca1, ir, l2, i) = rm;

	    if (rm != zero) {
	      for (l = icp1; l < *nov; ++l) {
		l1 = l - icp1 + 3;
		l3 = ARRAY2D(icf11, l, i) - 1;
		ARRAY3D(ca1, ir, l3, i) -= rm * buf[l1 + 1];
	      }
	      for (l = 0; l < *nov; ++l) {
		l1 = iprow + l;
		l2 = l1 + *nov;
		ARRAY3D(s1, ir, l, (*na - 1)) -= rm *buf[l1 + 1];
		ARRAY3D(a2, ir, l, (*na - 1)) -= rm *buf[l2 + 1];
	      }
	      for (l = 0; l < *ncb; ++l) {
		l1 = iprown2 + l;
		ARRAY3D(bb, ir, l, (*na - 1)) -= rm *buf[l1 + 1];
	      }
	    }
	  }

	  for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
	    l1 = ir - *nbc;
	    l2 = ib1 + l1;
	    rm = buf[l2 + 1] / buf[1];
	    if (rm != zero) {
	      for (l = 0; l < *nov; ++l) {
		l3 = iprown + l;
		ARRAY3D(cc, l, ir, (nap1 - 1)) -= rm * buf[l3 + 1];
	      }
	    }

	  }

	}

	if (oddc[i]) {
	  crecv();
	  for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
	    ir1 = ir - *nbc;
	    l1 = *nov + nrcmnbc;
	    l2 = *nov + ir1;
	    rm = buf[l2 + 1] / buf[l1 + 1];
	    if (rm != zero) {
	      for (l = 0; l < *nov; ++l) {
		ARRAY3D(cc, l, ir, 0) -= rm * buf[l + 1];
	      }
	    }
	  }
	}

      }

      /* L4: */
    }

    /* Global sum for D by recursive doubling */
    {
      integer tmp= (*nrc - *nbc) * *ncb;
      rd0(iam, kwt, &ARRAY2D(dd, 0, (nbcp1 - 1)), &tmp);
    }
	
  }
    
  return 0;
} /* reduce_ */


/*     ---------- ------ */
/* Subroutine */ int 
redrhs(integer *iam, integer *kwt, logical *par, doublereal *a1, doublereal *a2, doublereal *cc, doublereal *faa, doublereal *fc, integer *na, integer *nov, integer *ncb, integer *nrc, doublereal *ca1, integer *icf1, integer *icf2, integer *icf11, integer *ipr, integer *nbc)
{
  /* System generated locals */
  integer icf1_dim1, icf2_dim1, icf11_dim1, 
    a1_dim1, a1_dim2, a2_dim1, a2_dim2, 
    cc_dim1, cc_dim2, faa_dim1, 
    ca1_dim1, ca1_dim2, ipr_dim1;

  /* Local variables */
  integer niam, nlev;
  real xkwt;
  integer nbcp1, ipiv1, ipiv2, i;

  integer i1, i2, k1, l1, ic, ir;
  doublereal rm;
  logical master[KREDO];
  integer myleft[KREDO];
  logical worker[KREDO];
  doublereal buf[2];
  integer ism[KREDO], irm[KREDO];
  doublereal tmp;
  logical notsend;
  integer nap1, nam1, myright[KREDO], icp1;

  
    /* Parameter adjustments */
    /*--fc;*/
  ipr_dim1 = *nov;
  icf11_dim1 = *nov;
  icf2_dim1 = *nov;
  icf1_dim1 = *nov;
  ca1_dim1 = *nov;
  ca1_dim2 = *nov;
  faa_dim1 = *nov;
  a2_dim1 = *nov;
  a2_dim2 = *nov;
  a1_dim1 = *nov;
  a1_dim2 = *nov;
  cc_dim1 = *nov;
  cc_dim2 = *nrc;

    
  nbcp1 = *nbc + 1;
  nap1 = *na + 1;
  nam1 = *na - 1;
  xkwt = (real) (*kwt);
  {
    real tmp = r_lg10(xkwt) / r_lg10(2.0);
    nlev = i_nint(&tmp);
  }
  notsend = TRUE_;

/* At each recursive level determine the master node (holding the pivot */
/* row after swapping), which will send the pivot row to the worker node 
*/
/* at distance 2**(K-1) from the master. Here K is the recursion level. */

  if (*par) {
    for (i = 0; i < nlev; ++i) {
      master[i] = FALSE_;
      worker[i] = FALSE_;
      k1 = pow_ii(2, i);
      niam = *iam / k1;
      if (notsend) {
	if (niam % 2 == 0) {
	  master[i] = TRUE_;
	  notsend = FALSE_;
	  ism[i] = (i + 1) + *iam + 10000;
	  irm[i] = ism[i] + k1;
	  myright[i] = *iam + k1;
	} else {
	  worker[i] = TRUE_;
	  ism[i] = (i + 1) + *iam + 10000;
	  irm[i] = ism[i] - k1;
	  myleft[i] = *iam - k1;
	}
      }
    }
  }

  /* Reduce concurrently in each node */
  for (i1 = 0; i1 < nam1; ++i1) {
    i2 = i1 + 1;
    for (ic = 0; ic < *nov; ++ic) {
      icp1 = ic + 1;
      ipiv1 = ARRAY2D(ipr, ic, i1);
      if (ipiv1 <= *nov) {
	tmp = ARRAY2D(faa, ic, i1);
	ARRAY2D(faa, ic, i1) = ARRAY2D(faa, (ipiv1 - 1), i1);
	ARRAY2D(faa, (ipiv1 - 1), i1) = tmp;
      } else {
	l1 = (ipiv1 - *nov) - 1;
	tmp = ARRAY2D(faa, ic, i1);
	ARRAY2D(faa, ic, i1) = ARRAY2D(faa, l1, i2);
	ARRAY2D(faa, l1, i2) = tmp;
      }
      for (ir = icp1; ir < *nov; ++ir) {
	l1 = ARRAY2D(icf2, ic, i1) - 1;
	rm = ARRAY3D(a2, ir, l1, i1);
	ARRAY2D(faa, ir, i1) -= rm * ARRAY2D(faa, ic, i1);
      }
      for (ir = 0; ir < *nov; ++ir) {
	l1 = ARRAY2D(icf1, ic, i2) - 1;
	rm = ARRAY3D(a1, ir, l1, i2);
	ARRAY2D(faa, ir, i2) -= rm * ARRAY2D(faa, ic, i1);
      }
      for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
	l1 = ARRAY2D(icf2, ic, i1) - 1;
	rm = ARRAY3D(cc, l1, ir, i2);
	fc[ir] -= rm * ARRAY2D(faa, ic, i1);
      }
    }
  }

  /* Inter-node reduction needs communication between nodes */
  if (*par) {
    for (i = 0; i < nlev; ++i) {
      for (ic = 0; ic < *nov; ++ic) {
	icp1 = ic + 1;
	if (master[i]) {
	  ipiv1 = ARRAY2D(ipr, ic, (*na - 1));
	  if (ipiv1 <= *nov) {
	    buf[0] = ARRAY2D(faa, (ipiv1 - 1), (*na - 1));
	    ARRAY2D(faa, (ipiv1 - 1), *na) = ARRAY2D(faa, ic, (*na - 1));
	    ARRAY2D(faa, ic, (*na - 1)) = buf[0];
	    buf[1] = -1.;
	    csend();
	  } else {
	    buf[0] = ARRAY2D(faa, ic, (*na - 1));
	    buf[1] = (doublereal) (ARRAY2D(ipr, ic, (*na - 1)) - *nov);
	    csend();
	    crecv();
	  }

	  for (ir = icp1; ir < *nov; ++ir) {
	    l1 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
	    rm = ARRAY3D(a2, ir, l1, (*na - 1));
	    ARRAY2D(faa, ir, (*na - 1)) -= rm * ARRAY2D(faa, ic, (*na - 1));
	  }
	  for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
	    l1 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
	    rm = ARRAY3D(cc, l1, ir, (nap1 - 1));
	    fc[ir] -= rm * ARRAY2D(faa, ic, (*na - 1));
	  }
	}

	if (worker[i]) {
	  crecv();
	  ipiv2 = i_dnnt(&buf[1]);
	  if (ipiv2 < 0) {
	    tmp = buf[0];
	  } else {
	    tmp = ARRAY2D(faa, (ipiv2 - 1), (*na - 1));
	    ARRAY2D(faa, (ipiv2 - 1), (*na - 1)) = buf[0];
	    csend();
	  }

	  for (ir = 0; ir < *nov; ++ir) {
	    l1 = ARRAY2D(icf11, ic, i) - 1;
	    rm = ARRAY3D(ca1, ir, l1, i);
	    ARRAY2D(faa, ir, (*na - 1)) -= rm * tmp;
	  }
	}
      }
      /*           **Synchronization at each recursion level among all n
		   odes */


    }

    l1 = *nrc - *nbc;
    gdsum();

  }

  return 0;
} /* redrhs_ */


/*     ---------- ------ */
/* Subroutine */ int 
dimrge(integer *iam, integer *kwt, logical *par, doublereal *e, doublereal *cc, doublereal *d, doublereal *fc, integer *ifst, integer *na, integer *nrc, integer *nov, integer *ncb, integer *idb, integer *nllv, doublereal *fcc, doublereal *p0, doublereal *p1, doublereal *det, doublereal *s, doublereal *a2, doublereal *faa, doublereal *bb)
{
    

  /* System generated locals */
  integer e_dim1, cc_dim1, cc_dim2, d_dim1, 
    p0_dim1, p1_dim1, s_dim1, s_dim2, 
    faa_dim1, a2_dim1, a2_dim2, bb_dim1, 
    bb_dim2;

  /* Local variables */

  integer i, j, k;

  integer novpi, novpj, k1, k2;

  integer novpj2, kc, kr, ncrloc, msglen1, msglen2, nap1;

  double *xe;
  xe = (doublereal *)malloc(sizeof(doublereal)*(*nov + *nrc));


  /* Parameter adjustments */
  /*--fc;*/
  /*--xe;*/
  /*--fcc;*/
  faa_dim1 = *nov;
  a2_dim1 = *nov;
  a2_dim2 = *nov;
  s_dim1 = *nov;
  s_dim2 = *nov;
  p1_dim1 = *nov;
  p0_dim1 = *nov;
  cc_dim1 = *nov;
  cc_dim2 = *nrc;
  e_dim1 = *nov + *nrc;
  bb_dim1 = *nov;
  bb_dim2 = *ncb;
  d_dim1 = *ncb;
    
  nap1 = *na + 1;
  msglen1 = (*nrc * 8) * *nov;
  /* Computing 2nd power */
  msglen2 = (*nov + *nrc + ((*nov * *nov) * 2) + 1) * 8;
  ncrloc = *nrc + *nov;

  /* Send CC(1:NOV,1:NRC,1) in node 0 to node KWT-1 */

  if (*par) {
    if (*iam == 0) {
      csend();
    }
    if (*iam == *kwt - 1) {
      crecv();
    }
  }

  /* Copy */
  if (*iam == *kwt - 1) {
    for (i = 0; i < *nov; ++i) {
      for (j = 0; j < *nov; ++j) {
	novpj = *nov + j;
	ARRAY2D(e, i, j) = ARRAY3D(s, i, j, (*na - 1));
	ARRAY2D(p0, i, j) = ARRAY3D(s, i, j, (*na - 1));
	ARRAY2D(e, i, novpj) = ARRAY3D(a2, i, j, (*na - 1));
	ARRAY2D(p1, i, j) = ARRAY3D(a2, i, j, (*na - 1));
      }
      for (j = 0; j < *ncb; ++j) {
	novpj2 = (*nov * 2) + j;
	ARRAY2D(e, i, novpj2) = ARRAY3D(bb, i, j, (*na - 1));
      }
    }

    for (i = 0; i < *nrc; ++i) {
      novpi = *nov + i;
      for (j = 0; j < *nov; ++j) {
	novpj = *nov + j;
	ARRAY2D(e, novpi, j) = ARRAY3D(cc, j, i, 0);
	ARRAY2D(e, novpi, novpj) = ARRAY3D(cc, j, i, (nap1 - 1));
      }
      for (j = 0; j < *ncb; ++j) {
	novpj2 = (*nov * 2) + j;
	ARRAY2D(e, novpi, novpj2) = ARRAY2D(d, j, i);
      }
    }

    for (i = 0; i < *nov; ++i) {
      xe[i] = ARRAY2D(faa, i, (*na - 1));
    }

    for (i = 0; i < *nrc; ++i) {
      novpi = *nov + i;
      xe[novpi] = fc[i];
    }

    if (*idb >= 3) {
      fprintf(fp9," Residuals of reduced system:\n");	
	  
      fprintf(fp9," ");
      for (i = 0; i < ncrloc; ++i) {
	fprintf(fp9,"%11.3E",xe[i]);	
	if((i+ 1)%10==0)
	  fprintf(fp9,"\n ");
	    
      }
      fprintf(fp9,"\n");	
    }

    if (*idb >= 4) {

      fprintf(fp9," Reduced Jacobian matrix:\n");	
	      
      for (i = 0; i < ncrloc; ++i) {
	int total_printed = 0;
	for (j = 0; j < ncrloc; ++j) {
	  if((total_printed != 0)&&(total_printed % 10 == 0))
	    fprintf(fp9,"\n");	
	  fprintf(fp9," %11.3E",ARRAY2D(e, i, j));	
	  total_printed++;
	}
	fprintf(fp9,"\n");	
      }
    }

    /* Solve for FCC */
    if (*nllv == 0) {
      ge(ncrloc, ncrloc, e, 1, ncrloc, fcc, 
	 ncrloc, xe, det);
    } else if (*nllv > 0) {
      nlvc(ncrloc, ncrloc, *nllv, e, fcc);
    } else {
      for (i = 0; i < ncrloc - 1; ++i) {
	xe[i] = 0.;
      }
      xe[-1 + ncrloc] = 1.;
      ge(ncrloc, ncrloc, e, 1, ncrloc, fcc, 
	 ncrloc, xe, det);
    }
    if (*idb >= 4) {
      fprintf(fp9," Solution vector:\n");	
	  
      for (i = 0; i < ncrloc; ++i) {
	if((i!=0)&&(i%7==0))
	  fprintf(fp9,"\n");	
	fprintf(fp9," %11.3E",fcc[i]);	
      }
      fprintf(fp9,"\n");	
    }

    k1 = ncrloc;
    /* Computing 2nd power */
    k2 = k1 + (*nov) * (*nov);
    for (kr = 0; kr < *nov; ++kr) {
      for (kc = 0; kc < *nov; ++kc) {
	k = kr * *nov + kc;
	fcc[k1 + k] = ARRAY2D(p0, kr, kc);
	fcc[k2 + k] = ARRAY2D(p1, kr, kc);
      }
    }
    /* Computing 2nd power */
    fcc[ncrloc + ((*nov) * (*nov) * 2)] = *det;

  }

  /* Broadcast FCC from node KWT-1. The matrices P0 and P1 are */
  /* buffered in the tail of FCC so all nodes receive them. */
  if (*par) {
    if (*iam == *kwt - 1) {
      csend();
    } else {
      crecv();
    }
  }

  for (i = 0; i < *nrc; ++i) {
    fc[i] = fcc[*nov + i];
  }

  if (*iam < *kwt - 1) {
    k1 = ncrloc;
    /* Computing 2nd power */
    k2 = k1 + (*nov) * (*nov);
    for (kr = 1; kr <= *nov; ++kr) {
      for (kc = 1; kc <= *nov; ++kc) {
	k = kr * *nov + kc;
	ARRAY2D(p0, kr, kc) = fcc[k1 + k];
	ARRAY2D(p1, kr, kc) = fcc[k2 + k];
      }
    }
    /* Computing 2nd power */
    *det = fcc[ncrloc + ((*nov) * (*nov) * 2)];
  }
  /* free the memory*/
  /* Not the we have modified these parameter before, so
     we undo the modifications here and then free them. */
  /*xe \+= 1;*/
  free(xe);

  return 0;
} /* dimrge_ */


/*     ---------- ------ */
/* Subroutine */ int 
bcksub(integer *iam, integer *kwt, logical *par, doublereal *s1, doublereal *s2, doublereal *a2, doublereal *bb, doublereal *faa, doublereal *fc, doublereal *fcc, doublereal *sol1, doublereal *sol2, doublereal *sol3, integer *na, integer *nov, integer *ncb, integer *icf2)
{
  /* System generated locals */
  integer icf2_dim1, s1_dim1, s1_dim2, s2_dim1, 
    s2_dim2, a2_dim1, a2_dim2, bb_dim1, bb_dim2,
    sol1_dim1, sol2_dim1, 
    sol3_dim1, faa_dim1;

    /* Local variables */
  integer niam, ibuf;
  logical even = FALSE_;
  integer nlev;
  logical hasright;
  doublereal xkwt;
  integer rmsgtype, smsgtype, i, k, l;

  integer nlist[2], itest, l1, l2;
  doublereal sm;
  integer msglen;

  logical master[KREDO];
  integer myleft, kp1;
  logical odd = FALSE_;
  integer ism, irm;
  logical hasleft, notsend;
  integer nam1, myright, nov2, nov3;
  double *buf=NULL;


    /* Parameter adjustments */
    /*--fc;*/
    /*--fcc;*/
  icf2_dim1 = *nov;
  sol3_dim1 = *nov;
  sol2_dim1 = *nov;
  sol1_dim1 = *nov;
  faa_dim1 = *nov;
  a2_dim1 = *nov;
  a2_dim2 = *nov;
  s2_dim1 = *nov;
  s2_dim2 = *nov;
  s1_dim1 = *nov;
  s1_dim2 = *nov;
  bb_dim1 = *nov;
  bb_dim2 = *ncb;
    
  xkwt = (doublereal) (*kwt);
  {
    doublereal tmp = d_lg10(&xkwt) / r_lg10(2.0);
    nlev = i_dnnt(&tmp);
  }
  nov2 = *nov * 2;
  nov3 = *nov * 3;
  ibuf = (nov3 + 1) * 8;

  /* The backsubstitution in the reduction process is recursive. */
  notsend = TRUE_;

  /*At each recursion level determine the sender nodes (called MASTER here).
*/
  if (*par) {
    for (i = 0; i < nlev; ++i) {
      master[i] = FALSE_;
      niam = *iam / pow_ii(2, i);
      if (notsend) {
	if (niam % 2 == 0) {
	  master[i] = TRUE_;
	  notsend = FALSE_;
	}
      }
    }
  }

  if (*par) {

    /*Initialization for the master or sender node at the last recursion l
      evel.*/
    if (master[nlev - 1]) {
      for (l = 0; l < *nov; ++l) {
	ARRAY2D(sol1, l, (*na - 1)) = fcc[l];
	ARRAY2D(sol3, l, (*na - 1)) = fc[l];
      }
    }

    for (i = nlev - 1; i >= 0; --i) {
      if (master[i]) {
	ism = i + nlev + (*kwt * 4);
	irm = ism + 1;
	k = pow_ii(2, i - 1);
	/*              **Compute the ID of the receiving node */
	nlist[0] = *iam - k;
	nlist[1] = *iam + k;
	/*              **Receive solutions from previous level */
	if ((i + 1) < nlev) {
	  crecv();
	  niam = i_dnnt(&buf[nov3 + 1]);
	  if (*iam < niam) {
	    for (l = 0; l < *nov; ++l) {
	      ARRAY2D(sol1, l, (*na - 1)) = buf[l + 1];
	      ARRAY2D(sol3, l, (*na - 1)) = buf[*nov + l + 1];
	    }
	  } else {
	    for (l = 0; l < *nov; ++l) {
	      ARRAY2D(sol1, l, (*na - 1)) = buf[*nov + l + 1];
	      ARRAY2D(sol3, l, (*na - 1)) = buf[nov2 + l + 1];
	    }
	  }
	}
	/*              **Backsubstitute */
	for (k = *nov - 1; k >= 0; --k) {
	  kp1 = k + 1;
	  sm = 0.;
	  for (l = 0; l < *nov; ++l) {
	    sm += ARRAY3D(s1, k, l, (*na - 1)) * ARRAY2D(sol1, l, (*na - 1));
	    sm += ARRAY3D(s2, k, l, (*na - 1)) * ARRAY2D(sol3, l, (*na - 1));
	  }
	  for (l = 0; l < *ncb; ++l) {
	    sm += ARRAY3D(bb, k, l, (*na - 1)) * fc[*nov + l];
	  }
	  for (l = kp1; l < *nov; ++l) {
	    l1 = ARRAY2D(icf2, l, (*na - 1)) - 1;
	    sm += ARRAY2D(sol2, l1, (*na - 1)) * ARRAY3D(a2, k, l1, (*na - 1));
	  }
	  l2 = ARRAY2D(icf2, k, (*na - 1)) - 1;
	  ARRAY2D(sol2, l2, (*na - 1)) = (ARRAY2D(faa, k, (*na - 1)) - sm) / ARRAY3D(a2, k, l2, (*na - 1));
	}
	/*              **Send solutions to the next level */
	if (i + 1 > 1) {
	  for (l = 0; l < *nov; ++l) {
	    buf[l + 1] = ARRAY2D(sol1, l, (*na - 1));
	    buf[*nov + l + 1] = ARRAY2D(sol2, l, (*na - 1));
	    buf[nov2 + l + 1] = ARRAY2D(sol3, l, (*na - 1));
	  }
	  buf[nov3 + 1] = (doublereal) (*iam);
	  gsendx();
	}
      }
      /*           **Synchronization at each recursion level */


    }

    /* Define odd and even nodes */
    if (*iam % 2 == 0) {
      even = TRUE_;
    } else {
      odd = TRUE_;
    }

    /* Determine whether I have a right neighbor */
    if (*iam == *kwt - 1) {
      hasright = FALSE_;
    } else {
      hasright = TRUE_;
    }

    /* Determine whether I have a left neighbor */
    if (*iam == 0) {
      hasleft = FALSE_;
    } else {
      hasleft = TRUE_;
    }

    /* Define send message type */
    smsgtype = *iam + 1000;

    /* Define receive message type */
    rmsgtype = smsgtype - 1;

    /* Define my right neighbor */
    myleft = *iam - 1;
    myright = *iam + 1;
    msglen = *nov << 3;

    /* May only need odd sends to even */
    itest = 0;
    if (itest == 1) {
      if (odd && hasright) {
	csend();
      }
      if (even && hasleft) {
	crecv();
      }
    }

    /* Even nodes send and odd nodes receive */
    if (even && hasright) {
      csend();
    }
    if (odd && hasleft) {
      crecv();
    }

  } else {

    for (l = 0; l < *nov; ++l) {
      ARRAY2D(sol1, l, (*na - 1)) = fcc[l];
      ARRAY2D(sol2, l, (*na - 1)) = fc[l];
    }

  }

  if (*iam == *kwt - 1) {
    for (l = 0; l < *nov; ++l) {
      ARRAY2D(sol2, l, (*na - 1)) = fc[l];
    }
  }

  if (*na > 1) {
    for (l = 0; l < *nov; ++l) {
      ARRAY2D(sol1, l, (*na - 2)) = ARRAY2D(sol1, l, (*na - 1));
      ARRAY2D(sol3, l, (*na - 2)) = ARRAY2D(sol2, l, (*na - 1));
    }
  }

  /* Backsubstitution process; concurrently in each node. */
  nam1 = *na - 1;
  for (i = nam1 - 1; i >= 0; --i) {
    for (k = *nov - 1; k >= 0; --k) {
      sm = 0.;
      for (l = 0; l < *nov; ++l) {
	sm += ARRAY2D(sol1, l, i) * ARRAY3D(s1, k, l, i);
	sm += ARRAY2D(sol3, l, i) * ARRAY3D(s2, k, l, i);
      }
      for (l = 0; l < *ncb; ++l) {
	sm += fc[*nov + l] * ARRAY3D(bb, k, l, i);
      }
      for (l = k + 1; l < *nov; ++l) {
	l1 = ARRAY2D(icf2, l, i) - 1;
	sm += ARRAY2D(sol2, l1, i) * ARRAY3D(a2, k, l1, i);
      }
      l2 = ARRAY2D(icf2, k, i) - 1;
      ARRAY2D(sol2, l2, i) = (ARRAY2D(faa, k, i) - sm) / ARRAY3D(a2, k, l2, i);
    }
    for (l = 0; l < *nov; ++l) {
      ARRAY2D(sol1, l, (i + 1)) = ARRAY2D(sol2, l, i);
      if (i + 1 > 1) {
	ARRAY2D(sol3, l, (i - 1)) = ARRAY2D(sol2, l, i);
	ARRAY2D(sol1, l, (i - 1)) = ARRAY2D(sol1, l, i);
      }
    }
  }

  return 0;
} /* bcksub_ */


/*     ---------- ------ */
/* Subroutine */ int 
infpar(integer *iam, logical *par, doublereal *a, doublereal *b, doublereal *fa, doublereal *sol1, doublereal *sol2, doublereal *fc, integer *na, integer *nov, integer *nra, integer *nca, integer *ncb, integer *irf, integer *icf)
{
  /* System generated locals */
  integer irf_dim1, icf_dim1, a_dim1, a_dim2, 
    b_dim1, b_dim2, fa_dim1, sol1_dim1,
    sol2_dim1;

    /* Local variables */
  integer nram, icfj1, i, j;
  doublereal *x;
  integer nrapj, irfir, j1, novpj, icfnovpir, ir;
  doublereal sm;
  integer novpir, irp1;

  x = (doublereal *)malloc(sizeof(doublereal)*(*nra));


/* Determine the local varables by backsubstitition. */

    /* Parameter adjustments */
    /*--fc;*/
  sol2_dim1 = *nov;
  sol1_dim1 = *nov;
  irf_dim1 = *nra;
  fa_dim1 = *nra;
  icf_dim1 = *nca;
  a_dim1 = *nca;
  a_dim2 = *nra;
  b_dim1 = *ncb;
  b_dim2 = *nra;
    
  nram = *nra - *nov;

/* Backsubstitution in the condensation of parameters; no communication. 
*/
  for (i = 0; i < *na; ++i) {
    for (ir = nram - 1; ir >= 0; --ir) {
      irp1 = ir + 1;
      sm = 0.;
      irfir = ARRAY2D(irf, ir, i) - 1;
      for (j = 0; j < *nov; ++j) {
	nrapj = *nra + j;
	sm += ARRAY3D(a, j, irfir, i) * ARRAY2D(sol1, j, i);
	sm += ARRAY3D(a, nrapj, irfir, i) * ARRAY2D(sol2, j, i);
      }
      for (j = 0; j < *ncb; ++j) {
	novpj = *nov + j;
	sm += ARRAY3D(b, j, irfir , i) * fc[novpj];
      }
      for (j = irp1; j < nram; ++j) {
	j1 = j + *nov;
	icfj1 = ARRAY2D(icf, j1, i) - 1;
	sm += ARRAY3D(a, icfj1, irfir, i) * x[icfj1];
      }
      novpir = *nov + ir;
      icfnovpir = ARRAY2D(icf, novpir, i) - 1;
      x[icfnovpir] = (ARRAY2D(fa, irfir, i) - sm) / ARRAY3D(a, icfnovpir, irfir, i);
    }
    /*        **Copy SOL1 and X into FA */
    for (j = 0; j < *nov; ++j) {
      ARRAY2D(fa, j, i) = ARRAY2D(sol1, j, i);
    }
    for (j = *nov; j < *nra; ++j) {
      ARRAY2D(fa, j, i) = x[j];
    }
  }
  free(x);

  return 0;
} /* infpar_ */


/*     ---------- --- */
/* Subroutine */ int 
rd0(integer *iam, integer *kwt, doublereal *d, integer *nrc)
{

  /* Local variables */
  integer niam;
  logical even[KREDO];
  doublereal xkwt;
  integer i, n;

  integer nredo, msglen, rmtype[KREDO], smtype[KREDO];
  logical odd[KREDO];

  doublereal *buf;

  logical notsend;
  integer myright[KREDO];

  buf = (doublereal *)malloc(sizeof(doublereal)*(*nrc));

/*     RECURSIVE DOUBLING PROCEDURE TO GET */
/*     THE GLOBAL SUM OF VECTORS FROM */
/*     EACH NODE. THE GLOBAL SUM IS ONLY AVAILABLE */
/*     IN THE LAST NODE */

/* Copying */
    /* Parameter adjustments */
    /*--d;*/

    
  xkwt = (doublereal) (*kwt);

  /* Determine the recursion level */
  {
    doublereal tmp = log(xkwt) / log((double)2.);
    nredo = i_dnnt(&tmp);
  }

/* At each recursion level determine the odd and even nodes */
  notsend = TRUE_;
  for (n = 0; n < nredo; ++n) {
    smtype[n] = n + 1000 + *iam + 1;
    rmtype[n] = smtype[n] - pow_ii(2, n);
    myright[n] = *iam + pow_ii(2, n);
    even[n] = FALSE_;
    odd[n] = FALSE_;
    niam = *iam / pow_ii(2, n);
    if (notsend) {
      if (niam % 2 == 0) {
	even[n] = TRUE_;
	notsend = FALSE_;
      } else {
	odd[n] = TRUE_;
      }
    }
  }

  niam = *nrc;
  msglen = niam * 8;
  for (n = 0; n < nredo; ++n) {
    /*        **Even nodes send and odd nodes receive from left to right 
     */
    if (even[n]) {
      csend();
    }
    if (odd[n]) {
      crecv();
      /*          ** Accumulate the partial sum in the current receiving
		  node */
      for (i = 0; i < niam; ++i) {
	d[i] += buf[i];
      }
    }
  }
  free(buf);
  return 0;
} /* rd0_ */

/*     ---------- ------ */
/* Subroutine */ int 
print1(integer *nov, integer *na, integer *nra, integer *nca, integer *ncb, integer *nrc, doublereal *a, doublereal *b, doublereal *c, doublereal *d, doublereal *fa, doublereal *fc)
{
    

  /* System generated locals */
  integer a_dim1, a_dim2, b_dim1, b_dim2, c_dim1, 
    c_dim2, d_dim1, fa_dim1;

  /* Local variables */
  integer i, ic, ir;

  /* Parameter adjustments */
  /*--fc;*/
  fa_dim1 = *nra;
  a_dim1 = *nca;
  a_dim2 = *nra;
  d_dim1 = *ncb;
  b_dim1 = *ncb;
  b_dim2 = *nra;
  c_dim1 = *nca;
  c_dim2 = *nrc;
    
  fprintf(fp9,"AA , BB , FA (Full dimension) :\n");	
  /* should be 10.3f*/
  for (i = 0; i < *na; ++i) {
    fprintf(fp9,"I=%3ld\n",i + 1);
    for (ir = 0; ir < *nra; ++ir) {
      int total_written = 0;
      for (ic = 0; ic < *nca; ++ic) {
	if((total_written != 0) && (total_written%12 == 0))
	  fprintf(fp9,"\n");
	fprintf(fp9," %10.3E",ARRAY3D(a, ic, ir, i));
	total_written++;
      }
      for (ic = 0; ic < *ncb; ++ic) {
	if((total_written != 0) && (total_written%12 == 0))
	  fprintf(fp9,"\n");
	fprintf(fp9," %10.3E",ARRAY3D(b, ic, ir, i));	
	total_written++;
      }
      if((total_written != 0) && (total_written%12 == 0))
	fprintf(fp9,"\n");
      fprintf(fp9," %10.3E",ARRAY2D(fa, ir, i));	
      fprintf(fp9,"\n");	
    }
  }

  fprintf(fp9,"CC (Full dimension) :\n");	

  for (i = 0; i < *na; ++i) {
    fprintf(fp9,"I=%3ld\n",i + 1);	
    for (ir = 0; ir < *nrc; ++ir) {
      int total_written = 0;
      for (ic = 0; ic < *nca; ++ic) {
	if((total_written != 0) && (total_written%12 == 0))
	  fprintf(fp9,"\n");
	fprintf(fp9," %10.3E",ARRAY3D(c, ic, ir, i));	
	total_written++;
      }
      fprintf(fp9,"\n");	
    }
  }

  fprintf(fp9,"DD , FC\n");	

  for (ir = 0; ir < *nrc; ++ir) {
    int total_written = 0;
    for (ic = 0; ic < *ncb; ++ic) {
      if((total_written != 0) && (total_written%12 == 0))
	fprintf(fp9,"\n");
      fprintf(fp9," %10.3E",ARRAY2D(d, ic, ir));	
      total_written++;
    }
    fprintf(fp9," %10.3E\n",fc[ir]);	
  }


  return 0;
} /* print1_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*         Dummy Routines for the Sequential Version */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
integer 
mynode(void)
{
  integer ret_val;
  ret_val = 0;
  return ret_val;
}

integer 
numnodes(void)
{
  integer ret_val;
  ret_val = 1;
  return ret_val;
}


/* Subroutine */ int 
gsync(void)
{
  return 0;
} /* gsync_ */

doublereal 
dclock(void)
{
  real ret_val;

  ret_val = (double)0.;
  return ret_val;
} 


/* Subroutine */ int 
csend(void)
{
  return 0;
} /* csend_ */


/* Subroutine */ int 
crecv(void)
{
  return 0;
} /* crecv_ */


/* Subroutine */ int 
gdsum(void)
{
  return 0;
} /* gdsum_ */


/* Subroutine */ int 
gsendx(void)
{
  return 0;
} /* gsendx_ */


/* Subroutine */ int 
gcol(void)
{
  return 0;
} /* gcol_ */


/* Subroutine */ int 
led(void)
{
  return 0;
} /* led_ */


/* Subroutine */ int 
setiomode(void)
{
  return 0;
} /* setiomode_ */






































