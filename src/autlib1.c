#include "auto_f2c.h"
#include "auto_c.h"
#include "xAuto.h"
extern XAUTO xAuto;
extern int NODE;
extern int RestartLabel;
int restart_flag=0;
/* The memory for these are taken care of in main, and setubv for the
   mpi parallel case.  These are global since the they are used many times
   in the wrapper functions in autlib3.c (and autlib5.c) and the cost
   of allocating and deallocating them is prohibitive. */
struct {
  doublereal *dfu, *dfp, *uu1, *uu2, *ff1, *ff2;
} global_scratch={NULL,NULL,NULL,NULL,NULL,NULL};

/* The memory for these are taken care of in main, and setubv for the
   mpi parallel case.  These are global since they only need to be
   computed once for an entire run, so we do them at the
   beginning to save the cost later on. */
struct {
  integer irtn;
  integer *nrtn;
} global_rotations = {0,NULL};

/* There are used to short circuit the code.  getp is a user callable function
   that allows certain parameters to be returned.  Unfortunately, the
   data that this function works on is NOT user accessible, so cannot
   be part of its calling sequence.  Accordingly, this global structure is
   filled in with the necessary data so that getp has access to it when the
   user calls that routine. */
struct {
  rap_type *rav;
  iap_type *iav;
  doublereal *dtv;
} global_parameters = {NULL,NULL,NULL};

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*                    Initialization */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
FILE *fp8;
int fp8_is_open=0;
extern char fort8[200],fort3[200];
void send_eigen();
void send_mult();
int byeauto_();
void addbif(iap_type *iap, rap_type *rap, integer ntots, integer ibrs, double *par,integer *icp,int labw, double *a, double *uhigh, double *ulow, double *u0, double *ubar);

int init(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *thl, doublereal **thu_pointer, integer *iuz, doublereal *vuz)
{
  /* Local variables */
  doublereal hbff, biff;
  integer nbif;
  doublereal fldf;
  integer iads, ndim, nicp, ncol, mxbf, ithl;
  doublereal epsl;
  integer nthl, nfpr, nins, iplt, nint, jtmp;
  doublereal epsu, epss;
  integer nthu, itmx, itnw, ithu, ntot, ipos, nwtn;
  doublereal tivp;
  integer ntst, iuzr;
  doublereal spbf;
  integer nuzr, i;
  doublereal dsold, dsmin, dsmax, a0, a1;
  integer istop, itpst;
  doublereal ds, rl0, rl1;
  integer iad, jac, lab, nbc, iid, ibr, ndm;
  doublereal amp, det;
  integer ilp, nit, ips, isp, irs, itp, npr, isw, nmx, nbc0, nnt0;
  doublereal *thu;


  
  for (i = 0; i < NPARX; ++i) {
    icp[i] = i;
    jtmp = NPARX;
    icp[jtmp + i] = 0;
    par[i] = 0.;
    par[jtmp + i] = 0.;
    thl[icp[i]] = 1.;
  }

 ndim=xAuto.ndim;
  ips=xAuto.ips;
  irs=xAuto.irs;
  ilp=xAuto.ilp;

  thu = *thu_pointer = (doublereal *)malloc(sizeof(doublereal)*8*ndim);

  for (i = 0; i < ndim * 8; ++i) {
    thu[i] = 1.;
  }

  jtmp = NPARX;
  nicp=xAuto.nicp;

  for (i = 0; i < nicp; ++i) {
    icp[jtmp+i]=xAuto.icp[i];

  }
  
  if (nicp > 0) {
    for (i = 0; i < nicp; ++i) {
      jtmp = NPARX;
      icp[i] = icp[jtmp + i];
    }
  } else {
    nicp = 1;
    jtmp = NPARX;
    icp[jtmp] = icp[0];
  }
      
  ntst=xAuto.ntst;
  ncol=xAuto.ncol;
  iad=xAuto.iad;
  isp=xAuto.isp;
  isw=xAuto.isw;
  iplt=xAuto.iplt;
  nbc=xAuto.nbc;
  nint=xAuto.nint;
  nmx=xAuto.nmx;
  rl0=xAuto.rl0;
  rl1=xAuto.rl1;
  a0=xAuto.a0;
  a1=xAuto.a1;  
  npr=xAuto.npr;
  mxbf=xAuto.mxbf;
  iid=xAuto.iid;
  itmx=xAuto.itmx;
  itnw=xAuto.itnw;
  nwtn=xAuto.nwtn;
  jac=xAuto.jac;
  epsl=xAuto.epsl;
  epss=xAuto.epss;
  epsu=xAuto.epsu;
  ds=xAuto.ds;
  dsmin=xAuto.dsmin;
  dsmax=xAuto.dsmax;
  iads=xAuto.iads;
  
  
  if(dsmin < 0.0) {
    printf("Warning : DSMIN less then 0.0, will use absolute value instead.");
    dsmin = fabs(dsmin);
  }

  if(dsmax < 0.0) {
    printf("Warning : DSMAX less then 0.0, will use absolute value instead.");
    dsmax = fabs(dsmax);
  }
  nthl=xAuto.nthl;
  
  if (nthl > 0) {
    for (i = 0; i < nthl; ++i) 
      thl[xAuto.ithl[i]]=xAuto.thl[i];
    
  }
    
  nthu=0;      
  nuzr=xAuto.nuzr;
  for(i=0;i<nuzr;i++)
    {    
      iuz[i]=xAuto.iuz[i];
      vuz[i]=xAuto.vuz[i];
    }

  iap->ndim = ndim;

  iap->ips = ips;

  iap->irs = irs;

  iap->ilp = ilp;

  iap->ntst = ntst;

  iap->ncol = ncol;

  iap->iad = iad;

  iap->iads = iads;

  iap->isp = isp;

  iap->isw = isw;

  iap->iplt = iplt;

  iap->nbc = nbc;

  iap->nint = nint;

  iap->nmx = nmx;

  iap->nuzr = nuzr;

  iap->npr = npr;

  iap->mxbf = mxbf;

  iap->iid = iid;

  iap->itmx = itmx;

  iap->itnw = itnw;

  iap->nwtn = nwtn;

  iap->jac = jac;


  ndm = ndim;

  if (nbc != 0) {

    nbc0 = nbc;

  } else {

    nbc0 = 1;

  }

  if (nint != 0) {

    nnt0 = nint;

  } else {

    nnt0 = 1;

  }

  iuzr = 1;

  itp = 0;

  itpst = 0;

  nfpr = 1;

  ibr = 1;

  nit = 0;

  ntot = 0;

  nins = 0;

  istop = 0;

  nbif = 0;

  ipos = 1;

  lab = 0;


  iap->ndm = ndm;

  iap->nbc0 = nbc0;

  iap->nnt0 = nnt0;

  iap->iuzr = iuzr;

  iap->itp = itp;

  iap->itpst = itpst;

  iap->nfpr = nfpr;

  iap->ibr = ibr;

  iap->nit = nit;

  iap->ntot = ntot;

  iap->nins = nins;

  iap->istop = istop;

  iap->nbif = nbif;

  iap->ipos = ipos;

  iap->lab = lab;

  iap->nicp = nicp;


  rap->ds = ds;

  rap->dsmin = dsmin;

  rap->dsmax = dsmax;

  dsold = ds;

  rap->dsold = dsold;

  rap->rl0 = rl0;

  rap->rl1 = rl1;

  rap->a0 = a0;

  rap->a1 = a1;


  amp = 0.;

  det = 0.;

  tivp = 0.;

  fldf = 0.;

  hbff = 0.;

  biff = 0.;

  spbf = 0.;


  rap->amp = amp;

  rap->epsl = epsl;

  rap->epsu = epsu;

  rap->epss = epss;

  rap->det = det;

  rap->tivp = tivp;

  rap->fldf = fldf;

  rap->hbff = hbff;

  rap->biff = biff;

  rap->spbf = spbf;




  return 0;

} 





/*     ---------- ----- */
/* Subroutine */ int 
chdim(iap_type *iap)
{
  /* Local variables */
  integer npar;


/* Check dimensions. */


    
  npar = iap->nfpr;

  if (npar > NPARX) {
    if (iap->mynode == 0) {
      printf("Dimension exceeded : NPAR=%5ld  maximum=%5d (Increase NPARX in auto.h and recompile AUTO",npar,NPARX);
    }
    exit(0);
  }

  return 0;
} /* chdim_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*               The leading subroutines of AUTO */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

int autoae(iap_type *iap, rap_type *rap, doublereal *par, 
	   integer *icp, 
	   FUNI_TYPE((*funi)),
	   STPNT_TYPE_AE((*stpnt)),
	   PVLI_TYPE_AE((*pvli)),
	   doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)

{
  /* This is the entry subroutine for algebraic systems. */

  cnrlae(iap, rap, par, icp, funi, stpnt, pvli, thl, thu, iuz, vuz);
  return 0;
} 


int 
autobv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), STPNT_TYPE_BVP((*stpnt)), PVLI_TYPE_BVP((*pvli)), doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{


  /* THIS IS THE ENTRY ROUTINE FOR GENERAL BOUNDARY VALUE PROBLEMS. */

  cnrlbv(iap, rap, par, icp, funi, bcni, icni, stpnt, pvli, thl, thu, iuz, vuz);
  return 0;
} 

/*     ---------- ----- */
/* Subroutine */ int 
init1(iap_type *iap, rap_type *rap, integer *icp, doublereal *par)
{

    /* Local variables */
  integer ndim, ncol, nicp;

  integer nneg, nfpr, nint, jtmp, nuzr, i;
  doublereal dsmin, dsmax, fc;
  integer ic, jc;
  int ii;
  doublereal ds;
  integer nxploc, jac, nbc, ndm, ict[NPARX], ilp, ips, isp, irs, itp,
    isw, nmx;

  /* General initialization. Redefinition of constants. */


    
  ndim = iap->ndim;
  ips = iap->ips;
  irs = iap->irs;
  ilp = iap->ilp;
  ncol = iap->ncol;
  isp = iap->isp;
  isw = iap->isw;
  nbc = iap->nbc;
  nint = iap->nint;
  nmx = iap->nmx;
  nuzr = iap->nuzr;
  jac = iap->jac;
  itp = iap->itp;
  nfpr = iap->nfpr;
  nicp = iap->nicp;

  ds = rap->ds;
  dsmin = rap->dsmin;
  dsmax = rap->dsmax;

  if (isw == 0) {
    isw = 1;
  }

  /* Check and perturb pseudo arclength stepsize and steplimits. */
  /* (Perturbed to avoid exact computation of certain singular points). */

  if (ds == 0.) {
    ds = (double).1;
  }
  if (dsmin == 0.) {
    dsmin = fabs(ds) * 1e-4;
  }
  fc = HMACH1;
  ds = fc * ds;
  dsmin /= fc;
  dsmax = fc * dsmax;

  /* Redefinition for waves */
  if (ips == 11) {
    ips = 1;
    iap->ips = ips;
    ndim <<= 1;
    ndm = ndim;
    iap->ndm = ndm;
  } else if (ips == 12) {
    ips = 2;
    iap->ips = ips;
    ndim <<= 1;
    ndm = ndim;
    iap->ndm = ndm;
  }

  /* General Redefinition. */

  if (abs(ips) <= 1 && isw == 1) {
    /*        ** Algebraic Systems */
    nfpr = 1;

  } else if (ips == -2) {
    /*        ** Time integration */
    nfpr = 1;
    isp = 0;
    ilp = 0;
    icp[0] = 13;

  } else if (ips == 2 && abs(isw) == 1) {
    /*        ** Periodic Solutions */
    nbc = ndim;
    nint = 1;
    nfpr = nbc + nint - ndim + 1;
    /*        **ISW=1 when starting from a HB */
    if (itp == 3 || abs(itp) / 10 == 3) {
      isw = 1;
    }
    if (nicp == 1) {
      /*          **Variable period */
      icp[1] = 10;
    }

  } else if (ips == 4 && abs(isw) == 1) {
    /*        ** Boundary value problems */
    nfpr = nbc + nint - ndim + 1;

  } else if (ips == 7 && abs(isw) == 1) {
    /*        ** Boundary value problems */
    nfpr = nbc + nint - ndim + 1;

  } else if (ips == 9 && abs(isw) == 1) {
    /*        ** Homoclinic bifurcation analysis */
    /*        Redefine AUTO constants for homoclinic orbits */
    inho(iap, icp, par);
    ndim = iap->ndim;
    nbc = iap->nbc;
    nint = iap->nint;
    nuzr = iap->nuzr;
    nfpr = nbc + nint - ndim + 1;
    /*    printf("%d %d %d %d %d\n",ndim,nbc,nint,nuzr,nfpr);
     for(ii=0;ii<nfpr;ii++)
	 printf("%d %d \n",ii,icp[ii]);
    */
  } else if (ips == 14 || ips == 16) {
    /*        **Evolution calculations for Parabolic Systems */
    ndim <<= 1;
    nbc = ndim;
    nint = 0;
    nfpr = 1;
    ilp = 0;
    isp = 0;
    icp[0] = 13;

  } else if (ips == 17) {
    /*        **Stationary calculations for Parabolic Systems */
    ndim <<= 1;
    nbc = ndim;
    nint = 0;
    nfpr = 1;

  } else if (ips == 15) {
    /*          ** Optimization of periodic solutions */
    for (i = 0; i < nicp; ++i) {
      ict[i] = icp[i];
    }
    nfpr = 0;
    for (i = 0; i < nicp; ++i) {
      if (ict[i] >= 0) {
	icp[nfpr] = ict[i];
	++nfpr;
      }
    }
    icp[nfpr] = 9;
    icp[nfpr + 1] = 12;
    icp[nfpr + 2] = 13;
    nfpr += 3;
    ndim <<= 1;
    nbc = ndim;
    nint = nfpr - 1;
    /* overload to define optimality integrals */
    nneg = 0;
    for (i = 0; i < nicp; ++i) {
      ic = ict[i];
      jc = abs(ic) - 20;
      if (ic < 0 && jc > 0 && jc <= 11) {
	++nneg;
	icp[nfpr + nneg - 1] = jc;
      }
    }
    /* Set indices of output parameters */
    nicp = nfpr - 3;
    for (i = 0; i < nicp; ++i) {
      jtmp = NPARX;
      icp[jtmp + i] = icp[i];
    }

  } else if (ips == 5) {
    /*        ** Algebraic optimization Problems */
    if (iap->itp % 10 == 2 || iap->irs == 0) {
      iap->nfpr++;
    }
    nfpr = iap->nfpr;
    if (nfpr == 2) {
      ++ndim;
      icp[0] = 9;
    } else {
      ndim = (ndim << 1) + nfpr;
      icp[0] = 9;
    }

  } else if (irs > 0 && abs(isw) == 2) {
    /*        ** Continuation of singular points */

    if ((itp == 1 || abs(itp) / 10 == 1 || itp == 2 || abs(itp) / 10 == 2)
	&& abs(ips) <= 1) {
      /*          ** Fold continuation (Algebraic Problems) */
      ndim = (ndim << 1) + 1;
      nfpr = 2;

    } else if ((itp == 3 || abs(itp) / 10 == 3) && abs(ips) <= 1) {
      /*          ** Hopf bifurcation continuation (Maps, ODE, Waves) */
      ndim = ndim * 3 + 2;
      nfpr = 2;

    } else if ((itp == 5 || itp == 6) && ips == 2) {
      /*          ** Fold continuation (Periodic solutions); start */
      ndim <<= 1;
      nbc = ndim;
      nint = 3;
      nfpr = nbc + nint - ndim + 1;
      if (icp[2] == 10 || nicp == 2) {
	/*            ** Variable period */
	icp[1] = 12;
	icp[2] = 10;
	icp[3] = 11;
      } else {
	/*            ** Fixed period */
	icp[2] = 12;
	icp[3] = 11;
      }
      ilp = 0;
      isw = -2;
      isp = 0;
      nmx = 5;
      if (iap->mynode == 0) {
	printf("\nGenerating starting data :\n Restart at EP label below :\n");
	restart_flag=1;
      }

    } else if ((abs(itp) / 10 == 5 || abs(itp) / 10 == 6) && ips == 2) {
      /*          ** Fold continuation (Periodic solutions); restart */
      ndim <<= 1;
      nbc = ndim;
      nint = 3;
      nfpr = nbc + nint - ndim + 1;
      if (nicp == 2) {
	/*            ** Variable period */
	icp[2] = 10;
      }
      icp[3] = 11;

    } else if (itp == 7 && ips == 2) {
      /*          ** Continuation of period doubling bifurcations; start
 */
      ndim <<= 1;
      nbc = ndim;
      nint = 2;
      nfpr = nbc + nint - ndim + 1;
      if (icp[2] == 10 || nicp == 2) {
	/*            ** Variable period */
	icp[1] = 10;
	icp[2] = 12;
      } else {
	/*            ** Fixed period */
	icp[2] = 12;
      }
      ilp = 0;
      isw = -2;
      isp = 0;
      nmx = 5;
      if (iap->mynode == 0) {
	printf("\nGenerating starting data :\n Restart at EP label below :\n");
	restart_flag=1;
      }

    } else if (abs(itp) / 10 == 7 && ips == 2) {
      /*          ** Continuation of period doubling bifurcations; resta
rt */
      ndim <<= 1;
      nbc = ndim;
      nint = 2;
      nfpr = nbc + nint - ndim + 1;
      if (icp[2] == 10 || nicp == 2) {
	/*            ** Variable period */
	icp[2] = 10;
      }

    } else if (itp == 8 && ips == 2) {
      /*          ** Continuation of torus bifurcations; start */
      ndim *= 3;
      nbc = ndim;
      nint = 3;
      nfpr = nbc + nint - ndim + 1;
      icp[1] = 10;
      icp[2] = 11;
      icp[3] = 12;
      ilp = 0;
      isp = 0;
      isw = -2;
      nmx = 5;
      if (iap->mynode == 0) {
	printf("\nGenerating starting data :\n Restart at EP label below :\n");
	restart_flag=1;
      }

    } else if (abs(itp) / 10 == 8 && ips == 2) {
      /*          ** Continuation of torus bifurcations; restart */
      ndim *= 3;
      nbc = ndim;
      nint = 3;
      nfpr = nbc + nint - ndim + 1;
      icp[2] = 10;
      icp[3] = 11;

    } else if ((itp == 5 || itp == 6) && ips == 4) {
      /*          ** Continuation of folds (BVP; start) */
      ndim <<= 1;
      nbc <<= 1;
      nint = (nint << 1) + 1;
      nfpr = nbc + nint - ndim + 1;
      nxploc = nfpr / 2 - 1;
      if (nxploc > 0) {
	for (i = 0; i < nxploc; ++i) {
	  icp[nfpr / 2 + i] = i + 10;
	}
      }
      icp[nfpr / 2] = nfpr / 2 + 10;
      ilp = 0;
      isw = -2;
      isp = 0;
      nmx = 5; 
      if (iap->mynode == 0) {
	restart_flag=1;
	printf("\nGenerating starting data :\n Restart at EP label below :\n");
      }

    } else if ((abs(itp) / 10 == 5 || abs(itp) / 10 == 5) && ips == 4) {
      /*          ** Continuation of folds (BVP; restart) */
      ndim <<= 1;
      nbc <<= 1;
      nint = (nint << 1) + 1;
      nfpr = nbc + nint - ndim + 1;
      nxploc = nfpr / 2 - 1;
      if (nxploc > 0) {
	for (i = 0; i < nxploc; ++i) {
	  icp[nfpr / 2 + i] = i + 10;
	}
      }

    }

  }

  iap->ndim = ndim;
  iap->ips = ips;
  iap->irs = irs;
  iap->ilp = ilp;
  iap->ncol = ncol;
  iap->isp = isp;
  iap->isw = isw;
  iap->nbc = nbc;
  iap->nint = nint;
  iap->nmx = nmx;
  iap->nuzr = nuzr;
  iap->jac = jac;
  iap->nfpr = nfpr;
  iap->nicp = nicp;

  rap->ds = ds;
  rap->dsmin = dsmin;
  rap->dsmax = dsmax;
  /* printf("ready to go \n");
  printf("dim=%d ps=%d rs=%d lp=%d sp=%d sw=%d nbc=%d nint=%d nfpr=%d nicp=%d\n",
  ndim,ips,irs,ilp,isp,isw,nbc,nint,nfpr,nicp); */
  return 0;
} /* init1 */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*                    Algebraic Problems */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ------ */
/* Subroutine */ int 
cnrlae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), STPNT_TYPE_AE((*stpnt)), PVLI_TYPE_AE((*pvli)), doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{

    /* Local variables */
  integer nbfc;
  doublereal *dfdp;
  integer nbif, iads;
  doublereal *dfdu;
  integer mxbf;

  doublereal *uold, stla[NBIFX], stld[NBIFX];
  integer nins;
  doublereal *udot;
  integer ipos;
  doublereal *stud;	
  integer ntot, iuzr, nuzr;
  doublereal *f;
  integer i, k;
  doublereal *u, dsold;

  doublereal rlold[NPARX], rldot[NPARX], rlcur[NPARX];
  integer istop, itpst;

  doublereal *aa;

  doublereal ds, *du;

  integer lab, ibr, ndm;
  doublereal rbp, rds;
  integer ips;
  doublereal *rhs;
  integer irs, isp;
  doublereal rev;
  integer isw;
  doublereal rlp;
  integer nit, itp;
  doublereal *stu, *uzr;

  integer aa_first_dimension = iap->ndim + 1;

  dfdp = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*NPARX);
  dfdu = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*(iap->ndim));
  uold = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  udot = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  stud =  (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*NBIFX);
  f = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  u = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  aa = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim+1)*(iap->ndim+1));
  du = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim+1));
  rhs = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim+1));
  stu =  (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*NBIFX);
  uzr =  (doublereal *)malloc(sizeof(doublereal)*(iap->nuzr));
    
  /* Controls the bifurcation analysis of algebraic problems */


/* Local */


  
  ips = iap->ips;
  irs = iap->irs;
  iads = iap->iads;
  isp = iap->isp;
  isw = iap->isw;
  nuzr = iap->nuzr;
  mxbf = iap->mxbf;
  ndm = iap->ndm;
  itpst = iap->itpst;
  ibr = iap->ibr;

  ds = rap->ds;

  nins = 0;
  iap->nins = nins;
  rbp = 0.;
  rev = 0.;
  rlp = 0.;
  if (nuzr > 0) {
    for (i = 0; i < nuzr; ++i) {
      uzr[i] = 0.;
    }
  }
  rds = ds;
  dsold = ds;
  rap->dsold = dsold;
  nit = 0;
  iap->nit = nit;
  nbif = 0;
  iap->nbif = nbif;
  nbfc = 0;
  ipos = 1;
  iap->ipos = ipos;
  ntot = 0;
  iap->ntot = ntot;
  lab = 0;
  iap->lab = lab;

  for (i = 0; i < iap->ndim; ++i) {
    u[i] = 0.;
    du[i] = 0.;
    udot[i] = 0.;
    uold[i] = 0.;
    f[i] = 0.;
  }

  /* Generate the starting point */

  (*stpnt)(iap, rap, par, icp, u);
  (*pvli)(iap, rap, u, par);

/* Determine a suitable starting label and branch number */

  newlab(iap, rap);

  /* Write constants */

  sthd(iap, rap, par, icp, thl, thu);

  /* Write plotting data for the starting point */

  istop = 0;
  iap->istop = istop;
  if (irs == 0) {
    itp = itpst * 10 + 9;
  } else {
    itp = 0;
  }
  iap->itp = itp;
  rlcur[0] = par[icp[0]];
  stplae(iap, rap, par, icp, rlcur, u);
  istop = iap->istop;
  if (istop == 1) {
    goto L6;
  }

  /* Starting procedure  (to get second point on first branch) : */

  stprae(iap, rap, par, icp, funi, &rds, &aa_first_dimension, aa, 
	 rhs, rlcur, rlold, rldot, u, du, uold, udot, f, dfdu, dfdp, thl, 
	 thu);
  istop = iap->istop;
  if (istop == 1) {
    goto L5;
  }
  itp = 0;
  iap->itp = itp;
  goto L3;

/* Initialize computation of the next bifurcating branch. */

 L2:
  swpnt(iap, rap, par, icp, &rds, NBIFX, stud, stu, 
	stla, stld, rlcur, rlold, rldot, u, udot);

  ipos = iap->ipos;
  if (ipos == 1) {
    --nbif;
    iap->nbif = nbif;
    ++nbfc;
  }

  rbp = 0.;
  rev = 0.;
  rlp = 0.;
  if (nuzr > 0) {
    for (i = 0; i < nuzr; ++i) {
      uzr[i] = 0.;
    }
  }
  if (ipos == 0 || mxbf < 0) {
    ++ibr;
  }
  iap->ibr = ibr;

  ntot = 0;
  iap->ntot = ntot;
  istop = 0;
  iap->istop = istop;
  itp = 0;
  iap->itp = itp;
  nit = 0;
  iap->nit = nit;
  dsold = rds;
  rap->dsold = dsold;

/* Store plotting data for first point on the bifurcating branch */

  stplae(iap, rap, par, icp, rlcur, u);
  istop = iap->istop;
  if (istop == 1) {
    goto L6;
  }

  /* Determine the second point on the bifurcating branch */

  swprc(iap, rap, par, icp, funi, &aa_first_dimension, aa, rhs, 
	rlcur, rlold, rldot, u, du, uold, udot, f, dfdu, dfdp, &rds, thl, 
	thu);
  istop = iap->istop;
  if (istop == 1) {
    goto L5;
  }

  /* Store plotting data for second point : */

  stplae(iap, rap, par, icp, rlcur, u);
  istop = iap->istop;
  if (istop == 1) {
    goto L6;
  }
  rbp = 0.;
  rev = 0.;
  rlp = 0.;

/* Provide initial approximation to the next point on the branch */

 L3:
  contae(iap, rap, &rds, rlcur, rlold, rldot, u, uold, udot);

/* Find the next solution point on the branch */

  solvae(iap, rap, par, icp, funi, &rds, &aa_first_dimension, aa, 
	 rhs, rlcur, rlold, rldot, u, du, uold, udot, f, dfdu, dfdp, thl, 
	 thu);
  istop = iap->istop;
  if (istop == 1) {
    goto L5;
  }

  /* Check for user supplied parameter output parameter-values. */

  if (nuzr > 0) {
    for (iuzr = 0; iuzr < nuzr; ++iuzr) {
      iap->iuzr = iuzr;
      lcspae(iap, rap, par, icp, 
	     fnuzae, funi, &aa_first_dimension, aa, rhs, rlcur, rlold, rldot, u, du, uold, udot,
	     f, dfdu, dfdp, &uzr[iuzr], thl, thu, iuz, vuz);
      istop = iap->istop;
      if (istop == 1) {
	goto L5;
      }
      itp = iap->itp;
      if (itp == -1) {
	if (iuz[iuzr] >= 0) {
	  itp = -4 - itpst * 10;
	  iap->itp = itp;
	  for (k = 0; k < nuzr; ++k) {
	    uzr[k] = 0.;
	  }
	} else {
	  istop = -1
;
	  iap->istop = istop;
	}
      }
    }
  }

  /* Check for fold */

  if (iap->ilp == 1) {
    lcspae(iap, rap, par, icp, fnlpae, 
	   funi, &aa_first_dimension, 
	   aa, rhs, rlcur, rlold, rldot, u, du, uold, udot, f, dfdu, 
	   dfdp, &rlp, thl, thu, iuz, vuz);
    itp = iap->itp;
    if (itp == -1) {
      itp = itpst * 10 + 2;
      iap->itp = itp;
      rlp = 0.;
      rbp = 0.;
      rev = 0.;
    }
  }

  /* Check for branch point, and if so store data : */

  if (isp != 0) {
    lcspae(iap, rap, par, icp, 
	   fnbpae, funi, &aa_first_dimension, 
	   aa, rhs, rlcur, rlold, rldot, u, du, uold, udot, f, dfdu, 
	   dfdp, &rbp, thl, thu, iuz, vuz);
    istop = iap->istop;
    if (istop == 1) {
      goto L5;
    }
    itp = iap->itp;
    if (itp == -1) {
      itp = itpst * 10 + 1;
      iap->itp = itp;
      ++nbif;
      iap->nbif = nbif;
      stbif(iap, rap, par, icp, &aa_first_dimension, aa,
	    NBIFX, stud, stu, stla, stld, rlcur, rlold, rldot, u, 
	    du, udot, dfdu, dfdp, thl, thu);
      rlp = 0.;
      rbp = 0.;
      rev = 0.;
    }
  }

  /* Check for Hopf bifurcation */

  if (abs(ips) == 1) {
    lcspae(iap, rap, par, icp, fnhbae, 
	   funi, &aa_first_dimension, 
	   aa, rhs, rlcur, rlold, rldot, u, du, uold, udot, f, dfdu, 
	   dfdp, &rev, thl, thu, iuz, vuz);
    istop = iap->istop;
    if (istop == 1) {
      goto L5;
    }
    itp = iap->itp;
    if (itp == -1) {
      itp = itpst * 10 + 3;
      iap->itp = itp;
      rev = 0.;
    }
  }

  /* Store plotting data on unit 7 : */

 L5:
  stplae(iap, rap, par, icp, rlcur, u);

/* Adapt the stepsize along the branch */

  itp = iap->itp;
  ntot = iap->ntot;
  if (iads != 0 && ntot % iads == 0 && (itp % 10 == 0 || itp % 10 == 4)) {
    adptds(iap, rap, &rds);
  }

 L6:
  itp = 0;
  iap->itp = itp;
  istop = iap->istop;
  if (istop == 0) {
    goto L3;
  }

  nbif = iap->nbif;
  if (nbif != 0 && nbfc < abs(mxbf)) {
    goto L2;
  }

  free(dfdp);    
  free(dfdu);
  free(uold);
  free(udot);
  free(stud);
  free(f);
  free(u);
  free(aa);
  free(du);
  free(rhs);
  free(stu);
  free(uzr);

  return 0;
} /* cnrlae_ */


/*     ---------- ------ */
/* Subroutine */ int 
stpnus(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *u)
{
  integer ndim;



  

  /* Gets the starting data from user supplied STPNT */


  
  ndim = iap->ndim;

  stpnt(ndim, 0.0, u, par);

  return 0;
} /* stpnus_ */


/*     ---------- ------ */
/* Subroutine */ int 
stpnae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *u)
{
  logical found;
  integer nfprs;

  integer irs;

  


  /* Gets the starting data from unit 3 */



  
  irs = iap->irs;
  findlb(iap, rap, irs, &nfprs, &found);
  readlb(iap, rap, u, par);
  return 0;
} /* stpnae_ */


/*     ---------- ------ */
/* Subroutine */ int 
stprae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), doublereal *rds, integer *m1aaloc, doublereal *aa, doublereal *rhs, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *du, doublereal *uold, doublereal *udot, doublereal *f, doublereal *dfdu, doublereal *dfdp, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer aa_dim1;

  /* Local variables */
  integer ndim;

  doublereal sign;
  integer i, k;

  doublereal sc;
  doublereal ss;

  integer iid;





/* Finds the second point on the initial solution branch. */

  aa_dim1 = *m1aaloc;

  
  ndim = iap->ndim;
  iid = iap->iid;

  rlold[0] = par[icp[0]];
  for (i = 0; i < ndim; ++i) {
    uold[i] = u[i];
  }

  /* Determine the direction of the branch at the starting point */

  (*funi)(iap, rap, ndim, u, uold, icp, par, 2, f, dfdu, dfdp);
  for (i = 0; i < ndim; ++i) {
    rhs[i] = f[i];
    ARRAY2D(aa, i, ndim) = dfdp[(icp[0]) * ndim + i];
    ARRAY2D(aa, ndim, i) = 0.;
    for (k = 0; k < ndim; ++k) {
      ARRAY2D(aa, i, k) = dfdu[k * ndim + i];
    }
  }
  rhs[ndim] = 0.;
  ARRAY2D(aa, ndim, ndim) = 0.;

  if (iid >= 3) {
    integer tmp = ndim + 1;
    wrjac(iap, &tmp, m1aaloc, aa, rhs);
  }
  nlvc(ndim + 1, *m1aaloc, 1, aa, du);

  /* Scale and make sure that the PAR(ICP(1))-dot is positive. */

  ss = 0.;
  for (i = 0; i < ndim; ++i) {
    /* Computing 2nd power */
    ss += thu[i] * (du[i] * du[i]);
  }
  /* Computing 2nd power */
  ss += thl[icp[0]] * (du[ndim] * du[ndim]);

  sign = 1.;
  if (du[ndim] < 0.) {
    sign = -1.;
  }
  sc = sign / sqrt(ss);
  for (i = 0; i < ndim + 1; ++i) {
    du[i] = sc * du[i];
  }

  for (i = 0; i < ndim; ++i) {
    udot[i] = du[i];
  }
  rldot[0] = du[ndim];

/* Set initial approximations to the second point on the branch */

  for (i = 0; i < ndim; ++i) {
    u[i] = uold[i] + *rds * udot[i];
  }
  rlcur[0] = rlold[0] + *rds * rldot[0];

  solvae(iap, rap, par, icp, funi, rds, m1aaloc, aa,
	 rhs, rlcur, rlold, rldot, u, du, uold, 
	 udot, f, dfdu, dfdp, thl, thu);

  return 0;
} /* stprae_ */


/*     ---------- ------ */
/* Subroutine */ int 
contae(iap_type *iap, rap_type *rap, doublereal *rds, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *uold, doublereal *udot)
{
  /* Local variables */
  integer ndim, i;
  doublereal dsold;
  integer ips;


  

/* This subroutine determines an initial approximation to the next */
/* solution on a branch by extrapolating from the two preceding points. */
/* The step used in the preceding step has been stored in DSOLD. */



    
  ndim = iap->ndim;
  ips = iap->ips;

  dsold = rap->dsold;

  rldot[0] = (rlcur[0] - rlold[0]) / dsold;
  for (i = 0; i < ndim; ++i) {
    udot[i] = (u[i] - uold[i]) / dsold;
  }

  rlold[0] = rlcur[0];
  rlcur[0] += *rds * rldot[0];
  for (i = 0; i < ndim; ++i) {
    uold[i] = u[i];
    u[i] += udot[i] * *rds;
  }
  /*      Save old time for time integration */
  if (ips == -2) {
    rap->tivp = rlold[0];
  }

  return 0;
} /* contae_ */


/*     ---------- ----- */
/* Subroutine */ int 
solvae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), doublereal *rds, integer *m1aaloc, doublereal *aa, doublereal *rhs, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *du, doublereal *uold, doublereal *udot, doublereal *f, doublereal *dfdu, doublereal *dfdp, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer aa_dim1;

  /* Local variables */
  integer iads, ndim;
  doublereal drlm;
  integer ndmr;
  doublereal epsl, epsu, dumx;
  integer ntop, itnw, ntot, i, k;
  doublereal dsold, dsmin;

  doublereal dsmax;
  doublereal rdrlm, rdumx;
  integer istop;

  doublereal au;
  doublereal delref=0.0, ss, delmax;

  integer iid;
  doublereal dds;
  integer ibr, ndm;
  doublereal det, adu;
  integer nit, mxt;
  doublereal umx;
  integer nit1;
  static integer last_ntop=0;

/* This is the subroutine for computing solution branches. It solves */
/* the equations for finding the next point on the branch at distance DS 
*/
/* from the current point. An initial approximation to the new point */
/* ( i.e. to PAR(ICP(1)) and U ) has been supplied by CONT. */


/* Local */

  aa_dim1 = *m1aaloc;

  
  ndim = iap->ndim;
  iads = iap->iads;
  iid = iap->iid;
  itnw = iap->itnw;
  ndm = iap->ndm;
  ibr = iap->ibr;

  dsmin = rap->dsmin;
  dsmax = rap->dsmax;
  epsl = rap->epsl;
  epsu = rap->epsu;

 L1:
  dsold = *rds;
  rap->dsold = dsold;
  dds = 1. / *rds;
  nit = 0;
  iap->nit = nit;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;
  ndmr = ndim;
  if (ndmr > 6) {
    ndmr = 6;
  }
  if (iid >= 2 && iap->mynode == 0) {
    if(last_ntop != ntop) {
      fprintf(fp9,"========================================");
      fprintf(fp9,"========================================\n");
      last_ntop = ntop;
    }
    if (nit == 0) {
      fprintf(fp9,"  BR    PT  IT\n");
    }

    fprintf(fp9,"%4li%6li%4li    %14.6E              ",ibr,ntop,
	    nit,rlcur[0]);
    for (i = 0; i < ndmr; ++i) {
      fprintf(fp9,"%14.6E",u[i]);
    }
    fprintf(fp9,"\n");
  }

  /* Call user-supplied FUNC to evaluate the right hand side of the */
  /* differential equation and its derivatives : */

  for (nit1 = 1; nit1 <= itnw; ++nit1) {

    nit = nit1;
    iap->nit = nit;
    par[icp[0]] = rlcur[0];
    (*funi)(iap, rap, ndim, u, uold, icp, par, 2, f, dfdu, dfdp);

    /* Set up the Jacobian matrix and the right hand side : */

    for (i = 0; i < ndim; ++i) {
      ARRAY2D(aa, i, ndim) = dfdp[(icp[0]) * ndim + i];
      rhs[i] = -f[i];
      for (k = 0; k < ndim; ++k) {
	ARRAY2D(aa, i, k) = dfdu[k * ndim + i];
      }
    }
    for (k = 0; k < ndim; ++k) {
      ARRAY2D(aa, ndim, k) = thu[k] * 2. * (u[k] - uold[k]) * dds;
    }
    ARRAY2D(aa, ndim, ndim) = thl[icp[0]] * 2. * (rlcur[0] - 
						      rlold[0]) * dds;
    ss = 0.;
    for (i = 0; i < ndim; ++i) {
      /* Computing 2nd power */
      ss += thu[i] * ((u[i] - uold[i]) * (u[i] - uold[i]));
    }
    /* Computing 2nd power */
    rhs[ndim] = *rds - dds * ss - thl[icp[0]] * dds * ((rlcur[0] - rlold[0]) * (rlcur[0] - rlold[0]));

    /* Use Gauss elimination with pivoting to solve the linearized system 
       : */

    if (iid >= 5) {
      integer tmp = ndim + 1;
      wrjac(iap, &tmp, m1aaloc, aa, rhs);
    }

    ge(ndim + 1, *m1aaloc, aa, 1, ndim + 1, du, 
       ndim + 1, rhs, &det);
    rap->det = det;
    drlm = du[ndim];

    /* Add the Newton increments : */

    for (i = 0; i < ndim; ++i) {
      u[i] += du[i];
    }
    rlcur[0] += drlm;
    dumx = 0.;
    umx = 0.;
    for (i = 0; i < ndim; ++i) {
      adu = fabs(du[i]);
      au = fabs(u[i]);
      if (au > umx) {
	umx = au;
      }
      if (adu > dumx) {
	dumx = adu;
      }
    }

    if (iid >= 2 && iap->mynode == 0) {
      if (iap->mynode == 0) {

	fprintf(fp9,"%4li%6li%4li    %14.6E              ",ibr,ntop,
		nit,rlcur[0]);
	for (i = 0; i < ndmr; ++i) {
	  fprintf(fp9,"%14.6E",u[i]);
	}
	fprintf(fp9,"\n");
      }
    }

    rdrlm = fabs(drlm) / (fabs(rlcur[0]) + 1.);
    rdumx = dumx / (umx + 1.);
    if (rdrlm <= epsl && rdumx <= epsu) {
      pvlsae(iap, rap, u, par);
      if (iid >= 2) {
	fprintf(fp9,"\n");
      }
      return 0;
    }

    /* Check whether relative error has reached user-supplied tolerance : 
     */

    if (nit == 1) {
      delref = max(rdrlm,rdumx) * 20;
    } else {
      delmax = max(rdrlm,rdumx);
      if (delmax > delref) {
	goto L3;
      }
    }

    /* L2: */
  }

  /* Maximum number of iterations has been reached */

 L3:
  if (iads == 0 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li NOTE:No convergence with fixed step size\n",ibr,ntop);
  }
  if (iads == 0) {
    goto L5;
  }

  /* Reduce stepsize and try again */

  mxt = itnw;
  iap->nit = mxt;
  adptds(iap, rap, rds);
  if (fabs(*rds) < dsmin) {
    goto L4;
  }
  rlcur[0] = rlold[0] + *rds * rldot[0];
  for (i = 0; i < ndim; ++i) {
    u[i] = uold[i] + *rds * udot[i];
  }
  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9," NOTE:Retrying step\n");
  }
  goto L1;

/* Minimum stepsize reached */

 L4:
  if (iap->mynode == 0) {
    fprintf(fp9,"%4li%6li NOTE:No convergence using minimum step size\n",ibr,ntop);
  }
 L5:
  rlcur[0] = rlold[0];
  par[icp[0]] = rlcur[0];
  for (i = 0; i < ndim; ++i) {
    u[i] = uold[i];
  }
  istop = 1;
  iap->istop = istop;
  return 0;
} /* solvae_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*               Detection of Singular Points */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ------ */
/* Subroutine */ int 
lcspae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FNCS_TYPE_AE((*fncs)), FUNI_TYPE((*funi)), integer *m1aaloc, doublereal *aa, doublereal *rhs, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *du, doublereal *uold, doublereal *udot, doublereal *f, doublereal *dfdu, doublereal *dfdp, doublereal *q, doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{

  /* Local variables */
  logical chng;
  doublereal epss, rrds;
  integer ntop, itmx;
  doublereal rtmp;
  integer ntot;
  doublereal s, dsold, dsmax;
  integer istop;
  doublereal q0, q1, s0, s1, dq, ds, pq;

  integer itlcsp, iid, ibr;
  doublereal rds;
  integer itp;

  /* This subroutine uses the secant method to accurately locate special */
  /* points (branch points, folds, Hopf bifurcations, user zeroes). */
  /* These are characterized as zeroes of the function FNCS supplied in the 
*/
/* call. */
/* This subroutine calls CONT and SOLVAE with varying stepsize RDS. */
/* The special point is assumed to have been found with sufficient */
/* accuracy if the ratio between RDS and the user supplied value of */
/* DS is less than the user-supplied toler du. */

    
  iid = iap->iid;
  itmx = iap->itmx;
  ibr = iap->ibr;

  ds = rap->ds;
  dsmax = rap->dsmax;
  dsold = rap->dsold;
  epss = rap->epss;

  /* Check whether FNCS has changed sign (FNCS is EXTERNAL). */

  q0 = *q;
  q1 = (*fncs)(iap, rap, par, icp, &chng, funi, m1aaloc, 
	       aa, rlcur, rlold, rldot, u, uold, udot, 
	       rhs, dfdu, dfdp, iuz, vuz);
  pq = q0 * q1;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;
  if (pq >= 0. || ! chng) {
    *q = q1;
    return 0;
  }

  /* Use the secant method for the first step: */

  s0 = 0.;
  s1 = dsold;
  itlcsp = 0;
  dq = q0 - q1;
  rds = q1 / dq * (s1 - s0);
  rtmp = HMACH1;
 L1:
  rds = rtmp * rds;
  s = s1 + rds;

  /* Return if relative tolerance has been met : */

  rrds = fabs(rds) / (sqrt(fabs(ds * dsmax)) + 1);
  if (rrds < epss) {
    itp = -1;
    iap->itp = itp;
    *q = 0.;
    fprintf(fp9," ==> Location of special point :  Convergence.    Stepsize =%11.3E\n",rds);
    return 0;
  }

  /* If requested write additional output on unit 9 : */

  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9," ==> Location of special point :  Iteration %3li   Stepsize =%11.3E\n",itlcsp,rds);	
  }

  contae(iap, rap, &rds, rlcur, rlold, rldot, u, uold, 
	 udot);
  solvae(iap, rap, par, icp, funi, &rds, m1aaloc, aa, 
	 rhs, rlcur, rlold, rldot, u, du, uold, udot, f,
	 dfdu, dfdp, thl, thu);
  istop = iap->istop;
  if (istop == 1) {
    *q = 0.;
    return 0;
  }

  *q = (*fncs)(iap, rap, par, icp, &chng, funi, m1aaloc, 
	       aa, rlcur, rlold, rldot, u, uold, udot,
	       rhs, dfdu, dfdp, iuz, vuz);
  ++itlcsp;
  if (itlcsp <= itmx) {
    /*        Use Mueller's method with bracketing for subsequent steps */
    mueller(&q0, &q1, q, &s0, &s1, &s, &rds);
    goto L1;
  } else {
    if (iap->mynode == 0) {
      fprintf(fp9,"%4li%6li NOTE:Possible special point\n",ibr,(ntot + 1) % 10000);	
    }
    *q = 0.;
    return 0;
  }

} /* lcspae_ */


/*     ---------- ------- */
/* Subroutine */ int 
mueller(doublereal *q0, doublereal *q1, doublereal *q, doublereal *s0, doublereal *s1, doublereal *s, doublereal *rds)
{
    /* Local variables */
  doublereal a, b, c, d, r, h0, h1, dq;


  /* Mueller's method with bracketing */




  h0 = *s0 - *s;
  h1 = *s1 - *s;
  d = h0 * h1 * (h1 - h0);
  /* Computing 2nd power */

  a = (h1 * h1 * (*q0 - *q) - h0 *  h0 * (*q1 - *q)) / d;
  b = (-h1 * (*q0 - *q) + h0 * (*q1 - *q)) / d;
  if (fabs(b) <= RSMALL) {
    *rds = -(*q) / a;
  } else {
    c = a / (b * 2);
    /* Computing 2nd power */
    r = sqrt(c * c - *q / b);
    if (c < 0.) {
      *rds = -c - r;
    } else {
      *rds = -c + r;
    }
  }

  dq = *q1 * *q;
  if (dq < 0.) {
    *q0 = *q1;
    *s0 = *s1;
  }
  *q1 = *q;
  *s1 = *s;

  return 0;
} /* mueller_ */


/*     ------ --------- -------- ------ */
doublereal 
fnbpae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, logical *chng, FUNI_TYPE((*funi)), integer *m1aaloc, doublereal *aa, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *uold, doublereal *udot, doublereal *rhs, doublereal *dfdu, doublereal *dfdp, integer *iuz, doublereal *vuz)
{
  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  integer ntop, ntot, iid, ibr;
  doublereal det;


    
  iid = iap->iid;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  det = rap->det;
  ret_val = det;
  *chng = TRUE_;

  /* If requested write additional output on unit 9 : */

  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li        BP   Function %14.6E\n",ibr,ntop,ret_val);
  }

  return ret_val;
} /* fnbpae_ */


/*     ------ --------- -------- ------ */
doublereal 
fnlpae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, logical *chng, FUNI_TYPE((*funi)), integer *m1aaloc, doublereal *aa, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *uold, doublereal *udot, doublereal *rhs, doublereal *dfdu, doublereal *dfdp, integer *iuz, doublereal *vuz)
{
  /* System generated locals */
  integer aa_dim1;
  doublereal ret_val;

    

    /* Local variables */
  integer ndim, ntop, ntot, i, k;


  doublereal *ud;
  integer iid, ibr;
  doublereal det;

  ud = (double *)malloc(sizeof(doublereal)*(iap->ndim + 1));


  /* Local */


  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  /*--u;*/
  /*--uold;*/
  /*--udot;*/
  /*--rhs;*/
  /*--dfdu;*/
  /*--dfdp;*/
  aa_dim1 = *m1aaloc;

  
  ndim = iap->ndim;
  iid = iap->iid;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  par[icp[0]] = rlcur[0];
  (*funi)(iap, rap, ndim, u, uold, icp, par, 2,
	  rhs, dfdu, dfdp);
  for (i = 0; i < ndim; ++i) {
    ARRAY2D(aa, i, ndim) = dfdp[(icp[0]) * ndim + i];
    for (k = 0; k < ndim; ++k) {
      ARRAY2D(aa, i, k) = dfdu[k * ndim + i];
    }
  }
  for (k = 0; k < ndim; ++k) {
    ARRAY2D(aa, ndim, k) = udot[k];
    rhs[k] = 0.;
  }
  ARRAY2D(aa, ndim, ndim) = rldot[0];
  rhs[ndim] = 1.;

  ge(ndim + 1, *m1aaloc, aa, 1, ndim + 1, ud, ndim + 1, rhs, &det);
  rap->det = det;
  {
    integer tmp = ndim + 1;
    nrmlz(&tmp, ud);
  }
  ret_val = ud[ndim];
  rap->fldf = ret_val;
  *chng = TRUE_;

  /* If requested write additional output on unit 9 : */

  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li        Fold Function %14.6E\n",abs(ibr),ntop,ret_val);
  }
  free(ud);
  return ret_val;
} /* fnlpae_ */


/*     ------ --------- -------- ------ */
doublereal 
fnhbae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, logical *chng, FUNI_TYPE((*funi)), integer *m1aaloc, doublereal *aa, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *uold, doublereal *udot, doublereal *rhs, doublereal *dfdu, doublereal *dfdp, integer *iuz, doublereal *vuz)
{
  /* System generated locals */
  integer aa_dim1;
  doublereal ret_val;

  /* Local variables */
  integer ndim;
  doublereal arev, rmax;
  integer nins, ntop, ntot;
  doublecomplex ztmp;
  integer nins1, i, j;
  doublereal rimhb, ar;
  integer ntotp1;

  doublecomplex *ev;
  doublereal rp;

  integer iid, ibr, ndm, ier, loc=0, ips;
  doublereal rev;
  integer isw;

  ev = (doublecomplex *)malloc(sizeof(doublecomplex)*(iap->ndim));


  /* Local */

  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--u;*/
  /*--uold;*/
  /*--udot;*/
  /*--rhs;*/
  aa_dim1 = *m1aaloc;
  
  ndim = iap->ndim;
  ndm = iap->ndm;
  ips = iap->ips;
  isw = iap->isw;
  iid = iap->iid;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

/* INITIALIZE */

  *chng = FALSE_;

  /* Compute the eigenvalues of the Jacobian */

  eig(iap, &ndm, &ndim, dfdu, ev, &ier);
  if (ips == -1) {
    for (i = 0; i < ndm; ++i) {
      if (ev[i].r != -1. || d_imag(&ev[i]) != 0.) {
	doublecomplex in,out;
	in.r = ev[i].r + 1., in.i = ev[i].i;
	z_log(&out, &in);
	ev[i].r = out.r, ev[i].i = out.i;
      } else {
	ev[i].r = -RLARGE, ev[i].i = 0.;
      }
    }
  }
  /* here is where to put send_eigenvalue   */
  send_eigen(ibr,ntot+1,ndim,(doublecomplex*)&ev[0]);
  /* Order the eigenvalues by real part. */

  for (i = 0; i < ndm - 1; ++i) {
    rmax = -RLARGE;
    for (j = i; j < ndm; ++j) {
      rp = ev[j].r;
      if (rp >= rmax) {
	rmax = rp;
	loc = j;
      }
    }
    if (loc != i) {
      ztmp.r = ev[loc].r, ztmp.i = ev[loc].i;
      ev[loc].r = ev[i].r, ev[loc].i = ev[i].i;
      ev[i].r = ztmp.r, ev[i].i = ztmp.i;
    }
  }

  /* Compute the smallest real part. */

  rimhb = 0.;
  arev = RLARGE;
  rev = 0.;
  for (i = 0; i < ndm; ++i) {
    if (d_imag(&ev[i]) != 0.) {
      ar = fabs(ev[i].r);
      if (ar <= arev) {
	arev = ar;
	rev = ev[i].r;
	rimhb = fabs(d_imag(&ev[i]));
	if (rimhb != 0.) {
	  par[10] = pi(2.0) / rimhb;
	}
      }
    }
  }

  /* Compute the number of eigenvalues with negative real part. */

  nins1 = 0;
  if (isw != 2) {
    for (i = 0; i < ndm; ++i) {
      if (ev[i].r <= 0.) {
	++nins1;
      }
    }
  } else {
    for (i = 0; i < ndm; ++i) {
      if (ev[i].r <= HMACH) {
	++nins1;
      }
    }
  }

  if (isw == 2) {
    ret_val = 0.;
  } else {
    ret_val = rev;
  }
  rap->hbff = ret_val;
  nins = iap->nins;
  if (nins1 != nins) {
    *chng = TRUE_;
  }
  nins = nins1;
  iap->nins = nins;

  ntot = iap->ntot;
  ntotp1 = ntot + 1;
  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li        Hopf Function %14.6E\n",abs(ibr),ntop,ret_val);	
  }
  if (nins1 == ndm) {
    ntotp1 = -ntotp1;
  }

  if (iap->mynode == 0) {
    fprintf(fp9,"%4li%6li        Eigenvalues:                                 Stable:%3li\n",abs(ibr),ntop,nins);	
    if (ips == -1) {
      for (i = 0; i < ndm; ++i) {
	doublecomplex tmp;
	z_exp(&tmp, &ev[i]);
	fprintf(fp9,"%4li%6li        Eigenvalue%3li %14.6E%14.6E\n",abs(ibr),ntop,i+1,tmp.r,tmp.i);	
      }
    } else {
      for (i = 0; i < ndm; ++i) {
	fprintf(fp9,"%4li%6li        Eigenvalue%3li %14.6E%14.6E\n",abs(ibr),ntop,i+1,ev[i].r,ev[i].i);
      }
    }
  }

  free(ev);
  return ret_val;
} /* fnhbae_ */


/*     ------ --------- -------- ------ */
doublereal 
fnuzae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, logical *chng, FUNI_TYPE((*funi)), integer *m1aaloc, doublereal *aa, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *uold, doublereal *udot, doublereal *rhs, doublereal *dfdu, doublereal *dfdp, integer *iuz, doublereal *vuz)
{
  /* System generated locals */
  doublereal ret_val;

    

  /* Local variables */
  integer ntop, ntot, iuzr, iid, ibr;

  






  /* Parameter adjustments */
  /*--vuz;*/
  /*--iuz;*/
  /*--par;*/
    

  
  iid = iap->iid;
  iuzr = iap->iuzr;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  ret_val = par[abs(iuz[iuzr])] - vuz[iuzr];
  *chng = TRUE_;

  if (iid >= 3) {
    fprintf(fp9,"%4li%6li        User Func. %3li %16.6E\n",abs(ibr),ntop,iuzr,ret_val);	
  }

  return ret_val;
} /* fnuzae_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*                   Branch Switching for Algebraic Problems */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ----- */
/* Subroutine */ int 
stbif(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, integer *m1aaloc, doublereal *aa, integer m1sbloc, doublereal *stud, doublereal *stu, doublereal *stla, doublereal *stld, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *du, doublereal *udot, doublereal *dfdu, doublereal *dfdp, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer aa_dim1, stud_dim1, stu_dim1;

  /* Local variables */
  integer nbif, ndim;

  integer ntop, ntot, i, j;
  doublereal sc;
  doublereal ss;
  integer ibr;


  /* Stores branching data in the following arrays : */
  /*        STU    ( the solution vector U ) */
  /*        STUD   ( U-dot ) */
  /*        STLA   ( PAR(ICP(1)) ) */
  /*        STLD  ( PAR(ICP(1))-dot ) */
  /* Here the vector ( PAR(ICP(1))-dot , U-dot ) lies in the 2-d nullspace 
   */
  /* at branch point and is perpendicular to the direction vector of */
  /* known branch at this point. */

/* Local */

  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--stla;*/
  /*--stld;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  /*--u;*/
  /*--du;*/
  /*--udot;*/
  /*--dfdu;*/
  /*--dfdp;*/
  /*--thl;*/
  /*--thu;*/
  aa_dim1 = *m1aaloc;
  stu_dim1 = m1sbloc;
  stud_dim1 = m1sbloc;

  
  ndim = iap->ndim;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;
  nbif = iap->nbif;

  /* Keep track of the number of branch points stored. */

  if (nbif == NBIFX && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li NOTE:No more branch points can be stored\n",ibr,ntop);	
  }
  if (nbif > NBIFX) {
    nbif = NBIFX;
    iap->nbif = nbif;
    return 0;
  }

  for (i = 0; i < ndim; ++i) {
    for (j = 0; j < ndim; ++j) {
      ARRAY2D(aa, i, j) = dfdu[j * ndim + i];
    }
  }

  for (i = 0; i < ndim; ++i) {
    ARRAY2D(aa, i, ndim) = dfdp[(icp[0]) * ndim + i];
    ARRAY2D(aa, ndim, i) = udot[i];
  }
  ARRAY2D(aa, ndim, ndim) = rldot[0];

  nlvc(ndim + 1, *m1aaloc, 1, aa, du);

  ss = 0.;
  for (i = 0; i < ndim; ++i) {
    /* Computing 2nd power */
    ss += thu[i] * (du[i] * du[i]);
  }
  /* Computing 2nd power */
  ss += thl[icp[0]] * (du[ndim] * du[ndim]);
  sc = 1. / sqrt(ss);

  for (i = 0; i < ndim + 1; ++i) {
    du[i] = sc * du[i];
  }

  nbif = iap->nbif;
  stld[-1 + nbif] = du[ndim];
  for (i = 0; i < ndim; ++i) {
    ARRAY2D(stu, -1 + nbif, i) = u[i];
    ARRAY2D(stud, -1 + nbif, i) = du[i];
  }
  stla[-1 + nbif] = rlcur[0];

  return 0;
} /* stbif_ */


/*     ---------- ----- */
/* Subroutine */ int 
swpnt(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *rds, integer m1sbloc, doublereal *stud, doublereal *stu, doublereal *stla, doublereal *stld, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *udot)
{
  /* System generated locals */
  integer stud_dim1, stu_dim1;

    /* Local variables */
  integer nbif, ndim, mxbf, ipos, i, i1;
  doublereal ds;
  integer isw;


  

  /* This subroutine retrieves the branching data U, U-dot, PAR(ICP(1)), */
  /* PAR(ICP(1))-dot. If this initialization corresponds to the computation */
  /* of the bifurcating branch in opposite direction, then only the sign of */
  /*  the stepsize ( DS ) along the branch is reversed. */


  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--stla;*/
  /*--stld;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  /*--u;*/
  /*--udot;*/
  stu_dim1 = m1sbloc;
  stud_dim1 = m1sbloc;
  
  ndim = iap->ndim;
  isw = iap->isw;
  mxbf = iap->mxbf;
  nbif = iap->nbif;
  ipos = iap->ipos;

  ds = rap->ds;

  *rds = ds;
  if (ipos == 0) {
    *rds = -ds;
  }
  rlcur[0] = stla[0];
  par[icp[0]] = rlcur[0];
  rldot[0] = stld[0];
  for (i = 0; i < ndim; ++i) {
    u[i] = ARRAY2D(stu, 0, i);
    udot[i] = ARRAY2D(stud, 0, i);
  }
  if (abs(isw) == 2) {
    par[icp[1]] = u[-1 + ndim];
  }

  if (mxbf >= 0) {
    ipos = 1 - ipos;
    iap->ipos = ipos;
  }
  if (ipos == 0) {
    return 0;
  }

  for (i = 0; i < nbif; ++i) {
    stla[i] = stla[i + 1];
    stld[i] = stld[i + 1];
    for (i1 = 0; i1 < ndim; ++i1) {
      ARRAY2D(stu, i, i1) = ARRAY2D(stu, i + 1, i1);
      ARRAY2D(stud, i, i1) = ARRAY2D(stud, i + 1, i1);
    }
  }

  return 0;
} /* swpnt_ */


/*     ---------- ----- */
/* Subroutine */ int 
swprc(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), integer *m1aaloc, doublereal *aa, doublereal *rhs, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal *u, doublereal *du, doublereal *uold, doublereal *udot, doublereal *f, doublereal *dfdu, doublereal *dfdp, doublereal *rds, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer aa_dim1;

    /* Local variables */
  integer iads, ndim;
  doublereal drlm;
  integer ndmr;
  doublereal epsl, epsu, dumx;
  integer ntop, itnw, ntot, i, k;
  doublereal dsold, dsmin;

  doublereal dsmax, rdrlm, rdumx;
  integer istop;
  doublereal *u1;

  doublereal au;
  doublereal ss;

  integer iid;
  doublereal adu;
  integer ibr;
  doublereal det;
  integer nit, mxt;
  doublereal umx, rlm1;
  integer nit1;

  u1 = (double *)malloc(sizeof(doublereal)*(iap->ndim + 1));

  

  /* Controls the computation of the second point on a bifurcating branch. */
/* This point is required to lie in a hyper-plane at distance DS from the */
/* branch point. This hyper-plane is parallel to the tangent of the */
/* known branch at the branch point. */


/* Local */

  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--rhs;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  /*--u;*/
  /*--du;*/
  /*--uold;*/
  /*--udot;*/
  /*--f;*/
  /*--dfdu;*/
  /*--dfdp;*/
  /*--thl;*/
  /*--thu;*/
  aa_dim1 = *m1aaloc;
  
  ndim = iap->ndim;
  iads = iap->iads;
  iid = iap->iid;
  itnw = iap->itnw;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  dsmin = rap->dsmin;
  dsmax = rap->dsmax;
  epsl = rap->epsl;
  epsu = rap->epsu;

  /* Initialize and provide initial guess : */

  rlold[0] = rlcur[0];
  rlcur[0] = rlold[0] + *rds * rldot[0];
  for (i = 0; i < ndim; ++i) {
    uold[i] = u[i];
    u[i] = uold[i] + *rds * udot[i];
  }

 L2:
  dsold = *rds;
  rap->dsold = dsold;
  nit = 0;
  iap->nit = nit;

/* Write additional output on unit 9 if requested : */

  ndmr = ndim;
  if (ndmr > 6) {
    ndmr = 6;
  }
  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9," Branch %2ld N=%5ld IT=%2ld PAR(%2ld)=%11.3E U=",ibr,ntop,
	    nit,icp[0],rlcur[0]);	
    for (i = 0; i < ndmr; ++i) {
      fprintf(fp9,"%11.3E",u[i]);	
    }
    fprintf(fp9,"\n");	

  }

  rlm1 = rlcur[0];
  for (i = 0; i < ndim; ++i) {
    u1[i] = u[i];
  }

  for (nit = 0; nit < itnw; ++nit1) {

    iap->nit = nit + 1;
    par[icp[0]] = rlcur[0];
    (*funi)(iap, rap, ndim, u, uold, icp, par, 2, f, dfdu, dfdp);
    for (i = 0; i < ndim; ++i) {
      ARRAY2D(aa, i, ndim) = dfdp[(icp[0]) * ndim + i];
      rhs[i] = -f[i];
      for (k = 0; k < ndim; ++k) {
	ARRAY2D(aa, i, k) = dfdu[k * ndim + i];
      }
    }
    for (k = 0; k < ndim; ++k) {
      ARRAY2D(aa, ndim, k) = thu[k] * udot[k];
    }
    ARRAY2D(aa, ndim, ndim) = thl[icp[0]] * rldot[0];
    ss = 0.;
    for (i = 0; i < ndim; ++i) {
      ss += thu[i] * (u[i] - u1[i]) * udot[i];
    }
    rhs[ndim] = -ss - thl[icp[0]] * (rlcur[0] - rlm1) * rldot[0];

    /* Use Gauss elimination with pivoting to solve the linearized system 
       : */

    if (iid >= 5) {
      integer tmp = ndim + 1;
      wrjac(iap, &tmp, m1aaloc, aa, rhs);
    }
    ge(ndim + 1, *m1aaloc, aa, 1, ndim + 1, du, 
       ndim + 1, rhs, &det);
    rap->det = det;
    drlm = du[ndim];

    /* Add the Newton increments : */

    for (i = 0; i < ndim; ++i) {
      u[i] += du[i];
    }
    rlcur[0] += drlm;
    dumx = 0.;
    umx = 0.;
    for (i = 0; i < ndim; ++i) {
      adu = fabs(du[i]);
      if (adu > dumx) {
	dumx = adu;
      }
      au = fabs(u[i]);
      if (au > umx) {
	umx = au;
      }
    }

    if (iid >= 2 && iap->mynode == 0) {
      fprintf(fp9," Branch %2ld N=%5ld IT=%2ld PAR(%2ld)=%11.3E U=",ibr,ntop,
	      nit+1,icp[0],rlcur[0]);	
      for (i = 0; i < ndmr; ++i) {
	fprintf(fp9,"%11.3E",u[i]);	
      }
      fprintf(fp9,"\n");	
    }

    /* Check whether relative error has reached user-supplied tolerance : 
     */

    rdrlm = fabs(drlm) / (fabs(rlcur[0]) + 1.);
    rdumx = dumx / (umx + 1.);
    if (rdrlm < epsl && rdumx < epsu) {
      return 0;
    }
    /* L3: */
  }

  /* Maximum number of iterations reached. Reduce stepsize and try again. */

  if (iads == 0 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li NOTE:No convergence when switching branches with fixed step size\n",ibr,ntop);	
  }
  if (iads == 0) {
    goto L5;
  }

  mxt = itnw;
  iap->nit = mxt;
  adptds(iap, rap, rds);
  if (fabs(*rds) < dsmin) {
    goto L4;
  }
  rlcur[0] = rlold[0] + *rds * rldot[0];
  for (i = 0; i < ndim; ++i) {
    u[i] = uold[i] + *rds * udot[i];
  }
  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li NOTE:Retrying step\n",ibr,ntop);	
  }
  goto L2;

/* Minimum stepsize reached. */

 L4:
  if (iap->mynode == 0) {
    fprintf(fp9,"%4li%6li NOTE:No convergence when switching branches with minimum step size\n",ibr,ntop);	
  }
 L5:
  rlcur[0] = rlold[0];
  par[icp[0]] = rlcur[0];
  for (i = 0; i < ndim; ++i) {
    u[i] = uold[i];
  }
  istop = 1;
  iap->istop = istop;

  free(u1);

  return 0;
	
} /* swprc_ */

	  
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*                    Output (Algebraic Problems) */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ---- */
/* Subroutine */ int 
sthd(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *thl, doublereal *thu)
{
  /* Local variables */
  integer ndim, ncol, mxbf, nicp;
  doublereal epsl;
  integer nfpr, iplt, nint;
  doublereal epsu, epss;
  integer jtmp, itmx, itnw, nwtn, ntst, nuzr, i;
  doublereal dsmin, dsmax, a0, a1, ds, rl0, rl1;
  integer iad, jac, nbc, iid, ilp, ips, isp, irs, npr, isw, nmx;

  /* Write the values of the user defined parameters on unit 7. */
  /* This identifying information is preceded by a '   0' on each line. */
  /* The first line in the file contains the (generally) user-supplied */
  /* limits of the bifurcation diagram, viz. RL0,RL1,A0 and A1. */
  /* These are often convenient for an initial plot of the diagram. */
        
  /* Parameter adjustments */
  /*--thu;*/
  /*--thl;*/
  /*--icp;*/
  /*--par;*/
  
  ndim = iap->ndim;
  ips = iap->ips;
  irs = iap->irs;
  ilp = iap->ilp;
  ntst = iap->ntst;
  ncol = iap->ncol;
  iad = iap->iad;
  isp = iap->isp;
  isw = iap->isw;
  iplt = iap->iplt;
  nbc = iap->nbc;
  nint = iap->nint;
  nmx = iap->nmx;
  nuzr = iap->nuzr;
  npr = iap->npr;
  mxbf = iap->mxbf;
  iid = iap->iid;
  itmx = iap->itmx;
  itnw = iap->itnw;
  nwtn = iap->nwtn;
  jac = iap->jac;
  nfpr = iap->nfpr;
  nicp = iap->nicp;
  
  ds = rap->ds;
  dsmin = rap->dsmin;
  dsmax = rap->dsmax;
  rl0 = rap->rl0;
  rl1 = rap->rl1;
  a0 = rap->a0;
  a1 = rap->a1;
  epsl = rap->epsl;
  epsu = rap->epsu;
  epss = rap->epss;
  
  if (iap->mynode > 0) {
    return 0;
  }
  

  fprintf(fp7,"   0 %12.4E%12.4E%12.4E%12.4E\n",rl0,rl1,a0,a1);
  fprintf(fp7,"   0   EPSL=%11.4E  EPSU =%11.4E  EPSS =%11.4E\n",epsl,epsu,epss);
  fprintf(fp7,"   0   DS  =%11.4E  DSMIN=%11.4E  DSMAX=%11.4E\n",ds,dsmin,dsmax);
  fprintf(fp7,"   0   NDIM=%4li   IPS =%4li   IRS =%4li   ILP =%4li\n",ndim,ips,irs,ilp); 
  fprintf(fp7,"   0   NTST=%4li   NCOL=%4li   IAD =%4li   ISP =%4li\n",ntst,ncol,iad,isp);
  fprintf(fp7,"   0   ISW =%4li   IPLT=%4li   NBC =%4li   NINT=%4li\n",isw,iplt,nbc,nint);
  fprintf(fp7,"   0   NMX=%5ld   NPR =%4li   MXBF=%4li   IID =%4li\n",nmx,npr,mxbf,iid);
  fprintf(fp7,"   0   ITMX=%4li   ITNW=%4li   NWTN=%4li   JAC=%4li  NUZR=%4li\n",itmx,itnw,nwtn,jac,nuzr);

  if (nicp == 1) {
    jtmp = NPARX;
    fprintf(fp7,"   0   User-specified parameter:       ");
    for (i = 0; i < nicp; ++i) {
      fprintf(fp7,"%4li",icp[jtmp + i]);
    }
    fprintf(fp7,"\n");
  } else {
    jtmp = NPARX;
    fprintf(fp7,"   0   User-specified parameters:      ");
    for (i = 0; i < nicp; ++i) {
      fprintf(fp7,"%4li",icp[jtmp + i]);
    }
    fprintf(fp7,"\n");
  }

  if (nfpr == 1) {
    fprintf(fp7,"   0   Active continuation parameter:  ");
	
    for (i = 0; i < nfpr; ++i) {
      fprintf(fp7,"%4li",icp[i]);
    }
    fprintf(fp7,"\n");
  } else {
    fprintf(fp7,"   0   Active continuation parameters:  ");
    for (i = 0; i < nfpr; ++i) {
      fprintf(fp7,"%4li",icp[i]);
    }
    fprintf(fp7,"\n");
  }


  return 0;
} /* sthd_ */


/*     ---------- ------ */
/* Subroutine */ int 
headng(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, integer iunit, integer *n1, integer *n2)
{
  /* Local variables */
  integer iplt, itmp, i, j;
  char col[9][14+1];
  integer ndm, ips, itp;
  integer isw;

  /* Prints headings above columns on unit 6 and 7. */

  /* Parameter adjustments */
  /*--icp;*/
  /*--par;*/
    

    
  ips = iap->ips;
  isw = iap->isw;
  iplt = iap->iplt;
  ndm = iap->ndm;
  itp = iap->itp;

  /*initialize strings*/
  for (i = 0; i < 9; ++i) {
    sprintf(col[i], "              ");
  }

  if (iap->mynode == 0) {
    if (iunit == 6) {
      printf(" \n");
    }
    if (iunit == 7) {
      fprintf(fp7,"   0\n");
    }
    if (iunit == 9) {
      fprintf(fp9," \n");	
    }
  }

  j = 0;
  for (i = 0; i < *n1; ++i) {
    ++j;
    if ((double) j == (double)2.) {
      j = j + 1 + *n2;
    }   
    if (icp[i] > 9) {
      sprintf(col[j-1],"   PAR(%ld)    ",icp[i]);
    } else {
      sprintf(col[j-1],"   PAR(%ld)     ",icp[i]);
    }
	
  }

  if (iplt > ndm && iplt <= ndm << 1) {
    sprintf(col[1]," INTEGRAL U(%ld)",iplt-ndm);
  } else if (iplt > ndm << 1 && iplt <= ndm * 3) {
    sprintf(col[1]," L2-NORM U(%ld) ",iplt - (ndm * 2));
  } else if (iplt > 0 && iplt <= ndm) {
    if (abs(ips) <= 1 || ips == 5) {
      sprintf(col[1],"     U(%ld)     ",-iplt);
    } else {
      sprintf(col[1],"   MAX U(%ld)   ",iplt);
    }
  } else if (iplt < 0 && iplt >= -ndm) {
    if (abs(ips) <= 1 || ips == 5) {
      sprintf(col[1],"     U(%ld)     ",-iplt);
    } else {
      sprintf(col[1],"   MIN U(%ld)   ",-iplt);
    }
  } else {
    sprintf(col[1],"   L2-NORM    ");
  }

  if (*n2 > 0) {
    for (i = 0; i < *n2; ++i) {
      sprintf(col[i+2],"     U(%ld)     ",i + 1);
      itmp = i;
    }
    if ((ips >= 2 && ips <= 4) || (ips >= 6 && ips <= 9) || (ips >= 12 && ips <= 17)) {
      for (i = 3; i <= *n2 + 2; ++i) {
	col[i-1][3]='M';
	col[i-1][4]='A';
	col[i-1][5]='X';
      }
    }
  }

  for (i = 0; i < *n1 + *n2 + 1; ++i) {
    if (strcmp(col[i],"   PAR(10)    ") == 0 && ips > 0 && ips != 4) {
      sprintf(col[i],"    PERIOD    ");
    } else if (strcmp(col[i],"   PAR(9)    ") == 0 && (ips == 5 || ips == 15)) {
      sprintf(col[i],"     FOPT     ");
    } else if (strcmp(col[i],"   PAR(13)    ") == 0 && (ips == 14 || ips == 16)) {
      sprintf(col[i],"     TIME     ");
    }
  }

  if (iap->mynode == 0) {
    if (iunit == 6) {
      printf("  BR    PT  TY LAB ");
      for (i = 0; i < *n1 + *n2 + 1; ++i) {
	printf("%s",col[i]);
      }
      printf("\n");
      fflush(stdout);
    } else if (iunit == 7) {
      fprintf(fp7,"   0    PT  TY LAB ");
      for (i = 0; i < *n1 + *n2 + 1; ++i) {
	fprintf(fp7,"%s",col[i]);
      }
      fprintf(fp7,"\n");
    } else if (iunit == 9) {
      fprintf(fp9,"  BR    PT  TY LAB ");	
      for (i = 0; i < *n1 + *n2 + 1; ++i) {
	fprintf(fp9,"%s",col[i]);
      }
      fprintf(fp9,"\n");	

    }
  }

  return 0;
} /* headng_ */


/*     ---------- ------ */
/* Subroutine */ int 
stplae(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *rlcur, doublereal *u)
{

    /* Local variables */
  integer labw, ndim, nins, iplt, ntot, i;
  doublereal a0, a1;
  integer istop, itpst, ntots;
  doublereal ss;

  doublereal rl0, rl1;
  integer iab, lab, ibr, ndm;
  doublereal amp;
  integer ips, itp, npr, isw, nmx;
  int iflag=0;




/* Stores the bifurcation diagram on unit 7 (Algebraic Problems). */
/* Every line written contains, in order, the following: */

/*  IBR    : The label of the branch. */
/*  NTOT   : The index of the point on the branch. */
/*           (Points are numbered consecutively along a branch). */
/*           If IPS=1 or -1, then the sign of NTOT indicates stability : 
*/
/*            - = stable , + = unstable, unknown, or not relevant. */
/*  ITP    : An integer indicating the type of point : */

/*             1  (BP)  :   Branch point. */
/*             2  (LP)  :   Fold. */
/*             3  (HB)  :   Hopf bifurcation point. */
/*             4  (  )  :   Output point (Every NPR steps along branch). 
*/
/*            -4  (UZ)  :   Output point (Zero of user function). */
/*             9  (EP)  :   End point of branch, normal termination. */
/*            -9  (MX)  :   End point of branch, abnormal termination. */

/*  LAB        : The label of a special point. */
/*  PAR(ICP(1)): The principal parameter. */
/*  A          : The L2-norm of the solution vector, or other measure of 
*/
/*               the solution (see the user-supplied parameter IPLT). */
/*  U          : The first few components of the solution vector. */
/*  PAR(ICP(*)): Further free parameters (if any). */


  /* Parameter adjustments */
  /*--u;*/
  /*--rlcur;*/
  /*--icp;*/
    
    

  
  ndim = iap->ndim;
  ips = iap->ips;
  isw = iap->isw;
  iplt = iap->iplt;
  nmx = iap->nmx;
  npr = iap->npr;
  ndm = iap->ndm;
  itp = iap->itp;
  itpst = iap->itpst;
  ibr = iap->ibr;

  rl0 = rap->rl0;
  rl1 = rap->rl1;
  a0 = rap->a0;
  a1 = rap->a1;

  ntot = iap->ntot;
  ++ntot;
  iap->ntot = ntot;

  pvlsae(iap, rap, u, par);

/* ITP is set to 4 every NPR steps along a branch, and the entire */
/* solution is written on unit 8. */

  if (npr != 0) {
    if (ntot % npr == 0 && itp % 10 == 0) {
      itp = itpst * 10 + 4;
    }
    iap->itp = itp;
  }

  /* CHECK WHETHER LIMITS OF THE BIFURCATION DIAGRAM HAVE BEEN REACHED : */

  iab = abs(iplt);

  if (iab <= ndim && iab > 0) {
    amp = u[-1 + iab];
  } else if (iplt > ndim && iplt <= (ndim*2)) {
    amp = u[-1 + iplt - ndim];
  } else if (iplt > (ndim*2) && iplt <= ndim * 3) {
    amp = u[-1 + iplt - (ndim*2)];
  } else {
    ss = 0.;
    for (i = 0; i < ndm; ++i) {
      ss += u[i] * u[i];
    }
    amp = sqrt(ss);
  }
  rap->amp = amp;
  byeauto_(&iflag); 
  istop = iap->istop;
  if (istop == 1) {
    /*        Maximum number of iterations reached somewhere. */
    itp = -9 - itpst * 10;
    iap->itp = itp;
  } else if (istop == -1) {
    /*        ** UZR endpoint */
    itp = itpst * 10 + 9;
    iap->itp = itp;
  } else {
    if (rlcur[0] < rl0 || rlcur[0] > rl1 || amp < a0 || amp > a1 || ntot 
	== nmx|| iflag==1) {
      istop = 1;
      iap->istop = istop;
      itp = itpst * 10 + 9;
      iap->itp = itp;
    }
  }

  labw = 0;
  if (itp % 10 != 0) {
    lab = iap->lab;
    ++lab;
    iap->lab = lab;
    labw = lab;
  }

  /* Determine stability and print output on units 6 and 7. */

  ntots = ntot;
  nins = iap->nins;
  if (abs(ips) == 1 && abs(isw) != 2 && ntot > 1) {
    if (nins == ndim) {
      ntots = -ntot;
    }
  }
  addbif(iap,rap,ntots,iap->ibr,par,icp,labw,&amp, u, u, u, u);
  /* addbif_(ibr, ntot, itp, labw, 
   npar,amp, u, u, u, u,ndim);  */
  wrline(iap, rap, par, icp, &icp[NPARX], &ibr, &ntots, 
	 &labw, &amp, u);
  

  /* Write restart information for multi-parameter analysis : */

  if (labw != 0) {
    wrtsp8(iap, rap, par, icp, &labw, rlcur, u);
  }

  return 0;
} /* stplae_ */


/*     ---------- ------ */
/* Subroutine */ int 
wrline(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, integer *icu, integer *ibr, integer *ntot, integer *lab, doublereal *vaxis, doublereal *u)
{

  /* Local variables */
  integer nicp, nfpr, mtot, i;
  char atype[3];
  integer n1, n2;

  integer nt, ndm, ips, itp, isw,lb;
  lb=*lab;
  if((restart_flag==1)&&(lb!=0)){
    
    restart_flag=0;
    RestartLabel=*lab;
  }







/* Write one line of output on unit 6 and 7. */


  /* Parameter adjustments */
  /*--u;*/
  /*--icu;*/
  /*--icp;*/
  /*--par;*/
    

    
  ips = iap->ips;
  isw = iap->isw;
  ndm = iap->ndm;
  itp = iap->itp;
  nfpr = iap->nfpr;
  nicp = iap->nicp;

  n1 = nicp;
  n2 = ndm;
  nt = n1 + n2;

  if (n1 > 7) {
    n1 = 7;
    n2 = 0;
  } else if (nt > 7) {
    n2 = 7 - n1;
  }

  /* Write a heading above the first line. */

  if (abs(*ntot) == 1) {
    headng(iap, rap, par, icu, 6, &n1, &n2);
  }
  if (abs(*ntot) == 1) {
    headng(iap, rap, par, icu, 7, &n1, &n2);
  }
  headng(iap, rap, par, icu, 9, &n1, &n2);

  if (itp % 10 == 1) {
    strcpy(atype,"BP");
  } else if (itp % 10 == 2) {
    strcpy(atype, "LP");
  } else if (itp % 10 == 3) {
    strcpy(atype, "HB");
  } else if (itp % 10 == 4) {
    strcpy(atype, "  ");
  } else if (itp % 10 == -4) {
    strcpy(atype, "UZ");
  } else if (itp % 10 == 5) {
    strcpy(atype, "LP");
  } else if (itp % 10 == 6) {
    strcpy(atype, "BP");
  } else if (itp % 10 == 7) {
    strcpy(atype, "PD");
  } else if (itp % 10 == 8) {
    strcpy(atype, "TR");
  } else if (itp % 10 == 9) {
    strcpy(atype, "EP");
  } else if (itp % 10 == -9) {
    strcpy(atype, "MX");
  } else {
    strcpy(atype, "  ");
  }

  if (iap->mynode > 0) {
    return 0;
  }

  mtot = *ntot % 10000;
  if (n2 == 0) {
    if (itp % 10 != 0) {
      printf("%4li%6li  %c%c%4li",(*ibr),mtot,atype[0],atype[1],*(lab));
      printf("%14.6E",par[icu[0]]);
      printf("%14.6E",(*vaxis));
      for (i = 1; i < n1; ++i) {
	printf("%14.6E",par[icu[i]]);
      }
      printf("\n");
      fflush(stdout);
    }
    fprintf(fp7,"%4li%6li%4li%4li",(*ibr),mtot,itp,(*lab));
    fprintf(fp7,"%14.5E",par[icu[0]]);
    fprintf(fp7,"%14.5E",(*vaxis));
    for (i = 1; i < n1; ++i) {
      fprintf(fp7,"%14.5E",par[icu[i]]);
    }
    fprintf(fp7,"\n");
    fprintf(fp9,"%4li%6li  %c%c%4li",(*ibr),mtot,atype[0],atype[1],*(lab));	
    fprintf(fp9,"%14.6E",par[icu[0]]);
    fprintf(fp9,"%14.6E",(*vaxis));
    for (i = 1; i < n1; ++i) {
      fprintf(fp9,"%14.6E",par[icu[i]]);
    }
    fprintf(fp9,"\n");
  } else {
    if (n1 == 1) {
      if (itp % 10 != 0) {
	printf("%4li%6li  %c%c%4li",abs(*ibr),abs(mtot),atype[0],atype[1],*(lab));
	printf("%14.6E",par[icu[0]]);
	printf("%14.6E",(*vaxis));
	for (i = 0; i < n2; ++i) {
	  printf("%14.6E",u[i]);
	}
	printf("\n");
	fflush(stdout);
      }
      fprintf(fp7,"%4li%6li%4li%4li",(*ibr),mtot,itp,(*lab));
      fprintf(fp7,"%14.5E",par[icu[0]]);
      fprintf(fp7,"%14.5E",(*vaxis));
      for (i = 0; i < n2; ++i) {
	fprintf(fp7,"%14.5E",u[i]);
      }
      fprintf(fp7,"\n");
      fprintf(fp9,"%4li%6li  %c%c%4li",(*ibr),mtot,atype[0],atype[1],*(lab));
      fprintf(fp9,"%14.6E",par[icu[0]]);
      fprintf(fp9,"%14.6E",(*vaxis));
      for (i = 0; i < n2; ++i) {
	fprintf(fp9,"%14.6E",u[i]);
      }
      fprintf(fp9,"\n");	

    } else {
      if (itp % 10 != 0) {
	printf("%4li%6li  %c%c%4li",abs(*ibr),abs(mtot),atype[0],atype[1],*(lab));
	printf("%14.6E",par[icu[0]]);
	printf("%14.6E",(*vaxis));
	for (i = 0; i < n2; ++i) {
	  printf("%14.6E",u[i]);
	}
	for (i = 1; i < n1; ++i) {
	  printf("%14.6E",par[icu[i]]);
	}
	printf("\n");
	fflush(stdout);      
      }
      fprintf(fp7,"%4li%6li%4li%4li",(*ibr),mtot,itp,(*lab));
      fprintf(fp7,"%14.5E",par[icu[0]]);
      fprintf(fp7,"%14.5E",(*vaxis));
      for (i = 0; i < n2; ++i) {
	fprintf(fp7,"%14.5E",u[i]);
      }
      for (i = 1; i < n1; ++i) {
	fprintf(fp7,"%14.5E",par[icu[i]]);
      }
      fprintf(fp7,"\n");
      fprintf(fp9,"%4li%6li  %c%c%4li",(*ibr),mtot,atype[0],atype[1],*(lab));
      fprintf(fp9,"%14.6E",par[icu[0]]);
      fprintf(fp9,"%14.6E",(*vaxis));
      for (i = 0; i < n2; ++i) {
	fprintf(fp9,"%14.6E",u[i]);
      }
      for (i = 1; i < n1; ++i) {
	fprintf(fp9,"%14.6E",par[icu[i]]);
      }
      fprintf(fp9,"\n");	
    }
  }


  return 0;
} /* wrline_ */

/*     ---------- ------ */
/* Subroutine */ int 
wrtsp8(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, integer *lab, doublereal *rlcur, doublereal *u)
{
  /* Local variables */
  integer ndim, nfpr, ntpl, jtmp, mtot, ntot, i;
  doublereal t;
  integer nrowpr, ibr;
  doublereal amp;
  integer nar, itp, isw;

  /*   static FILE *fp8=NULL; */

  if(fp8_is_open==0){
    fp8 = fopen(fort8,"w");
    fp8_is_open=1;
    if(fp8 == NULL) {
      fprintf(stderr,"Error:  Could not open fort.8\n");
      exit(1);
    }
  }




  /* Write restart information on singular points, plotting points, etc., */
  /* on unit 8. */


  /* Parameter adjustments */
  /*--u;*/
  /*--rlcur;*/
  /*--icp;*/
  /*--par;*/
    
    

    
  ndim = iap->ndim;
  isw = iap->isw;
  itp = iap->itp;
  ibr = iap->ibr;
  nfpr = iap->nfpr;
  ntot = iap->ntot;

  ntpl = 1;
  nar = ndim + 1;
  jtmp = NPARX;
  nrowpr = ndim / 7 + 1 + (jtmp - 1) / 7 + 1;
  par[icp[0]] = rlcur[0];
  t = 0.;
  amp = 0.;
  rap->amp = amp;
  if (iap->mynode > 0) {
    return 0;
  }

  mtot = ntot % 10000;
  fprintf(fp8,"%5ld",ibr);
  fprintf(fp8,"%5ld",mtot);
  fprintf(fp8,"%5ld",itp);
  fprintf(fp8,"%5ld",(*lab));
  fprintf(fp8,"%5ld",nfpr);
  fprintf(fp8,"%5ld",isw);
  fprintf(fp8,"%5ld",ntpl);
  fprintf(fp8,"%5ld",nar);
  fprintf(fp8,"%5ld",nrowpr);
  fprintf(fp8,"%5d",0);
  fprintf(fp8,"%5d",0);
  fprintf(fp8,"%5d\n",NPARX);
  fprintf(fp8,"    %19.10E",t);
  for (i = 0; i < ndim; ++i) {
    if((i>0)&&((i+1)%7==0))
      fprintf(fp8,"\n    ");
    fprintf(fp8,"%19.10E",u[i]);
  }
  fprintf(fp8,"\n");
  for (i = 0; i < NPARX; ++i) {
    if(i==0)
      fprintf(fp8,"    ");
    if((i>0)&&(i%7==0))
      fprintf(fp8,"\n    ");
    fprintf(fp8,"%19.10E",par[i]);
  }
  fprintf(fp8,"\n");    
  fflush(fp8);

  return 0;
} /* wrtsp8_ */


/*     ---------- ------ */
/* Subroutine */ int 
wrjac(iap_type *iap, integer *n, integer *m1aaloc, doublereal *aa, doublereal *rhs)
{
  /* System generated locals */
  integer aa_dim1;

  /* Local variables */
  integer i, j;

  /* Parameter adjustments */
  /*--rhs;*/
  aa_dim1 = *m1aaloc;
    
  if (iap->mynode > 0) {
    return 0;
  }
  fprintf(fp9," Residual vector :\n");	

  for (i = 0; i < *n; ++i) {
    fprintf(fp9," %10.3E",rhs[i]);	
  }
  fprintf(fp9,"\n");	
  fprintf(fp9," Jacobian matrix :\n");	
  for (i = 0; i < *n; ++i) {
    for (j = 0; j < *n; ++j) {
      fprintf(fp9," %10.3E",ARRAY2D(aa, i, j));	
    }
    fprintf(fp9,"\n");	

  }


  return 0;
} /* wrjac_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*                    Mesh and Weight Generation */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- --- */
/* Subroutine */ int 
msh(const iap_type *iap, const rap_type *rap, doublereal *tm)
{

    /* Local variables */
  integer ntst, j;
  doublereal dt;




/* Generates a uniform mesh on [0,1]. */


  /* Parameter adjustments */
  /*--tm;*/
    

  
  ntst = iap->ntst;

  tm[0] = 0.;
  dt = 1. / ntst;
  for (j = 0; j < ntst; ++j) {
    tm[j + 1] = (j + 1) * dt;
  }

  return 0;
} /* msh_ */


/*     ---------- ------ */
/* Subroutine */ int 
genwts(const integer ncol, const integer n1, doublereal *wt, doublereal *wp)
{
  /* System generated locals */
  integer wt_dim1, wp_dim1;

    /* Local variables */
  doublereal d;
  integer i, k, l;
  doublereal p, denom;

  integer ib, ic;
  doublereal *xm, *zm, sum;
  integer ncp1;


  xm = (doublereal *)malloc(sizeof(doublereal)*(ncol + 1));
  zm = (doublereal *)malloc(sizeof(doublereal)*(ncol));

/* Generates weights of the collocation method. The user selected */
/* number of collocation points (ncol) must be one of { 2,...,7 }. */

/* The following weights are generated : */

/*         WT : for the function value, */
/*         WP : for the first derivative, */

/* Local */

/* Generate the collocation points : */
  /* Parameter adjustments */
  wp_dim1 = n1;
  wt_dim1 = n1;
    
  cpnts(ncol, zm);

  ncp1 = ncol + 1;
  d = 1. / ncol;
  for (i = 0; i < ncp1; ++i) {
    xm[i] = i * d;
  }

  /* Generate weights : */

  for (ib = 0; ib < ncp1; ++ib) {
    denom = 1.;
    for (k = 0; k < ncp1; ++k) {
      if (k != ib) {
	denom *= xm[ib] - xm[k];
      }
    }
    for (ic = 0; ic < ncol; ++ic) {
      /* Weights for the function values : */
      p = 1.;
      for (k = 0; k < ncp1; ++k) {
	if (k != ib) {
	  p *= zm[ic] - xm[k];
	}
      }
      ARRAY2D(wt, ib, ic) = p / denom;
      /* Weights for derivatives : */
      sum = 0.;
      for (l = 0; l < ncp1; ++l) {
	if (l != ib) {
	  p = 1.;
	  for (k = 0; k < ncp1; ++k) {
	    if (k != ib && k != l) {
	      p *= zm[ic] - xm[k];
	    }
	  }
	  sum += p;
	}
      }
      ARRAY2D(wp, ib, ic) = sum / denom;
    }
  }
  free(xm);
  free(zm);

  return 0;
} /* genwts_ */


/*     ---------- ----- */
/* Subroutine */ int 
cpnts(const integer ncol, doublereal *zm)
{
    /* Local variables */
  doublereal c, r, c1, c2, c3;


  /* Generates the collocation points with respect to [0,1]. */
  if (ncol > 7) {
    fprintf(stderr,"Dimension exceeded : NCOL=%5ld  maximum=7\n",ncol);
    fprintf(stderr,"AUTO does not contain weights for NCOL > 1\n");
    fprintf(stderr,"Please reset NCOL to 7 or smaller\n");
    exit(1);
  }

  switch ((int)(ncol - 1)) {
  case 1:  goto L2;
  case 2:  goto L3;
  case 3:  goto L4;
  case 4:  goto L5;
  case 5:  goto L6;
  case 6:  goto L7;
  }

 L2:
  c = .5 / sqrt(3.);
  zm[0] = .5 - c;
  zm[1] = c + .5;
  return 0;

 L3:
  c = sqrt(.6) * .5;
  zm[0] = .5 - c;
  zm[1] = .5;
  zm[2] = c + .5;
  return 0;

 L4:
  r = .8571428571428571;
  /* Computing 2nd power */
  c = sqrt(r * r - .34285714285714286) * .5;
  c1 = sqrt(c + .42857142857142855) * .5;
  c2 = sqrt(.42857142857142855 - c) * .5;
  zm[0] = .5 - c1;
  zm[1] = .5 - c2;
  zm[2] = c2 + .5;
  zm[3] = c1 + .5;
  return 0;

 L5:
  c1 = .45308992296933198;
  c2 = .26923465505284155;
  zm[0] = .5 - c1;
  zm[1] = .5 - c2;
  zm[2] = .5;
  zm[3] = c2 + .5;
  zm[4] = c1 + .5;
  return 0;

 L6:
  c1 = .46623475710157603;
  c2 = .33060469323313224;
  c3 = .11930959304159845;
  zm[0] = .5 - c1;
  zm[1] = .5 - c2;
  zm[2] = .5 - c3;
  zm[3] = c3 + .5;
  zm[4] = c2 + .5;
  zm[5] = c1 + .5;
  return 0;

 L7:
  c1 = .4745539956171379;
  c2 = .37076559279969723;
  c3 = .20292257568869859;
  zm[0] = .5 - c1;
  zm[1] = .5 - c2;
  zm[2] = .5 - c3;
  zm[3] = .5;
  zm[4] = c3 + .5;
  zm[5] = c2 + .5;
  zm[6] = c1 + .5;
  return 0;
} /* cpnts_ */


/*     ---------- ------ */
/* Subroutine */ int 
cntdif(integer *n, doublereal *d)
{

    /* Local variables */
  integer i, k, k1;
  doublereal sc;

  /* Generates the coefficients of the central difference formula for */
  /* Nth derivative at uniformly spaced points */
  /*              0 = x  < x  < ... < x  = 1. */
  /*                   0    1          N */


  /* Parameter adjustments */
  /*--d;*/

  
  d[0] = 1.;
  if (*n == 0) {
    return 0;
  }
 
  for (i = 0; i < *n; ++i) {
    d[i + 1] = 0.;
    for (k = 0; k < i + 1; ++k) {
      k1 = i + 1 - k;
      d[k1] = d[k1 - 1] - d[k1];
    }
    d[0] = -d[0];
  }

  /* Scale to [0,1]  : */

  sc = (doublereal) pow_ii(*n, *n);
  for (i = 0; i < *n + 1; ++i) {
    d[i] = sc * d[i];
  }

  return 0;
} /* cntdif_ */


/*     ---------- ---- */
/* Subroutine */ int 
wint(const integer n, doublereal *wi)
{
  doublereal c;

  /* Generates the weights for the integration formula based on polynomial */
  /* interpolation at N equally spaced points in [0,1]. */
  
  switch ((int)(n - 2)) {
  case 1:  goto L3;
  case 2:  goto L4;
  case 3:  goto L5;
  case 4:  goto L6;
  case 5:  goto L7;
  case 6:  goto L8;
  }

 L3:
  c = .16666666666666666;
  wi[0] = c;
  wi[1] = c * 4.;
  wi[2] = c;
  return 0;

 L4:
  c = .125;
  wi[0] = c;
  wi[1] = c * 3.;
  wi[2] = wi[1];
  wi[3] = c;
  return 0;

 L5:
  c = .011111111111111112;
  wi[0] = c * 7.;
  wi[1] = c * 32.;
  wi[2] = c * 12.;
  wi[3] = wi[1];
  wi[4] = wi[0];
  return 0;

 L6:
  wi[0] = .065972222222222224;
  wi[1] = .26041666666666669;
  wi[2] = .1736111111111111;
  wi[3] = wi[2];
  wi[4] = wi[1];
  wi[5] = wi[0];
  return 0;

 L7:
  wi[0] = .04880952380952381;
  wi[1] = .25714285714285712;
  wi[2] = .03214285714285714;
  wi[3] = .32380952380952382;
  wi[4] = wi[2];
  wi[5] = wi[1];
  wi[6] = wi[0];
  return 0;

 L8:
  wi[0] = .043460648148148151;
  wi[1] = .20700231481481482;
  wi[2] = .076562500000000006;
  wi[3] = .17297453703703702;
  wi[4] = wi[3];
  wi[5] = wi[2];
  wi[6] = wi[1];
  wi[7] = wi[0];

  return 0;
} /* wint_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*          Stepsize and Mesh Adaption */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ------ */
/* Subroutine */ int 
adptds(iap_type *iap, rap_type *rap, doublereal *rds)
{
  /* Local variables */
  doublereal ards;
  integer ntop, itnw, ntot;
  doublereal dsmax;
  integer n1, iid, ibr, nit;

/* The stepsize along the branch of solutions is adapted depending on the 
*/
/* number of Newton iterations in the previous step (called if IADS > 0). 
*/



    
  dsmax = rap->dsmax;
  iid = iap->iid;
  itnw = iap->itnw;
  ibr = iap->ibr;
  nit = iap->nit;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  if (itnw <= 3) {
    itnw = 3;
    n1 = 2;
  } else {
    n1 = itnw / 2;
  }

  if (nit <= 1) {
    *rds *= 2.;
  } else if (nit == 2) {
    *rds *= (double)1.5;
  } else if (nit > 2 && nit <= n1) {
    *rds *= (double)1.1;
  } else if (nit >= itnw) {
    *rds *= .5;
  }

  ards = fabs(*rds);
  if (ards > dsmax) {
    *rds = *rds * dsmax / ards;
  }

  fprintf(fp9,"%4li%6li        Iterations     %3li\n",abs(ibr),ntop - 1,nit);	
  fprintf(fp9,"%4li%6li        Stepsize      %14.6E\n",abs(ibr),ntop,(*rds));	

  return 0;
} /* adptds_ */


/*     ---------- ----- */
/* Subroutine */ int 
adapt(iap_type *iap, rap_type *rap, integer *nold, integer *ncold, integer *nnew, integer *ncnew, doublereal *tm, doublereal *dtm, integer *ndxloc, doublereal *ups, doublereal *vps)
{
  /* System generated locals */
  integer ups_dim1, vps_dim1;

    /* Local variables */
  integer ndim, iper;
  integer i, j, noldp1, nnewp1;

  integer nrwnew;
  integer ips, isw;
    
  doublereal *tint, *uint;
  doublereal *tm2;
  integer *itm;


  uint = (doublereal *)malloc(sizeof(doublereal)*(*ndxloc)*(iap->ndim * iap->ncol));
  tint = (doublereal *)malloc(sizeof(doublereal)*(*ndxloc));
  tm2 = (doublereal *)malloc(sizeof(doublereal)*(*ndxloc));
  itm = (integer *)malloc(sizeof(integer)*(*ndxloc));

  /* Adapts the distribution of the mesh points so that the increase of the */
  /* monotone function EQDF becomes approximately equidistributed over the */
  /* intervals. The functions UPS and VPS are interpolated on new mesh. */

/* Local */

  /* Parameter adjustments */
  /*--tm;*/
  /*--dtm;*/
  vps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
    
  ndim = iap->ndim;
  ips = iap->ips;
  isw = iap->isw;

  noldp1 = *nold + 1;
  nnewp1 = *nnew + 1;
  nrwnew = ndim * *ncnew;

  for (j = 0; j < (*ndxloc); ++j) {
    for (i = 0; i < (iap->ndim * iap->ncol); ++i) {
      uint[j + i * (*ndxloc)] = 0.;
    }
  }

  /* For periodic boundary conditions extrapolate by periodicity. */

  if (ips == 2 && abs(isw) != 2) {
    iper = 1;
  } else {
    iper = 0;
  }

  /* Generate the new mesh : */

  newmsh(iap, rap, ndxloc, ups, nold, 
	 ncold, tm, dtm, nnew, tint, &iper);

  /* Replace UPS by its interpolant on the new mesh : */

  interp(iap, rap, &ndim, &noldp1, ncold, tm, ndxloc, ups, &nnewp1, ncnew, tint, uint, tm2, itm);
  for (j = 0; j < nnewp1; ++j) {
    for (i = 0; i < nrwnew; ++i) {
      ARRAY2D(ups, j, i) = uint[j + i * (*ndxloc)];
    }
  }

  /* Replace VPS by its interpolant on the new mesh : */

  interp(iap, rap, &ndim, &noldp1, ncold, tm, ndxloc, 
	 vps, &nnewp1, ncnew, tint, uint, tm2, itm);
  for (j = 0; j < nnewp1; ++j) {
    for (i = 0; i < nrwnew; ++i) {
      ARRAY2D(vps, j, i) = uint[j + i * (*ndxloc)];
    }
  }

  /* Replace old mesh : */

  tm[0] = 0.;
  for (j = 0; j < *nnew; ++j) {
    dtm[j] = tint[j+1] - tint[j];
    tm[j + 1] = tint[j+1];
  }

  free(uint);
  free(tint);
  free(tm2);
  free(itm);

  return 0;
} /* adapt_ */


/*     ---------- ------ */
/* Subroutine */ int 
interp(iap_type *iap, rap_type *rap, integer *ndim, integer *n, integer *nc, doublereal *tm, integer *ndxloc, doublereal *ups, integer *n1, integer *nc1, doublereal *tm1, doublereal *ups1, doublereal *tm2, integer *itm1)
{
  /* System generated locals */
  integer ups_dim1, ups1_dim1;

  /* Local variables */

  doublereal d;
  integer i, j, k, l;
  doublereal *w, *x, z__;
  integer j1, k1, l1;
  doublereal ri;

  integer n1m1, ncp1;

  w = (doublereal *)malloc(sizeof(doublereal)*(*nc+1));
  x = (doublereal *)malloc(sizeof(doublereal)*(*nc+1));


/* Finds interpolant (TM(.) , UPS(.) ) on new mesh TM1. */

/* Local */

  /* Parameter adjustments */
  /*--tm;*/
  /*--tm1;*/
  /*--tm2;*/
  /*--itm1;*/
  ups1_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
    
  ncp1 = *nc + 1;
  n1m1 = *n1 - 1;

  for (i = 0; i < *nc1; ++i) {
    ri = (doublereal) (i);
    d = ri / *nc1;
    for (j1 = 0; j1 < n1m1; ++j1) {
      tm2[j1] = tm1[j1] + d * (tm1[j1 + 1] - tm1[j1]);
    }
    ordr(iap, rap, n, tm, &n1m1, tm2, itm1);
    for (j1 = 0; j1 <n1m1; ++j1) {
      j = itm1[j1];
      z__ = tm2[j1];
      d = (tm[j] - tm[-1 + j]) / *nc;
      for (l = 0; l < ncp1; ++l) {
	x[l] = tm[-1 + j] + l * d;
      }
      intwts(iap, rap, &ncp1, &z__, x, w);
      for (k = 0; k < *ndim; ++k) {
	k1 = i * *ndim + k;
	ARRAY2D(ups1, j1, k1) = w[ncp1 - 1] * ARRAY2D(ups, j, k);
	for (l = 0; l < *nc; ++l) {
	  l1 = k + l * *ndim;
	  ARRAY2D(ups1, j1, k1) += w[l] * ARRAY2D(ups, (j - 1), l1);
	}
      }
    }
  }

  for (i = 0; i < *ndim; ++i) {
    ARRAY2D(ups1, (*n1 - 1), i) = ARRAY2D(ups, (*n - 1), i);
  }
  free(w);
  free(x);

  return 0;
} /* interp_ */


/*     ---------- ------ */
/* Subroutine */ int 
newmsh(iap_type *iap, rap_type *rap, integer *ndxloc, doublereal *ups, integer *nold, integer *ncold, doublereal *tmold, doublereal *dtmold, integer *nnew, doublereal *tmnew, integer *iper)
{
    /* Local variables */

  integer ndim;

  integer j;
  doublereal x;
  integer j1, noldp1, nnewp1;
  doublereal dal;

  doublereal *uneq;
  doublereal *eqf;
  integer *ial;

  uneq=(doublereal *)malloc(sizeof(doublereal)*(*nnew + 1));
  eqf=(doublereal *)malloc(sizeof(doublereal)*(*nold + 1));
  ial=(integer *)malloc(sizeof(integer)*(*nnew + 1));

  /* Redistributes the mesh according to the function EQDF. */

/* Local */

  /* Parameter adjustments */
  /*--tmnew;*/
  /*--dtmold;*/
  /*--tmold;*/
    

  
  ndim = iap->ndim;

  /* Put the values of the monotonely increasing function EQDF in EQF. */

  eqdf(iap, rap, nold, &ndim, ncold, dtmold, ndxloc, ups, eqf, 
       iper);

  /* Uniformly divide the range of EQDF : */

  noldp1 = *nold + 1;
  nnewp1 = *nnew + 1;
  dal = eqf[noldp1 - 1] / *nnew;
  for (j = 0; j < nnewp1; ++j) {
    uneq[j] = (j) * dal;
  }

  ordr(iap, rap, &noldp1, eqf, &nnewp1, uneq, ial);

/* Generate the new mesh in TMNEW : */

  for (j1 = 0; j1 < nnewp1; ++j1) {
    j = ial[j1];
    x = (uneq[j1] - eqf[j - 1]) / (eqf[j] - eqf[j - 1]);
    tmnew[j1] = (1. - x) * tmold[-1 + j] + x * tmold[j];
  }

  free(uneq);
  free(eqf);
  free(ial);
  return 0;
} /* newmsh_ */


/*     ---------- ---- */
/* Subroutine */ int 
ordr(iap_type *iap, rap_type *rap, integer *n, doublereal *tm, integer *n1, doublereal *tm1, integer *itm1)
{

    /* Local variables */
  integer j, k0, j1, k1=0;


  

  /* TM and TM1 are two ascending arrays with values in [0,1]. On exit the 
*/
/* value of ITM1( i ) specifies the index of the TM-interval in which */
/* TM1(i) lies. */


  /* Parameter adjustments */
  /*--tm;*/
  /*--itm1;*/
  /*--tm1;*/

  
  k0 = 2;
  for (j1 = 0; j1 < *n1; ++j1) {
    for (j = k0; j <= *n; ++j) {
      k1 = j;
      if (tm1[j1] < tm[-1 + j]) {
	goto L1;
      }
    }
  L1:
    itm1[j1] = k1 - 1;
    k0 = k1;
  }

  return 0;
} /* ordr_ */


/*     ---------- ------ */
/* Subroutine */ int 
intwts(iap_type *iap, rap_type *rap, integer *n, doublereal *z__, doublereal *x, doublereal *wts)
{

    /* Local variables */
  integer k;
  doublereal p, denom;
  integer ib;


  

  /* Generates weights for Lagrange interpolation. */


  /* Parameter adjustments */
  /*--wts;*/
  /*--x;*/

    
  for (ib = 0; ib < *n; ++ib) {
    p = 1.;
    denom = 1.;
    for (k = 0; k < *n; ++k) {
      if (k != ib) {
	p *= *z__ - x[k];
	denom *= x[ib] - x[k];
      }
    }
    wts[ib] = p / denom;
  }

  return 0;
} /* intwts_ */


/*     ---------- ---- */
/* Subroutine */ int 
eqdf(iap_type *iap, rap_type *rap, integer *ntst, integer *ndim, integer *ncol, doublereal *dtm, integer *ndxloc, doublereal *ups, doublereal *eqf, integer *iper)
{
  /* System generated locals */
  integer ups_dim1;

  /* Local variables */
  doublereal dtav, e;
  integer i, j, k;
  logical small;
  integer k1;
  doublereal *hd, sc, *wh;

  integer jp1;
  doublereal pwr;

  hd = (doublereal *)malloc(sizeof(doublereal)*(*ntst + 1)*(*ndim * *ncol));
  wh = (doublereal *)malloc(sizeof(doublereal)*(*ncol + 1));

  /* Compute approximation to NCOL-th derivative : */
  /* Parameter adjustments */
  /*--dtm;*/
  /*--eqf;*/
  ups_dim1 = *ndxloc;

  cntdif(ncol, wh);

  small = TRUE_;
  for (j = 0; j < *ntst; ++j) {
    jp1 = j + 1;
    sc = 1. / pow_di(&dtm[j], ncol);
    for (i = 0; i < *ndim; ++i) {
      hd[j + i * (*ntst + 1)] = wh[*ncol] * ARRAY2D(ups, jp1, i);
      for (k = 0; k < *ncol; ++k) {
	k1 = i + k * *ndim;
	hd[j + i * (*ntst + 1)] += wh[k] * ARRAY2D(ups, j, k1);
      }
      hd[j + i * (*ntst + 1)] = sc * hd[j + i * (*ntst + 1)];
      if (fabs(hd[j + i * (*ntst + 1)]) > HMACH) {
	small = FALSE_;
      }
    }
  }

  /* Take care of "small derivative" case. */

  if (small) {
    for (i = 0; i < *ntst + 1; ++i) {
      eqf[i] = (doublereal) (i);
    }
    free(hd);
    free(wh);
    return 0;
  }

  if (*iper == 1) {
    /*        *Extend by periodicity : */
    for (i = 0; i < *ndim; ++i) {
      hd[*ntst + i * (*ntst + 1)] = hd[i * (*ntst + 1)];
    }
    dtm[*ntst] = dtm[0];
  } else {
    /*        *Extend by extrapolation : */
    for (i = 0; i < *ndim; ++i) {
      hd[*ntst + i * (*ntst + 1)] = hd[(*ntst -1) + i * (*ntst + 1)] * 2 - hd[(*ntst - 2) + i * (*ntst + 1)];
    }
    dtm[*ntst] = dtm[-1 + *ntst];
  }

  /* Compute approximation to (NCOL+1)-st derivative : */

  for (j = 0; j < *ntst; ++j) {
    jp1 = j + 1;
    dtav = (dtm[j] + dtm[j + 1]) * .5;
    sc = 1. / dtav;
    for (i = 0; i < *ndim; ++i) {
      hd[j + i * (*ntst + 1)] = sc * (hd[j + 1 + i * (*ntst + 1)] - hd[j + i * (*ntst + 1)]);
    }
  }

  /* Define the equidistribution function : */

  pwr = 1. / (*ncol + 1.);
  eqf[0] = 0.;
  for (j = 0; j < *ntst; ++j) {
    e = 0.;
    for (i = 0; i < *ndim; ++i) {
      doublereal tmp = fabs(hd[j + i * (*ntst + 1)]);
      e += pow_dd(&tmp, &pwr);
    }
    eqf[j + 1] = eqf[j] + dtm[j] * e;
  }
  free(hd);
  free(wh);
  return 0;
} /* eqdf_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*                    General Support Routines */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- --- */
/* Subroutine */ int 
eig(iap_type *iap, integer *ndim, integer *m1a, doublereal *a, doublecomplex *ev, integer *ier)
{
  /* System generated locals */
  integer a_dim1;

    /* Local variables */
  integer matz, ntop, ntot, i;
  doublereal *z__;

  doublereal *wi, *wr, *fv1;
  integer *iv1, ibr;

  z__ = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*(iap->ndim));
  wi  = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  wr  = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  fv1 = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  iv1 = (integer *)malloc(sizeof(integer)*(iap->ndim));

  /* This subroutine uses the EISPACK subroutine RG to compute the */
  /* eigenvalues of the general real matrix A. */
  /* NDIM is the dimension of A. */
  /* M1A is the first dimension of A as in the DIMENSION statement. */
  /* The eigenvalues are to be returned in the complex vector EV. */


/* Local */

  /* Parameter adjustments */
  /*--ev;*/
  a_dim1 = *m1a;
  
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  *ier = 0;
  matz = 0;

  rg(*m1a, *ndim, a, wr, wi, matz, z__, iv1, fv1, ier);

  for (i = 0; i < *ndim; ++i) {
    ev[i].r = wr[i], ev[i].i = wi[i];
  }

  if (*ier != 0) {
    *ier = 1;
  }
  if (*ier == 1) {
    fprintf(fp9,"%4li%6li NOTE:Error return from EISPACK routine RG\n",ibr,ntop);	
  }

  free(z__);
  free(wi );
  free(wr );
  free(fv1);
  free(iv1);
  return 0;
} /* eig_ */


/*     ---------- ---- */
/* Subroutine */ int 
nlvc(integer n, integer m, integer k, doublereal *a, doublereal *u)
{
  /* System generated locals */
  integer a_dim1;

    /* Local variables */
  integer ipiv, jpiv, i, j, l;
  doublereal p;
  integer i1, jj, kk;
  doublereal rm, sm;
  integer ip1, nmk;
  doublereal piv;
  integer jjp1;

  integer *ir,*ic;
  ir=(integer *)malloc(sizeof(integer)*(n));
  ic=(integer *)malloc(sizeof(integer)*(n));

/* Finds a null-vector of a singular matrix A. */
/* The null space of A is assumed to be K-dimensional. */

/* Parameters : */

/*     N : number of equations, */
/*     M : first dimension of A from DIMENSION statement, */
/*     K : dimension of nullspace, */
/*     A : N * N matrix of coefficients, */
/*     U : on exit U contains the null vector, */
/* IR,IC : integer arrays of dimension at least N. */



  /* Parameter adjustments */
  /*--u;*/
  /*--ir;*/
  /*--ic;*/
  a_dim1 = m;
  
  for (i = 0; i < n; ++i) {
    ic[i] = i;
    ir[i] = i;
  }

  /*   Elimination. */

  nmk = n - k;

  for (jj = 0; jj < nmk; ++jj) {
    ipiv = jj;
    jpiv = jj;
    piv = 0.;
    for (i = jj; i < n; ++i) {
      for (j = jj; j < n; ++j) {
	p = fabs(ARRAY2D(a, ir[i], ic[j]));
	if (p > piv) {
	  piv = p;
	  ipiv = i;
	  jpiv = j;
	}
      }
    }
    if (piv < RSMALL) {
      fprintf(fp9,"        NOTE:Pivot %3li < %10.3E  in NLVC : A null space may be multi-dimensional\n",jj,RSMALL);	

    }

    kk = ir[jj];
    ir[jj] = ir[ipiv];
    ir[ipiv] = kk;

    kk = ic[jj];
    ic[jj] = ic[jpiv];
    ic[jpiv] = kk;

    jjp1 = jj + 1;
    for (l = jjp1; l < n; ++l) {
      rm = ARRAY2D(a, ir[l], ic[jj]) / ARRAY2D(a, ir[jj], ic[jj]);
      if (rm != 0.) {
	for (i = jjp1; i < n; ++i) {
	  ARRAY2D(a, ir[l], ic[i]) -= rm * ARRAY2D(a, ir[jj], ic[i]);
	}
      }
    }
  }

  /*   Backsubstitution : */

  for (i = 0; i < k; ++i) {
    u[ic[-1 + n - i]] = 1.;
  }

  for (i1 = 0 ; i1 < nmk; ++i1) {
    i = nmk - i1 - 1;
    sm = 0.;
    ip1 = i + 1;
    for (j = ip1; j <  n; ++j) {
      sm += ARRAY2D(a, ir[i], ic[j]) * u[ic[j]];
    }
    u[ic[i]] = -sm / ARRAY2D(a, ir[i], ic[i]);
  }

  /* Since the pointers for ir and ic were adjusted we need to put the
       back here. */
  /*ir \+= 1;*/
  /*ic \+= 1;*/
  free(ir);
  free(ic);
  return 0;
} /* nlvc_ */


/*     ---------- ----- */
/* Subroutine */ int 
nrmlz(integer *ndim, doublereal *v)
{
    /* Local variables */
  doublereal c;
  integer i;
  doublereal ss;


  


  /* Scale the vector V so that its discrete L2-norm becomes 1. */

  /* Parameter adjustments */
  /*--v;*/

  
  ss = 0.;
  for (i = 0; i < *ndim; ++i) {
    ss += v[i] * v[i];
  }
  c = 1. / sqrt(ss);
  for (i = 0; i < *ndim; ++i) {
    v[i] *= c;
  }

  return 0;
} /* nrmlz_ */


/*     ------ --------- -------- */
doublereal 
pi(doublereal r)
{
  /* System generated locals */
  doublereal ret_val;

    





  ret_val = r * 4. * atan(1.);

  return ret_val;
} /* pi */


/*     ---------- -- */
/* Subroutine */ int 
ge(integer n, integer m1a, doublereal *a, integer nrhs, integer ndxloc, doublereal *u, integer m1f, doublereal *f, doublereal *det)
{
  /* System generated locals */
  integer a_dim1, u_dim1, f_dim1;

  /* Local variables */
  integer ipiv, jpiv, i, j, k, l;
  doublereal p;
  integer i1, jj;
  doublereal rm, sm;
  integer ip1, irh;
  doublereal piv;
  integer jjp1;

  integer *ic,*ir;
  ic=(integer *)malloc(sizeof(integer)*(n));
  ir=(integer *)malloc(sizeof(integer)*(n));


/* Solves the linear system  A U = F by Gauss elimination */
/* with complete pivoting. */

/* Parameters : */

/*   N   : number of equations, */
/*   M1A : first dimension of A from DIMENSION statement, */
/*   A   : N * N matrix of coefficients, */
/*   NRHS: 0   if no right hand sides (determinant only), */
/*         >0   if there are NRHS right hand sides, */
/*   ndxloc : first dimension of U from DIMENSION statement, */
/*   U   : on exit U contains the solution vector(s), */
/*   M1F : first dimension of F from DIMENSION statement, */
/*   F   : right hand side vector(s), */
/*  IR,IC: integer vectors of dimension at least N. */

/* The input matrix A is overwritten. */


  /* Parameter adjustments */
  /*--ir;*/
  /*--ic;*/
  a_dim1 = m1a;
  u_dim1 = ndxloc;
  f_dim1 = m1f;
  
  for (i = 0; i < n; ++i) {
    ic[i] = i;
    ir[i] = i;
  }

  /*   Elimination. */

  *det = 1.;

  for (jj = 0; jj < n - 1; ++jj) {
    ipiv = jj;
    jpiv = jj;
    piv = 0.;
    for (i = jj; i < n; ++i) {
      for (j = jj; j < n; ++j) {
	p = fabs(ARRAY2D(a, ir[i], ic[j]));
	if (p > piv) {
	  piv = p;
	  ipiv = i;
	  jpiv = j;
	}
      }
    }
    *det *= ARRAY2D(a, ir[ipiv], ic[jpiv]);
#define GE_PIVOTS_DEBUG
#ifdef GE_PIVOTS_DEBUG
    if(jj==0)  fprintf(fp9,"\n Pivots in GE");
    if((jj%6)==0) fprintf(fp9,"\n");
    fprintf(fp9," %4d %12.3e ",jj,fabs(ARRAY2D(a, ir[ipiv], ic[jpiv])));
#endif
    if (ipiv != jj) {
      *det = -(*det);
    }
    if (jpiv != jj) {
      *det = -(*det);
    }

    if (piv < RSMALL) {
      fprintf(fp9,"         NOTE:Pivot %3li < %10.3E, in GE\n",jj,RSMALL);	
    }

    k = ir[jj];
    ir[jj] = ir[ipiv];
    ir[ipiv] = k;

    k = ic[jj];
    ic[jj] = ic[jpiv];
    ic[jpiv] = k;

    jjp1 = jj + 1;
    for (l = jjp1; l < n; ++l) {
      
      rm = ARRAY2D(a, ir[l], ic[jj]) / ARRAY2D(a, ir[jj], ic[jj]);
      if (rm != 0.) {
	for (i = jjp1; i < n; ++i) {
	  ARRAY2D(a, ir[l], ic[i]) -= rm * ARRAY2D(a, ir[jj], ic[i]);
	}
	if (nrhs != 0) {
	  for (irh = 0; irh < nrhs; ++irh) {
	    ARRAY2D(f, ir[l], irh) -= rm * ARRAY2D(f, ir[jj], irh);
	  }
	}
      }
    }
  }
  *det *= ARRAY2D(a, ir[n - 1], ic[n - 1]);
#ifdef GE_PIVOTS_DEBUG
     if((jj%6)==0) fprintf(fp9,"\n");
     fprintf(fp9," %4d %12.3e \n",n-1,ARRAY2D(a, ir[n - 1], ic[n - 1]));
#endif

  if (nrhs == 0) {
    /* Since the pointers for ir and ic were adjusted we need to put the
       back here. */
    /*ir \+= 1;*/
    /*ic \+= 1;*/
    free(ir);
    free(ic);
    return 0;
  }

  /*   Backsubstitution : */

  for (irh = 0; irh < nrhs; ++irh) {
#ifndef FLOATING_POINT_TRAP
    if(ARRAY2D(a, ir[n - 1], ic[n - 1]) == 0) {
      printf("Division by Zero, exiting\n");
      exit(0);
    }
#endif
    ARRAY2D(u, ic[n - 1], irh) = ARRAY2D(f, ir[n - 1], irh) / ARRAY2D(a, ir[n - 1], ic[n - 1]);
    for (i1 = 0; i1 < n - 1; ++i1) {
      i = n - (i1 + 1) - 1;
      sm = 0.;
      ip1 = i + 1;
      for (j = ip1; j < n; ++j) {
	sm += ARRAY2D(a, ir[i], ic[j]) * ARRAY2D(u, ic[j], irh);
      }
      ARRAY2D(u, ic[i], irh) = (ARRAY2D(f, ir[i], irh) - sm) / ARRAY2D(a, ir[i], ic[i]);
    }
  }

  /* Since the pointers for ir and ic were adjusted we need to put the
       back here. */
  /*ir \+= 1;*/
  /*ic \+= 1;*/
  free(ir);
  free(ic);

  return 0;
} /* ge_ */


/*     ---------- ------ */
/* Subroutine */ int 
newlab(iap_type *iap, rap_type *rap)
{

    

  /* Local variables */
  integer mlab, ibrs, nars;

  integer labrs, nskip, nfprs, itprs, iswrs, ntplrs, ntotrs, lab, 
    ibr, mbr, ips, itp, irs, isw;
  logical eof3;

  /* Determine a suitable label when restarting. */

  
  ips = iap->ips;
  irs = iap->irs;
  isw = iap->isw;
  itp = iap->itp;

  mbr = 0;
  mlab = 0;
  rewind(fp3);

 L1:
  if (fscanf(fp3,"%ld",&ibrs) != 1) {
    goto L2;
  }
  if (fscanf(fp3,"%ld",&ntotrs) != 1) {
    goto L2;
  }
  if (fscanf(fp3,"%ld",&itprs) != 1) {
    goto L2;
  }
  if (fscanf(fp3,"%ld",&labrs) != 1) {
    goto L2;
  }
  if (fscanf(fp3,"%ld",&nfprs) != 1) {
    goto L2;
  }
  if (fscanf(fp3,"%ld",&iswrs) != 1) {
    goto L2;
  }
  if (fscanf(fp3,"%ld",&ntplrs) != 1) {
    goto L2;
  }
  if (fscanf(fp3,"%ld",&nars) != 1) {
    goto L2;
  }
  if (fscanf(fp3,"%ld",&nskip) != 1) {
    goto L2;
  }
  /*go to the end of the line*/
  while(fgetc(fp3)!='\n');

  if (ibrs > mbr) {
    mbr = ibrs;
  }
  if (labrs > mlab) {
    mlab = labrs;
  }
  skip3(&nskip, &eof3);
  if (! eof3) {
    goto L1;
  }

 L2:
  lab = mlab;
  iap->lab = lab;
  if (isw < 0 || irs == 0) {
    ibr = mbr + 1;
    iap->ibr = ibr;
  } else if ((abs(itp) < 10 && abs(isw) == 2) || (ips == 2 && itp == 3) || 
	     (ips == 4 && isw == 2 && abs(itp) < 10) || (ips == 5 && itp % 10 == 2)) {
    ibr = irs;
    iap->ibr = ibr;
  }

  return 0;
} /* newlab_ */


int 
findlb(iap_type *iap, const rap_type *rap, 
       integer irs, integer *nfpr, logical *found)
{
  /* Local variables */
  integer nars;

  integer labrs, nskip, itpst, iswrs, ntplrs, ntotrs, ibr, itp, 
    isw;
  logical eof3;

  
  /* Locates restart point with label IRS and determines type. */
  /* If the label can not be located on unit 3 then FOUND will be .FALSE. */
  

  
  *found = FALSE_;
  rewind(fp3);
  isw = iap->isw;
  
  while(1) {
    if (fscanf(fp3,"%ld",&ibr) != 1) {
      break;
    }
    if (fscanf(fp3,"%ld",&ntotrs) != 1) {
      break;
    }
    if (fscanf(fp3,"%ld",&itp) != 1) {
      break;
    }
    if (fscanf(fp3,"%ld",&labrs) != 1) {
      break;
    }
    if (fscanf(fp3,"%ld",&(*nfpr)) != 1) {
      break;
    }
    if (fscanf(fp3,"%ld",&iswrs) != 1) {
      break;
    }
    if (fscanf(fp3,"%ld",&ntplrs) != 1) {
      break;
    }
    if (fscanf(fp3,"%ld",&nars) != 1) {
      break;
    }
    if (fscanf(fp3,"%ld",&nskip) != 1) {
      break;
    }
    /*go to the end of the line*/
    while(fgetc(fp3)!='\n');
    iap->itp = itp;
    iap->ibr = ibr;
    if (labrs == irs) {
      *found = TRUE_;
      if (abs(isw) == 2) {
	if (abs(itp) < 10) {
	  itpst = abs(itp);
	  iap->itpst = itpst;
	} else {
	  itpst = abs(itp / 10);
	  iap->itpst = itpst;
	}
      } else {
        itpst = 0;
	iap->itpst = itpst;
      }
      fseek(fp3,-2,SEEK_CUR);
      while((fgetc(fp3)!='\n') && (ftell(fp3)!=1)){
	fseek(fp3,-2,SEEK_CUR);
      }
      return 0;
    } else {
      skip3(&nskip, &eof3);
      if (eof3) {
	break;
      }
    }
  }
  return 0;
} 


/*     ---------- ------ */
/* Subroutine */ int 
readlb(const iap_type *iap, const rap_type *rap, doublereal *u, doublereal *par)
{
  /* Local variables */
  integer labr, ndim, ibrr, itpr, iswr, i;
  doublereal t;
  integer nparr, nfprr, n1, n2, ntotr, nskipr, ntplrs, nar;

  /* Reads the restart data for algebraic problems. */

  fscanf(fp3,"%ld",&ibrr);
  fscanf(fp3,"%ld",&ntotr);
  fscanf(fp3,"%ld",&itpr);
  fscanf(fp3,"%ld",&labr);
  fscanf(fp3,"%ld",&nfprr);
  fscanf(fp3,"%ld",&iswr);
  fscanf(fp3,"%ld",&ntplrs);
  fscanf(fp3,"%ld",&nar);
  fscanf(fp3,"%ld",&nskipr);
  fscanf(fp3,"%ld",&n1);
  fscanf(fp3,"%ld",&n2);
  fscanf(fp3,"%ld",&nparr);
  ndim = nar - 1;
  fscanf(fp3,"%le",&t);
  for (i = 0; i < ndim; ++i) {
    fscanf(fp3,"%le",&u[i]);
  }
  if (nparr > NPARX) {
    nparr = NPARX;
    printf("Warning : NPARX too small for restart data :\n restart PAR(i) skipped for i > %3ld\n",nparr);
  }
  for (i = 0; i < nparr; ++i) {
    fscanf(fp3,"%le",&par[i]);
  }

  return 0;
} /* readlb_ */


/*     ---------- ----- */
/* Subroutine */ int 
skip3(integer *nskip, logical *eof3)
{

  /* Local variables */
  integer i;


  /* Skips the specified number of lines on unit 3. */


  *eof3 = FALSE_;
  for (i = 0; i < *nskip; ++i) {
    /* NOTE from Randy:  I am not 100% happy with this.  I am 
       not sure if this properly simulates the Fortran behavior */
    while(1) {
      int tmp = fgetc(fp3);
      if(tmp==EOF) {
	*eof3 = TRUE_;
	return 0;
      }
      if((char)tmp=='\n') {
	break;
      }
    }
  }
  return 0;
} /* skip3_ */


/*     ------ --------- -------- ----- */
doublereal 
rinpr(iap_type *iap, integer *ndim1, integer *ndxloc, doublereal *ups, doublereal *vps, doublereal *dtm, doublereal *thu)
{
  /* System generated locals */
  integer ups_dim1, vps_dim1;
  doublereal ret_val;

  /* Local variables */
  integer ndim, ncol;

  integer ntst, i, j, k;
  doublereal s;
  integer k1;
  doublereal sj, *wi;
  integer jp1;

  wi = (doublereal *)malloc(sizeof(doublereal)*(iap->ncol+1));

  /* Computes the L2 inner product of UPS and VPS. */
  /* (Using the first NDIM1 components only.) */

  /* Local */

  /* Parameter adjustments */
  /*--dtm;*/
  /*--thu;*/
  vps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
    
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;

  /* Weights for the integration formulae : */
  wint(ncol + 1, wi);

  s = 0.;
  for (j = 0; j < ntst; ++j) {
    jp1 = j + 1;
    sj = 0.;
    for (i = 0; i < *ndim1; ++i) {
      for (k = 0; k < ncol; ++k) {
	k1 = k * ndim + i;
	sj += wi[k] * thu[i] * ARRAY2D(ups, j, k1) * ARRAY2D(vps, j, k1);
      }
      sj += wi[ncol] * thu[i] * ARRAY2D(ups, jp1, i) * ARRAY2D(vps, jp1, i);
    }
    s += dtm[j] * sj;
  }

  ret_val = s;
  free(wi);

  return ret_val;
} /* rinpr_ */


/*     ------ --------- -------- ------ */
doublereal 
rnrmsq(iap_type *iap, integer *ndim1, integer *ndxloc, doublereal *ups, doublereal *dtm, doublereal *thu)
{
  /* System generated locals */
  doublereal ret_val;

    /* Local variables */






/* Finds the norm of UPS (first NDIM1 components are included only). */

  /* Parameter adjustments */
  /*--thu;*/

  
  ret_val = rinpr(iap, ndim1, ndxloc, ups, ups, dtm, thu);

  return ret_val;
} /* rnrmsq_ */


/*     ------ --------- -------- ----- */
doublereal 
rintg(iap_type *iap, integer *ndxloc, integer ic, doublereal *ups, doublereal *dtm)
{
  /* System generated locals */
  integer ups_dim1;
  doublereal ret_val;

  /* Local variables */
  integer ndim, ncol;

  integer ntst, j, k;
  doublereal s;
  integer k1;
  doublereal sj, *wi;
  integer jp1;

  wi = (doublereal *)malloc(sizeof(doublereal)*(iap->ncol+1));

  /* Computes the integral of the IC'th component of UPS. */

/* Local */

  /* Parameter adjustments */
  /*--dtm;*/
  ups_dim1 = *ndxloc;
  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;

  /* Weights for the integration formulae : */
  wint(ncol + 1, wi);
  s = 0.;
  for (j = 0; j < ntst; ++j) {
    jp1 = j + 1;
    sj = 0.;
    for (k = 0; k < ncol; ++k) {
      k1 = k * ndim + ic - 1;
      sj += wi[k] * ARRAY2D(ups, j, k1);
    }
    sj += wi[ncol] * ARRAY2D(ups, jp1, (ic - 1));
    s += dtm[j] * sj;
  }

  ret_val = s;

  free(wi);
  return ret_val;
} /* rintg_ */


/*     ------ --------- -------- ----- */
doublereal 
rnrm2(iap_type *iap, integer *ndxloc, integer *ic, doublereal *ups, doublereal *dtm)
{
  /* System generated locals */
  integer ups_dim1;
  doublereal ret_val;

  


  /* Local variables */
  integer ndim, ncol;

  integer ntst, j, k;
  doublereal s;
  integer k1;
  doublereal sj, *wi;
  integer jp1;

  wi = (doublereal *)malloc(sizeof(doublereal)*(iap->ncol+1));

  /* Computes the L2-norm of the IC'th component of UPS. */

/* Local */

  /* Parameter adjustments */
  /*--dtm;*/
  ups_dim1 = *ndxloc;
  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;

  /* Weights for the integration formulae : */
  wint(ncol + 1, wi);
  s = 0.;
  for (j = 0; j < ntst; ++j) {
    jp1 = j + 1;
    sj = 0.;
    for (k = 0; k < ncol; ++k) {
      k1 = k * ndim + *ic - 1;
      /* Computing 2nd power */
      sj += wi[k] * (ARRAY2D(ups, j, k1) * ARRAY2D(ups, j, k1));
    }
    /* Computing 2nd power */
    sj += wi[ncol] * (ARRAY2D(ups, jp1, (*ic - 1)) * ARRAY2D(ups, jp1, (*ic - 1)));
    s += dtm[j] * sj;
  }

  ret_val = sqrt(s);
  free(wi);

  return ret_val;
} /* rnrm2_ */


/*     ------ --------- -------- ------ */
doublereal 
rmxups(iap_type *iap, integer *ndxloc, integer *i, doublereal *ups)
{
  /* System generated locals */
  integer ups_dim1;
  doublereal ret_val;

  /* Local variables */
  integer ndim, ncol, ntst, j, k, k1;


  

  /* Computes the maximum of the I'th component of UPS. */


  /* Parameter adjustments */
  ups_dim1 = *ndxloc;
  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;

  ret_val = ARRAY2D(ups, 0, (*i - 1));

  for (j = 0; j < ntst; ++j) {
    for (k = 0; k < ncol; ++k) {
      k1 = k * ndim + *i - 1;
      if (ARRAY2D(ups, j, k1) > ret_val) {
	ret_val = ARRAY2D(ups, j, k1);
      }
    }
  }
  if (ARRAY2D(ups, ntst, (*i - 1)) > ret_val) {
    ret_val = ARRAY2D(ups, ntst, (*i - 1));
  }

  return ret_val;
} /* rmxups_ */


/*     ------ --------- -------- ------ */
doublereal 
rmnups(iap_type *iap, integer *ndxloc, integer *i, doublereal *ups)
{
  /* System generated locals */
  integer ups_dim1;
  doublereal ret_val;

  /* Local variables */
  integer ndim, ncol, ntst, j, k, k1;


  

  /* Computes the minimum of the I'th component of UPS. */


  /* Parameter adjustments */
  ups_dim1 = *ndxloc;
  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;

  ret_val = ARRAY2D(ups, 0, (*i - 1));

  for (j = 0; j < ntst; ++j) {
    for (k = 0; k < ncol; ++k) {
      k1 = k * ndim + (*i - 1);
      if (ARRAY2D(ups, j, k1) < ret_val) {
	ret_val = ARRAY2D(ups, j, k1);
      }
    }
  }
  if (ARRAY2D(ups, ntst, (*i - 1)) < ret_val) {
    ret_val = ARRAY2D(ups, ntst, (*i - 1 ));
  }

  return ret_val;
} /* rmnups_ */


/*     ---------- ------ */
/* Subroutine */ int 
scaleb(iap_type *iap, integer *icp, integer *ndxloc, doublereal *dvps, doublereal *rld, doublereal *dtm, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer dvps_dim1;

  /* Local variables */
  integer ndim, ncol, nfpr, nrow, ntst, i, j;
  doublereal sc, ss;

/* Scales the vector (DVPS,RLD) so its norm becomes 1. */


  /* Parameter adjustments */
  /*--icp;*/
  /*--rld;*/
  /*--dtm;*/
  /*--thl;*/
  /*--thu;*/
  dvps_dim1 = *ndxloc;

  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  nfpr = iap->nfpr;

  ss = rnrmsq(iap, &ndim, ndxloc, dvps, dtm, thu);
  for (i = 0; i < nfpr; ++i) {
    /* Computing 2nd power */
    ss += thl[icp[i]] * (rld[i] * rld[i]);
  }

  sc = 1. / sqrt(ss);

  nrow = ndim * ncol;
  for (j = 0; j < ntst; ++j) {
    for (i = 0; i < nrow; ++i) {
      ARRAY2D(dvps, j, i) *= sc;
    }
  }

  for (i = 0; i < ndim; ++i) {
    ARRAY2D(dvps, ntst, i) *= sc;
  }

  for (i = 0; i < nfpr; ++i) {
    rld[i] = sc * rld[i];
  }

  return 0;
} /* scaleb_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*                    General Boundary Value Problems */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ------ */
/* Subroutine */ int 
cnrlbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), STPNT_TYPE_BVP((*stpnt)), PVLI_TYPE_BVP((*pvli)), doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{
    /* Local variables */
  integer iads, ndim, ncol;
  integer ntot, ntst, iuzr, nuzr, i, j, k;
  doublereal dsold;
  integer nodir = 0;
  doublereal rlold[NPARX], rldot[NPARX], rlcur[NPARX];
  integer nitps, istop;
  integer itpst;
  doublereal ds;
  doublereal bp1;
  doublereal sp1;
  integer iad;
  integer ilp, ips, isp, irs;
  doublereal rds, rlp;
  integer isw, itp;
  doublereal *uzr;
    
  integer ntst_plus_one=iap->ntst + 1;

  doublereal *ups;
  doublereal *upoldp;
  doublereal *uoldps;
  doublereal *udotps;
  doublereal *dups;
  doublereal *fa;
  doublereal *dtm;
  doublereal *tm;
  doublereal *fc;
  doublereal *p0;
  doublereal *p1;
  doublecomplex *ev;
  ups=(doublereal *)malloc(sizeof(doublereal)*(iap->ntst + 1)*(iap->ndim * iap->ncol)); 
  upoldp=(doublereal *)malloc(sizeof(doublereal)*(iap->ntst + 1)*(iap->ndim * iap->ncol)); 
  uoldps=(doublereal *)malloc(sizeof(doublereal)*(iap->ntst + 1)*(iap->ndim * iap->ncol)); 
  udotps=(doublereal *)malloc(sizeof(doublereal)*(iap->ntst + 1)*(iap->ndim * iap->ncol)); 
  dups=(doublereal *)malloc(sizeof(doublereal)*(iap->ntst + 1)*(iap->ndim * iap->ncol)); 
  fa=(doublereal *)malloc(sizeof(doublereal)*(iap->ntst + 1)*(iap->ndim * iap->ncol)); 
  dtm=(doublereal *)malloc(sizeof(doublereal)*(iap->ntst + 1)); 
  tm=(doublereal *)malloc(sizeof(doublereal)*(iap->ntst + 1)); 
  fc=(doublereal *)malloc(sizeof(doublereal)*(iap->nbc+iap->nint+1)); 
  p0=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*(iap->ndim));
  p1=(doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*(iap->ndim));
  ev=(doublecomplex *)malloc(sizeof(doublecomplex)*(iap->ndim));
  uzr=(doublereal *)malloc(sizeof(doublereal)*(iap->nuzr));
    

  /* INITIALIZE COMPUTATION OF BRANCH */

  /* Parameter adjustments */
  /*--iuz;*/
  /*--icp;*/
  /*--par;*/
    
    

  
  ndim = iap->ndim;
  ips = iap->ips;
  irs = iap->irs;
  ilp = iap->ilp;
  ntst = iap->ntst;
  ncol = iap->ncol;
  iad = iap->iad;
  iads = iap->iads;
  isp = iap->isp;
  isw = iap->isw;
  nuzr = iap->nuzr;
  itpst = iap->itpst;

  ds = rap->ds;

  rds = ds;
  dsold = rds;
  rap->dsold = dsold;
  if (isp < 0) {
    isp = -isp;
    iap->isp = isp;
  }
  sp1 = 0.;
  bp1 = 0.;
  rlp = 0.;
  if (nuzr > 0) {
    for (i = 0; i < nuzr; ++i) {
      uzr[i] = 0.;
    }
  }
  nitps = 0;
  iap->nit = nitps;
  ntot = 0;
  iap->ntot = ntot;
  istop = 0;
  iap->istop = istop;

  for (j = 0; j < (iap->ntst + 1); ++j) {
    for (i = 0; i < (iap->ndim * iap -> ncol); ++i) {
      ups[j + i * (iap->ntst + 1)] = 0.;
      uoldps[j + i * (iap->ntst + 1)] = 0.;
      upoldp[j + i * (iap->ntst + 1)] = 0.;
      dups[j + i * (iap->ntst + 1)] = 0.;
      udotps[j + i * (iap->ntst + 1)] = 0.;
      fa[j + i * (iap->ntst + 1)] = 0.;
    }
  }
  
  for(j = 0; j < NPARX; ++j) {
    rldot[j] = 0.0;
    rlcur[j] = 0.0;
    rlold[j] = 0.0;
  }

  rsptbv(iap, rap, par, icp, funi, stpnt, &rds, rlcur, 
	 rlold, rldot, &ntst_plus_one, ups, uoldps, udotps, upoldp, dups, tm, 
	 dtm, ev, &nodir, thl, thu);
  (*pvli)(iap, rap, icp, dtm, &ntst_plus_one, ups, &ndim, p0, p1, par);

  setrtn(iap, &ntst, &ntst_plus_one, ups, par);

  if (nodir == 1 && isw > 0) {
    stdrbv(iap, rap, par, icp, funi, bcni, icni, rlcur, 
	   rlold, rldot, iap->ntst+1, ups, dups, uoldps, udotps, upoldp, 
	   fa, fc, dtm, 0, p0, p1, thl, thu);
  } else if (irs != 0 && isw < 0) {
    stdrbv(iap, rap, par, icp, funi, bcni, icni, rlcur, 
	   rlold, rldot, iap->ntst+1, ups, dups, uoldps, udotps, upoldp, 
	   fa, fc, dtm, 1, p0, p1, thl, thu);
  }

  /* Store plotting data for restart point : */

  sthd(iap, rap, par, icp, thl, thu);
  if (irs == 0) {
    itp = itpst * 10 + 9;
  } else {
    itp = 0;
  }
  iap->itp = itp;
  istop = 0;
  iap->istop = istop;
  (*pvli)(iap, rap, icp, dtm, &ntst_plus_one, ups, &ndim, p0, p1, par);
  stplbv(iap, rap, par, icp, rldot, &ntst_plus_one, ups, udotps,
	 tm, dtm, thl, thu);
  istop = iap->istop;
  if (istop == 1) {
    free(ups);
    free(upoldp);
    free(uoldps);
    free(udotps);
    free(dups);
    free(fa);
    free(dtm);
    free(tm);
    free(fc);
    free(p0);
    free(p1);
    free(ev);
    free(uzr);
    return 0;
  }

  extrbv(iap, rap, funi, &rds, rlcur, rlold, rldot, &ntst_plus_one, ups,
	 uoldps, udotps);

  itp = 0;
  iap->itp = itp;
  goto L2;

 L1:
  itp = 0;
  iap->itp = itp;
  ntot = iap->ntot;

/* Adapt the mesh to the solution. */

  if (iad != 0) {
    if (ntot % iad == 0) {
      adapt(iap, rap, &ntst, &ncol, &ntst, &ncol, tm, dtm, &
	    ntst_plus_one, ups, uoldps);
    }
  }

  /* Adapt the stepsize along the branch. */

  if (iads != 0) {
    if (ntot % iads == 0) {
      adptds(iap, rap, &rds);
    }
  }

  /* Provide initial approximation and determine next point. */

#define SECANT_GUESS
#ifdef SECANT_GUESS
  contbv(iap, rap, par, icp, funi, &rds, rlcur, rlold, 
	 rldot, &ntst_plus_one, ups, uoldps, udotps, upoldp, dtm, thl, thu);
#else
  {
    doublereal *uolddotps,*rlolddot;
    uolddotps = (doublereal *)malloc(sizeof(doublereal)*(iap->ntst+1)*(iap->ndim)*(iap->ncol));
    rlolddot = (doublereal *)malloc(sizeof(doublereal)*(iap->nfpr));

    for (i = 0; i < iap->nfpr; ++i) {
      rlolddot[i] = rldot[i];
    }
    for (j = 0; j < iap->ntst; ++j) {
      for (i = 0; i < iap->ncol * iap->ndim; ++i) {
	uolddotps[j + i * (iap->ntst+1)] = udotps[j + i * (iap->ntst+1)];
      }
    }
    for (i = 0; i < iap->ndim; ++i) {
      uolddotps[ntst + i * (iap->ntst+1)] = udotps[ntst + i * (iap->ntst+1)]; 
    }
    
    stdrbv(iap, rap, par, icp, funi, bcni, icni, rlcur, 
	   rlold, rldot, iap->ntst+1, ups, dups, uoldps, udotps, upoldp, 
	   fa, fc, dtm, 0, p0, p1, thl, thu);
	
    {
      doublereal dot_product = rinpr(iap,&(iap->ndim),&ntst_plus_one,udotps,uolddotps,dtm,thu);
      
      for (i = 0; i < iap->nfpr; ++i) {
	/* Computing 2nd power */
	//FIXME  No sure if this is right
	//This is the original
	//dot_product += thl[icp[-1+i]] * (rldot[i] * rlolddot[i]);
	dot_product += thl[icp[i]] * (rldot[i] * rlolddot[i]);
      }
      
      if(dot_product < 0) {
	for (i = 0; i < iap->nfpr; ++i) {
	  rldot[i] = -rldot[i];
	}
	for (j = 0; j < iap->ntst; ++j) {
	  for (i = 0; i < iap->ndim * iap->ncol; ++i) {
	    udotps[j + i * (iap->ntst+1)] = -udotps[j + i * (iap->ntst+1)];
	  }
	}
	for (i = 0; i < iap->ndim; ++i) {
	  udotps[ntst + i * (iap->ntst+1)] = -udotps[ntst + i * (iap->ntst+1)]; 
	}
      }
    }
    free(uolddotps);
    free(rlolddot);
  }	  
  extrbv(iap, rap, funi, &rds, rlcur, rlold, rldot, &ntst_plus_one, ups,
	 uoldps, udotps);

  stupbv(iap, rap, par, icp, funi, rlcur, rlold, rldot,
	 &ntst_plus_one, ups, uoldps, upoldp);
#endif
  
 L2:
  stepbv(iap, rap, par, icp, funi, bcni, icni, pvli, &rds, 
	 rlcur, rlold, rldot, &ntst_plus_one, ups, dups, uoldps, udotps, upoldp,
	 fa, fc, tm, dtm, p0, p1, thl, thu);
  istop = iap->istop;
  if (istop == 1) {
    goto L3;
  }

  /* Check for user supplied parameter output parameter-values. */

  if (nuzr > 0) {
    for (iuzr = 0; iuzr < nuzr; ++iuzr) {
      iap->iuzr = iuzr;
      lcspbv(iap, rap, par, icp, fnuzbv, funi, bcni, 
	     icni, pvli, &uzr[iuzr], rlcur, rlold, rldot, &
	     ntst_plus_one, ups, dups, uoldps, udotps, upoldp, fa, fc, tm, 
	     dtm, p0, p1, ev, thl, thu, iuz, vuz);
      istop = iap->istop;
      if (istop == 1) {
	goto L3;
      }
      itp = iap->itp;
      if (itp == -1) {
	if (iuz[iuzr] >= 0) {
	  itp = -4 - itpst * 10;
	  iap->itp = itp;
	  for (k = 0; k < nuzr; ++k) {
	    uzr[k] = 0.;
	  }
	} else {
	  istop = -1;
	  iap->istop = istop;
	}
      }
    }
  }

  /* Check for fold. */

  if (ilp == 1) {
    lcspbv(iap, rap, par, icp, fnlpbv, funi, bcni, icni,
	   pvli, &rlp, rlcur, rlold, rldot, &ntst_plus_one, ups, dups, 
	   uoldps, udotps, upoldp, fa, fc, tm, dtm, p0, p1, ev, thl, thu,
	   iuz, vuz);
    istop = iap->istop;
    if (istop == 1) {
      goto L3;
    }
    itp = iap->itp;
    if (itp == -1) {
      itp = itpst * 10 + 5;
      iap->itp = itp;
      rlp = 0.;
      bp1 = 0.;
      sp1 = 0.;
    }
  }

  /* Check for branch point. */

  if (isp >= 2) {
    lcspbv(iap, rap, par, icp, fnbpbv, funi, bcni, icni,
	   pvli, &bp1, rlcur, rlold, rldot, &ntst_plus_one, ups, dups, 
	   uoldps, udotps, upoldp, fa, fc, tm, dtm, p0, p1, ev, thl, thu,
	   iuz, vuz);
    istop = iap->istop;
    if (istop == 1) {
      goto L3;
    }
    itp = iap->itp;
    if (itp == -1) {
      itp = itpst * 10 + 6;
      iap->itp = itp;
      rlp = 0.;
      bp1 = 0.;
      sp1 = 0.;
    }
  }

  /* Check for period-doubling and torus bifurcation. */

  if ((isp == 1 || isp == 2) && (ips == 2 || ips == 7 || ips == 12)) {
    lcspbv(iap, rap, par, icp, fnspbv, funi, bcni, icni,
	   pvli, &sp1, rlcur, rlold, rldot, &ntst_plus_one, ups, dups, 
	   uoldps, udotps, upoldp, fa, fc, tm, dtm, p0, p1, ev, thl, thu,
	   iuz, vuz);
    istop = iap->istop;
    if (istop == 1) {
      goto L3;
    }
    itp = iap->itp;
    if (itp == -1) {
      /*          **Secondary periodic bifurcation: determine type */
      tpspbv(iap, rap, par, icp, ev);
      rlp = 0.;
      bp1 = 0.;
      sp1 = 0.;
    }
  }

  /* Store plotting data. */

 L3:
  (*pvli)(iap, rap, icp, dtm, &ntst_plus_one, ups, &ndim, p0, p1, par);
  stplbv(iap, rap, par, icp, rldot, &ntst_plus_one, ups, udotps,
	 tm, dtm, thl, thu);

  istop = iap->istop;
  if (istop == 0) {
    goto L1;
  } else {
    free(ups);
    free(upoldp);
    free(uoldps);
    free(udotps);
    free(dups);
    free(fa);
    free(dtm);
    free(tm);
    free(fc);
    free(p0);
    free(p1);
    free(ev);
    free(uzr);
    return 0;
  }
    

} /* cnrlbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
contbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), doublereal *rds, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *dtm, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer ups_dim1, udotps_dim1, uoldps_dim1;

  /* Local variables */
  integer ndim, ncol, nfpr, nrow, ntst, i, j;
  doublereal dsold;

  doublereal dds;


  

  /* Determines an initial approximation to the next solution point, */
  /* by a computation of the null space of the Jacobian. */
  /* The stepsize used in the preceding step has been stored in DSOLD. */


  /* Parameter adjustments */
  /*--dtm;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  udotps_dim1 = *ndxloc;
  uoldps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  nfpr = iap->nfpr;

  dsold = rap->dsold;

  /* Compute rate of change (along branch) of PAR(ICP(1)) and U : */

  dds = 1. / dsold;
  nrow = ndim * ncol;
  for (j = 0; j < ntst + 1; ++j) {
    for (i = 0; i < nrow; ++i) {
      ARRAY2D(udotps, j, i) = (ARRAY2D(ups, j, i) - ARRAY2D(uoldps, j, i)) * dds;
    }
  }
  for (i = 0; i < nfpr; ++i) {
    rldot[i] = (rlcur[i] - rlold[i]) * dds;
  }
  /*        Rescale, to set the norm of (UDOTPS,RLDOT) equal to 1. */
  scaleb(iap, icp, ndxloc, udotps, rldot, dtm, 
	 thl, thu);

  /* Extrapolate to get initial approximation to next solution point. */

  extrbv(iap, rap, funi, rds, rlcur, rlold, rldot, 
	 ndxloc, ups, uoldps, 
	 udotps);

  /* Store time-derivative. */

  stupbv(iap, rap, par, icp, funi, rlcur, rlold, rldot,
	 ndxloc, ups, uoldps, upoldp);

  return 0;
} /* contbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
extrbv(iap_type *iap, rap_type *rap, FUNI_TYPE((*funi)), doublereal *rds, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *uoldps, doublereal *udotps)
{
  /* System generated locals */
  integer ups_dim1, udotps_dim1, uoldps_dim1;

  /* Local variables */
  integer ndim, ncol, nfpr, nrow, ntst, i, j;


  

  /* Determines an initial approximation to the next solution by */
  /* a computation of the null space of the Jacobian. */
  /* The stepsize used in the preceding step has been stored in DSOLD. */



  /* Parameter adjustments */
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  udotps_dim1 = *ndxloc;
  uoldps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  nfpr = iap->nfpr;

  nrow = ndim * ncol;
  for (i = 0; i < nfpr; ++i) {
    rlold[i] = rlcur[i];
    rlcur[i] += *rds * rldot[i];
  }
  for (j = 0; j < ntst + 1; ++j) {
    for (i = 0; i < nrow; ++i) {
      ARRAY2D(uoldps, j, i) = ARRAY2D(ups, j, i);
      ARRAY2D(ups, j, i) += *rds * ARRAY2D(udotps, j, i);
    }
  }

  return 0;
} /* extrbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
stupbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *uoldps, doublereal *upoldp)
{
  /* System generated locals */
  integer ups_dim1, uoldps_dim1, upoldp_dim1;

  /* Local variables */
  integer ndim, ncol;
  integer nfpr, ntst;
  integer i, j, k;
  integer n1, ips;
  doublereal *dfdp,*dfdu,*uold,*f,*u;

  dfdp = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*NPARX);
  dfdu = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*(iap->ndim));
  uold = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  f    = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));
  u    = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));


/* Stores U-prime (derivative with respect to T) in UPOLDP. */


/* Local */

  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  upoldp_dim1 = *ndxloc;
  uoldps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
  
  ndim = iap->ndim;
  ips = iap->ips;
  ntst = iap->ntst;
  ncol = iap->ncol;
  nfpr = iap->nfpr;

  for (i = 0; i < nfpr; ++i) {
    par[icp[i]] = rlold[i];
  }

  for (j = 0; j < ntst + 1; ++j) {
    for (i = 0; i < ndim; ++i) {
      u[i] = ARRAY2D(uoldps, j, i);
      if (ips == 14 || ips == 16) {
	uold[i] = ARRAY2D(uoldps, j, i) * 2 - ARRAY2D(ups, j, i);
      } else {
	uold[i] = ARRAY2D(uoldps, j, i);
      }
    }
    (*funi)(iap, rap, ndim, u, uold, icp, par, 0, f, 
	    dfdu, dfdp);
    for (i = 0; i < ndim; ++i) {
      ARRAY2D(upoldp, j, i) = f[i];
    }
  }

  for (k = 1; k <= ncol - 1; ++k) {
    n1 = k * ndim;
    for (j = 0; j < ntst; ++j) {
      for (i = 0; i < ndim; ++i) {
	u[i] = ARRAY2D(uoldps, j, (n1 + i));
	if (ips == 14 || ips == 16) {
	  uold[i] = ARRAY2D(uoldps, j, (n1 + i)) * 2 
	    - ARRAY2D(ups, j, (n1 + i));
	} else {
	  uold[i] = ARRAY2D(uoldps, j, (n1 + i));
	}
      }
      (*funi)(iap, rap, ndim, u, uold, icp, par, 0, f, dfdu, dfdp);
      for (i = 0; i < ndim; ++i) {
	ARRAY2D(upoldp, j, (n1 + i)) = f[i];
      }
    }
  }

  for (i = 0; i < nfpr; ++i) {
    par[icp[i]] = rlcur[i];
  }
  free(dfdp);
  free(dfdu);
  free(uold);
  free(f   );
  free(u   );


  return 0;
} /* stupbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
stepbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), PVLI_TYPE_BVP((*pvli)), doublereal *rds, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *dups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *fa, doublereal *fc, doublereal *tm, doublereal *dtm, doublereal *p0, doublereal *p1, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer ups_dim1, uoldps_dim1, udotps_dim1, fa_dim1;

  /* Local variables */
  integer iads;
  doublereal adrl;
  logical done;
  integer ndim, ncol;
  doublereal epsl, rdrl;
  integer nfpr, ifst;
  doublereal epsu;
  integer ntop, itnw, nllv;
  doublereal dumx;
  integer ntot, nrow=0, nwtn, ntst, i, j;
  doublereal dsold, dsmin;
  integer nitps;
  doublereal rdumx;
  integer istop;
  doublereal au;

  doublereal delref=0.0, delmax;

  integer iid;
  doublereal adu;
  integer ibr, mxt;
  doublereal umx;
  integer nit1;

  /* Controls the solution of the nonlinear equations (by Newton's method) 
*/
/* for the next solution (PAR(ICP(*)) , U) on a branch of solutions. */



  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  /*--fc;*/
  /*--tm;*/
  /*--dtm;*/
  fa_dim1 = *ndxloc;
  udotps_dim1 = *ndxloc;
  uoldps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
    
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  iads = iap->iads;
  iid = iap->iid;
  itnw = iap->itnw;
  nwtn = iap->nwtn;
  nfpr = iap->nfpr;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  dsmin = rap->dsmin;
  epsl = rap->epsl;
  epsu = rap->epsu;

 L1:
  dsold = *rds;
  rap->dsold = dsold;
  nitps = 0;
  iap->nit = nitps;

  /* Write additional output on unit 9 if requested. */

  wrtbv9(iap, rap, par, icp, rlcur, ndxloc, 
	 ups, tm, dtm, thl, thu);

  /* Generate the Jacobian matrix and the right hand side. */

  for (nit1 = 1; nit1 <= itnw; ++nit1) {

    nitps = nit1;
    iap->nit = nitps;
    nllv = 0;

    ifst = 0;
    if (nitps <= nwtn) {
      ifst = 1;
    }

    solvbv(&ifst, iap, rap, par, icp, funi, bcni, icni, 
	   rds, &nllv, rlcur, rlold, rldot, ndxloc, 
	   ups, dups, uoldps, 
	   udotps, upoldp, dtm, fa, fc, p0, 
	   p1, thl, thu);
    /* Add Newton increments. */

    for (i = 0; i < ndim; ++i) {
      ARRAY2D(ups, ntst, i) += fc[i];
    }
    for (i = 0; i < nfpr; ++i) {
      rlcur[i] += fc[ndim + i];
      par[icp[i]] = rlcur[i];
    }

    dumx = 0.;
    umx = 0.;
    nrow = ndim * ncol;
    for (j = 0; j < ntst; ++j) {
      for (i = 0; i < nrow; ++i) {
	adu = fabs(ARRAY2D(fa, j, i));
	if (adu > dumx) {
	  dumx = adu;
	}
	au = fabs(ARRAY2D(ups, j, i));
	if (au > umx) {
	  umx = au;
	}
	ARRAY2D(ups, j, i) += ARRAY2D(fa, j, i);
      }
    }
    /* I am not sure why this is here, and it requires a wierd
       special case for the parsing of points which are "special"
       (i.e. bifurcations, folds, etc).  I may want to get rid
       of it at some time.
    */
    wrtbv9(iap, rap, par, icp, rlcur, ndxloc, ups, tm, dtm, thl, thu);

    /* Check whether user-supplied error tolerances have been met : */

    done = TRUE_;
    rdrl = 0.;
    for (i = 0; i < nfpr; ++i) {
      adrl = fabs(fc[ndim + i]) / (fabs(rlcur[i]) + 1.);
      if (adrl > epsl) {
	done = FALSE_;
      }
      if (adrl > rdrl) {
	rdrl = adrl;
      }
    }
    rdumx = dumx / (umx + 1.);
    if (done && rdumx < epsu) {
      (*pvli)(iap, rap, icp, dtm, ndxloc, ups, &ndim, p0, p1, par);
      if (iid >= 2) {
	fprintf(fp9," \n");	
      }
      return 0;
    }

    if (nitps == 1) {
      delref = max(rdrl,rdumx) * 20;
    } else {
      delmax = max(rdrl,rdumx);
      if (delmax > delref) {
	goto L3;
      }
    }

    /* L2: */
  }

  /* Maximum number of iterations reached. */

 L3:
  if (iads == 0 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li NOTE:No convergence with fixed step size\n",ibr,ntop);	
  }
  if (iads == 0) {
    goto L13;
  }

  /* Reduce stepsize and try again. */

  mxt = itnw;
  iap->nit = mxt;
  adptds(iap, rap, rds);
  if (fabs(*rds) < dsmin) {
    goto L12;
  }
  for (i = 0; i < nfpr; ++i) {
    rlcur[i] = rlold[i] + *rds * rldot[i];
  }
  for (j = 0; j < ntst + 1; ++j) {
    for (i = 0; i < nrow; ++i) {
      ARRAY2D(ups, j, i) = ARRAY2D(uoldps, j, i) + *rds * 
	ARRAY2D(udotps, j, i);
    }
  }
  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li NOTE:Retrying step\n",ibr,ntop);	

  }
  goto L1;

/* Minimum stepsize reached. */

 L12:
  if (iap->mynode == 0) {
    fprintf(fp9,"%4li%6li, NOTE:No convergence using minimum step size\n",ibr,ntop);	

  }
 L13:
  for (i = 0; i < nfpr; ++i) {
    rlcur[i] = rlold[i];
    par[icp[i]] = rlcur[i];
  }
  for (j = 0; j < ntst + 1; ++j) {
    for (i = 0; i < nrow; ++i) {
      ARRAY2D(ups, j, i) = ARRAY2D(uoldps, j, i);
    }
  }
  istop = 1;
  iap->istop = istop;


  return 0;
} /* stepbv_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*      Restart of Solution Branches ( Differential Equations ) */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ------ */
/* Subroutine */ int 
rsptbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), STPNT_TYPE_BVP((*stpnt)), doublereal *rds, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *dups, doublereal *tm, doublereal *dtm, doublecomplex *ev, integer *nodir, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer ups_dim1, uoldps_dim1;

    /* Local variables */
  integer ndim, ncol, nfpr, ntst, i, j;

  integer ntsrs;

  integer ncolrs;

  integer irs, isw;

  integer ntst_fort8,ncol_fort8,junk;
  

  /* Restarts computation of a branch of solutions at point labelled IRS. */
  /* The output written on unit 8 by a previous run is now expected as */
  /* input on unit 3. The label IRS, where computation is to resume, must */
  /* be specified in the user-supplied subroutine INIT. */
  /* If IRS=0 then the starting point must be provided analytically in the 
*/
/* user-supplied subroutine STPNT. */




  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  /*--tm;*/
  /*--dtm;*/
  /*--ev;*/
  uoldps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
    
  ndim = iap->ndim;
  irs = iap->irs;
  ntst = iap->ntst;
  ncol = iap->ncol;
  isw = iap->isw;
  nfpr = iap->nfpr;

  /* Get restart data : */

    /* First take a peak at the file to see if ntst and
       ncol are different then the values found in
       the parameter file fort.2.
    */
  if(iap->irs > 0) {
    findlb(iap, rap, iap->irs, &junk, &junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&junk);
    fscanf(fp3,"%ld",&ntst_fort8);
    fscanf(fp3,"%ld",&ncol_fort8);
  } else {
    ntst_fort8 = iap->ntst;
    ncol_fort8 = iap->ncol;
  }

  {
    integer ntst_used,ncol_used;
    doublereal *ups_new;
    doublereal *upoldp_new;
    doublereal *udotps_new;
    doublereal *tm_new;
    doublereal *dtm_new;
    integer ndxloc_orig=*ndxloc;

    /* use the bigger of the size defined in fort.2 and the one defined in fort.8 */
    if(ntst_fort8 > ntst)
      ntst_used=ntst_fort8;
    else
      ntst_used=ntst;

    if(ncol_fort8 > ncol)
      ncol_used=ncol_fort8;
    else
      ncol_used=ncol;
      
    *ndxloc=(ntst_used + 1)*4;
    ups_new    = (doublereal *)malloc(sizeof(doublereal)*(*ndxloc)*(iap->ndim * ncol_used));
    upoldp_new    = (doublereal *)malloc(sizeof(doublereal)*(*ndxloc)*(iap->ndim * ncol_used));
    udotps_new = (doublereal *)malloc(sizeof(doublereal)*(*ndxloc)*(iap->ndim * ncol_used));
    tm_new     = (doublereal *)malloc(sizeof(doublereal)*(*ndxloc));
    dtm_new    = (doublereal *)malloc(sizeof(doublereal)*(*ndxloc));

    /*initialize arrays*/
    for(i=0;i<*ndxloc;i++) {
      dtm_new[i]=0.0;
      tm_new[i]=0.0;
      for(j=0;j<ndim*ncol_used;j++) {
	ups_new[i+j*(*ndxloc)] = 0.0;
	upoldp_new[i+j*(*ndxloc)] = 0.0;
	udotps_new[i+j*(*ndxloc)] = 0.0;
      }
    }
    (*stpnt)(iap, rap, par, icp, &ntsrs, &ncolrs, rlcur, 
	     rldot, ndxloc, ups_new, udotps_new, upoldp_new, tm_new, 
	     dtm_new, nodir, thl, thu);
    /* Determine a suitable starting label and branch number. */
      
    newlab(iap, rap);
      
    for (j = 0; j < ntsrs; ++j) {
      dtm_new[j] = tm_new[j + 1] - tm_new[j];
    }
      
    /* Adapt mesh if necessary : */
      
    if (ntst != ntsrs || ncol != ncolrs) {
      adapt(iap, rap, &ntsrs, &ncolrs, &ntst, &ncol, tm_new, dtm_new, 
	    ndxloc, ups_new, udotps_new);
    }
    /* Copy from the temporary large arrays into the normal arrays. */
    for(i=0;i<ntst+1;i++) {
      dtm[-1 + i+1] = dtm_new[i];
      tm[-1 + i+1] = tm_new[i];
      for(j=0;j<ndim*ncol;j++) {
	ups[i+j*(ntst+1)] = ups_new[i+j*(*ndxloc)];
	upoldp[i+j*(ntst+1)] = upoldp_new[i+j*(*ndxloc)];
	udotps[i+j*(ntst+1)] = udotps_new[i+j*(*ndxloc)];
      }
    }
    *ndxloc=ndxloc_orig;
    free(ups_new);
    free(upoldp_new);
    free(udotps_new);
    free(tm_new);
    free(dtm_new);
  }

  /* Set UOLDPS, RLOLD. */

  for (i = 0; i < nfpr; ++i) {
    rlcur[i] = par[icp[i]];
    rlold[i] = rlcur[i];
  }

  for (i = 0; i < ndim * ncol; ++i) {
    for (j = 0; j < ntst + 1; ++j) {
      ARRAY2D(uoldps, j, i) = ARRAY2D(ups, j, i);
    }
  }

  /* Store U-prime (derivative with respect to time or space variable). */

  if (*nodir == -1) {
    /*        ** Restart from a Hopf bifurcation. */
    *nodir = 0;
    isw = 1;
  } else {
    /*        ** Restart from orbit. */
    stupbv(iap, rap, par, icp, funi, rlcur, rlold, rldot, ndxloc, ups, uoldps, 
	   upoldp);
  }

  return 0;
} /* rsptbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
stpnbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, integer *ntsrs, integer *ncolrs, doublereal *rlcur, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *udotps, doublereal *upoldp, doublereal *tm, doublereal *dtm, integer *nodir, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer ups_dim1, udotps_dim1;

    /* Local variables */
  integer ndim, nars;
  doublereal temp[7];
  integer nfpr;

  integer i, j, k;

  logical found;
  integer icprs[NPARX], nparr, nskip;

  integer nfprs, k1, k2, itprs, iswrs, nskip1, nskip2;

  integer ndimrd, ndimrs, ntplrs, ntotrs, lab, ibr, ips, irs, isw;
  logical eof3;


  /* This subroutine locates and retrieves the information required to */
  /* restart computation at the point with label IRS. */
  /* This information is expected on unit 3. */

/* Local */


  /* Parameter adjustments */
  /*--tm;*/
  /*--dtm;*/
  /*--par;*/
  /*--icp;*/
  /*--rlcur;*/
  /*--rldot;*/
  udotps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
    
  ndim = iap->ndim;
  ips = iap->ips;
  irs = iap->irs;
  isw = iap->isw;
  nfpr = iap->nfpr;

  findlb(iap, rap, irs, &nfprs, &found);
  fscanf(fp3,"%ld",&ibr);
  fscanf(fp3,"%ld",&ntotrs);
  fscanf(fp3,"%ld",&itprs);
  fscanf(fp3,"%ld",&lab);
  fscanf(fp3,"%ld",&nfprs);
  fscanf(fp3,"%ld",&iswrs);
  fscanf(fp3,"%ld",&ntplrs);
  fscanf(fp3,"%ld",&nars);
  fscanf(fp3,"%ld",&nskip);
  fscanf(fp3,"%ld",&(*ntsrs));
  fscanf(fp3,"%ld",&(*ncolrs));
  fscanf(fp3,"%ld",&nparr);
  iap->ibr = ibr;
  iap->lab = lab;

  ndimrs = nars - 1;
  nskip1 = (ndimrs + 1) / 8 - ndim / 7;
  nskip2 = (ndimrs + 1) / 9 - ndim / 8;
  if (ndim <= ndimrs) {
    ndimrd = ndim;
  } else {
    ndimrd = ndimrs;
  }

  for (j = 0; j < *ntsrs; ++j) {
    for (i = 0; i < *ncolrs; ++i) {
      k1 = i * ndim;
      k2 = k1 + ndimrd - 1;
      fscanf(fp3,"%le",&temp[i]);
      for (k = k1; k <= k2; ++k) {
	fscanf(fp3,"%lf",&ARRAY2D(ups, j, k));
      }
      /*go to the end of the line*/
      while(fgetc(fp3)!='\n');
	    
      if (nskip1 > 0) {
	skip3(&nskip1, &eof3);
      }
    }
    tm[j] = temp[0];
  }
  fscanf(fp3,"%le",&tm[*ntsrs]);
  for (k = 0; k < ndimrd; ++k) {
    fscanf(fp3,"%le",&ARRAY2D(ups, *ntsrs, k));
  }
  /*go to the end of the line*/
  while(fgetc(fp3)!='\n');
  if (nskip1 > 0) {
    skip3(&nskip1, &eof3);
  }

  for (i = 0; i < nfprs; ++i) {
    fscanf(fp3,"%ld",&icprs[i]);
  }
  for (i = 0; i < nfprs; ++i) {
    fscanf(fp3,"%le",&rldot[i]);
  }

  /* Read U-dot (deriv. with respect to arclength along solution branch). */

  for (j = 0; j < *ntsrs; ++j) {
    for (i = 0; i < *ncolrs; ++i) {
      k1 = i * ndim;
      k2 = k1 + ndimrd - 1;
      for (k = k1; k <= k2; ++k) {
	fscanf(fp3,"%le",&ARRAY2D(udotps, j, k));
      }
      /*go to the end of the line*/
      while(fgetc(fp3)!='\n');
      if (nskip2 > 0) {
	skip3(&nskip2, &eof3);
      }
    }
  }
  for (k = 0; k < ndimrd; ++k) {
    fscanf(fp3,"%le",&ARRAY2D(udotps, *ntsrs, k));
  }
  /*go to the end of the line*/
  while(fgetc(fp3)!='\n');
  if (nskip2 > 0) {
    skip3(&nskip2, &eof3);
  }

  /* Read the parameter values. */

  if (nparr > NPARX) {
    nparr = NPARX;
    printf("Warning : NPARX too small for restart data :\n restart PAR(i) skipped for i > %3ld\n",nparr);
  }
  for (i = 0; i < nparr; ++i) {
    fscanf(fp3,"%le",&par[i]);
  }
  for (i = 0; i < nfpr; ++i) {
    rlcur[i] = par[icp[i]];
  }

  /* Special case : Preprocess restart data in case of homoclinic */
  /* continuation */

  if (ips == 9) {
    preho(ndxloc, ntsrs, &ndimrd, &ndim, ncolrs, ups, udotps, tm, par);

    /* Special case : Preprocess restart data in case of branch switching 
     */
    /* at a period doubling bifurcation. */

  } else if ((ips == 2 || ips == 6) && isw == -1 && itprs == 7) {
    pdble(iap, rap, &ndim, ntsrs, ncolrs, ndxloc, ups, udotps, tm, par);
    return 0;
  }

  /* Take care of the case where the free parameters have been changed at */
  /* the restart point. */

  if (nfprs != nfpr) {
    *nodir = 1;
    return 0;
  }
  for (i = 0; i < nfpr; ++i) {
    if (icprs[i] != icp[i]) {
      *nodir = 1;
      return 0;
    }
  }

  return 0;
} /* stpnbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
stpnub(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, integer *ntsrs, integer *ncolrs, doublereal *rlcur, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *udotps, doublereal *upoldp, doublereal *tm, doublereal *dtm, integer *nodir, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer ups_dim1, udotps_dim1;

  /* Local variables */
  integer ndim, ncol, nfpr, ntst, ncol1, i, j, k;
  doublereal t, *u;
  integer k1, k2;

  doublereal dt;
  integer lab, ibr;

  u = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim));

/* Generates a starting point for the continuation of a branch of */
/* of solutions to general boundary value problems by calling the user */
/* supplied subroutine STPNT where an analytical solution is given. */

/* Local */

  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--rlcur;*/
  /*--rldot;*/
  /*--tm;*/
  /*--dtm;*/
  udotps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;
    
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  nfpr = iap->nfpr;

/* Generate the (initially uniform) mesh. */

  msh(iap, rap, tm);
  dt = 1. / (ntst * ncol);

  for (j = 0; j < ntst + 1; ++j) {
    if (j == ntst) {
      ncol1 = 1;
    } else {
      ncol1 = ncol;
    }
    for (i = 0; i < ncol1; ++i) {
      t = tm[j] + i * dt;
      k1 = i * ndim;
      k2 = (i + 1) * ndim;
      stpnt(ndim, t, u, par);
      for (k = k1; k < k2; ++k) {
	ARRAY2D(ups, j, k) = u[k - k1];
      }
    }
  }

  *ntsrs = ntst;
  *ncolrs = ncol;
  ibr = 1;
  iap->ibr = ibr;
  lab = 0;
  iap->lab = lab;

  for (i = 0; i < nfpr; ++i) {
    rlcur[i] = par[icp[i]];
  }

  *nodir = 1;

  free(u);
  return 0;
} /* stpnub_ */


/*     ---------- ------ */
/* Subroutine */ int 
setrtn(iap_type *iap, integer *ntst, integer *ndxloc, doublereal *ups, doublereal *par)
{
  /* System generated locals */
  integer ups_dim1;

  /* Local variables */
  integer i;

  integer nbc;

  /* Initialization for rotations */
  
  
  /* Parameter adjustments */
  /*--par;*/
  ups_dim1 = *ndxloc;
  
  par[18] = pi(2.0);
  nbc = iap->nbc;

  global_rotations.irtn = 0;
  for (i = 0; i < nbc; ++i) {
    doublereal tmp = (ARRAY2D(ups, *ntst, i) - ARRAY2D(ups, 0, i)) / pi(2.0);
    global_rotations.nrtn[i] = i_dnnt(&tmp);

    if (global_rotations.nrtn[i] != 0) {
      global_rotations.irtn = 1;
    }
  }

  return 0;
} /* setrtn_ */

/*     ---------- ------ */
/* Subroutine */ int 
stdrbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer ndxloc, doublereal *ups, doublereal *dups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *fa, doublereal *fc, doublereal *dtm, integer iperp, doublereal *p0, doublereal *p1, doublereal *thl, doublereal *thu)
{
  /* System generated locals */
  integer udotps_dim1, fa_dim1;

  /* Builtin Local */

    /* variables functions */
  integer ndim, ncol, nfpr, ifst, nllv;
  doublereal rdsz;
  integer nrow, ntst, i, j;

  integer iid;




/* Generates a direction vector (UDOTPS,RLDOT) that is needed to start */
/* the computation of a branch when no direction vector is given. */



/* Generate the Jacobian matrix with zero direction vector. */
/* (Then the last row of the Jacobian will be zero) */
/* in case the starting direction is to be determined. */

  /* Parameter adjustments */
  /*--par;*/
  /*--icp;*/
  /*--rlcur;*/
  /*--rlold;*/
  /*--rldot;*/
  /*--fc;*/
  /*--dtm;*/
  fa_dim1 = ndxloc;
  udotps_dim1 = ndxloc;
  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  iid = iap->iid;
  nfpr = iap->nfpr;

  nrow = ndim * ncol;
  if (iperp == 0) {
    for (j = 0; j < ntst + 1; ++j) {
      for (i = 0; i < nrow; ++i) {
	ARRAY2D(udotps, j, i) = 0.;
      }
    }
    for (i = 0; i < nfpr; ++i) {
      rldot[i] = 0.;
    }
  }

  rdsz = 0.;
  nllv = 1;
  ifst = 1;
  solvbv(&ifst, iap, rap, par, icp, funi, bcni, icni, &rdsz, 
	 &nllv, rlcur, rlold, rldot, &ndxloc, ups, dups, uoldps, 
	 udotps, upoldp, dtm, fa, fc, 
	 p0, p1, thl, thu);

  /* Compute the starting direction. */

  for (i = 0; i < ndim; ++i) {
    ARRAY2D(udotps, ntst, i) = fc[i];
  }
  for (i = 0; i < nfpr; ++i) {
    rldot[i] = fc[ndim + i];
    par[icp[i]] = rlcur[i];
  }

  for (j = 0; j < ntst; ++j) {
    for (i = 0; i < nrow; ++i) {
      ARRAY2D(udotps, j, i) = ARRAY2D(fa, j, i);
    }
  }

  /* Scale the starting direction. */

  scaleb(iap, icp, &ndxloc, udotps, rldot, dtm, thl, thu);

  /* Make sure that RLDOT(1) is positive (unless zero). */

  if (rldot[0] < 0.) {
    for (i = 0; i < nfpr; ++i) {
      rldot[i] = -rldot[i];
    }
    for (j = 0; j < ntst; ++j) {
      for (i = 0; i < nrow; ++i) {
	ARRAY2D(udotps, j, i) = -ARRAY2D(udotps, j, i);
      }
    }
#define BUG_FIX
#ifdef BUG_FIX
    for (i = 0; i < ndim; ++i) {
      ARRAY2D(udotps, ntst, i) = -ARRAY2D(udotps, ntst, i); 
    }
#endif
  }

  if (iap->mynode > 0) {
    return 0;
  }

  if (iid >= 2) {
    fprintf(fp9,"Starting direction of the free parameter(s) :\n");	

    for (i = 0; i < nfpr; ++i) {
      fprintf(fp9," PAR(%3ld) :%11.3E\n",icp[i],rldot[i]);	
    }
  }


  return 0;
} /* stdrbv_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*  Detection and Location of Branch Points in Boundary Value Problems */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ------ */
/* Subroutine */ int 
lcspbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, FNCS_TYPE_BVP((*fncs)), FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), PVLI_TYPE_BVP((*pvli)), doublereal *q, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *dups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *fa, doublereal *fc, doublereal *tm, doublereal *dtm, doublereal *p0, doublereal *p1, doublecomplex *ev, doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{
  /* Local variables */
  logical chng;
  doublereal epss;
  integer ntop, itmx;
  doublereal rtmp, rrds;
  integer ntot;
  doublereal s, dsold, dsmax;
  integer istop;
  doublereal q0, q1, s0, s1;
  integer nitsp1;
  doublereal dq, ds, pq;

  integer iid, ibr;
  doublereal rds;
  integer itp;

  /* This subroutine uses the Secant method to accurately locate folds */
  /* branch points, and zero(es) of user parameter values. */
  /* Such points are located as points on a solution branch where the */
  /* EXTERNAL function FNCS changes sign. */
  /* It involves calling the basic solution subroutines CONTBV and STEP */
  /* with decreasing values of RDS (stepsize along branch). */
  /* The point is assumed to have been found with sufficient accuracy if */
  /* the ratio between RDS and the user supplied value of DS is less than */
  /* the user-supplied tolerance EPSS. */
  /* This subroutine is called from CNRLB, which controls the computation */
  /* of branches of solutions to general boundary value problems. */

  /* Parameter adjustments */
  /*--ev;*/
  /*--dtm;*/
  /*--tm;*/
  /*--rldot;*/
  /*--rlold;*/
  /*--rlcur;*/

    
  iid = iap->iid;
  itmx = iap->itmx;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  ds = rap->ds;
  dsmax = rap->dsmax;
  dsold = rap->dsold;
  epss = rap->epss;

  /* Check for zero. */

  q0 = *q;
  q1 = (*fncs)(iap, rap, par, icp, &chng, funi, bcni, icni, p0, p1, 
	       ev, rlcur, rlold, rldot, ndxloc, ups, uoldps, 
	       udotps, upoldp, fa, fc, dups, tm, dtm, thl, thu, iuz, vuz)
    ;

  pq = q0 * q1;
  if (pq >= 0. || ! chng) {
    *q = q1;
    return 0;
  }

  /* Use the secant method for the first step: */

  s0 = 0.;
  s1 = dsold;
  nitsp1 = 0;
  dq = q0 - q1;
  rds = q1 / dq * (s1 - s0);
  rtmp = HMACH1;
 L1:
  rds = rtmp * rds;
  s = s1 + rds;

  /* Return if tolerance has been met : */

  rrds = fabs(rds) / (sqrt(fabs(ds * dsmax)) + 1);
  if (rrds < epss) {
    itp = -1;
    iap->itp = itp;
    /* xx???   Q=0.d0 */
    fprintf(fp9,"==> Location of special point : Convergence.    Stepsize =%11.3E\n",rds);	

    return 0;
  }

  /* If requested write additional output on unit 9 : */

  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9," ==> Location of special point :  Iteration %3ld   Stepsize =%11.3E\n",nitsp1,rds);	
  }

  contbv(iap, rap, par, icp, funi, &rds, rlcur, rlold, &
	 rldot[0], ndxloc, ups, uoldps, udotps, upoldp, dtm, thl, thu);
  stepbv(iap, rap, par, icp, funi, bcni, icni, pvli, &rds, &rlcur[-1 + 
								 1], rlold, rldot, ndxloc, ups, dups, uoldps, udotps, 
	 upoldp, fa, fc, tm, dtm, p0, p1, thl, thu);
  istop = iap->istop;
  if (istop != 0) {
    *q = 0.;
    return 0;
  }

  /* Check for zero. */

  *q = (*fncs)(iap, rap, par, icp, &chng, funi, bcni, icni, p0, p1, 
	       ev, rlcur, rlold, rldot, ndxloc, ups, uoldps, 
	       udotps, upoldp, fa, fc, dups, tm, dtm, thl, thu, iuz, vuz)
    ;

  ++nitsp1;
  if (nitsp1 <= itmx) {
    /*        Use Mueller's method with bracketing for subsequent steps */
    mueller(&q0, &q1, q, &s0, &s1, &s, &rds);
    goto L1;
  }

  if (iap->mynode > 0) {
    return 0;
  }

  fprintf(fp9,"%4li%6li NOTE:Possible special point\n",ibr,ntop);	
  *q = 0.;

  return 0;
} /* lcspbv_ */


/*     ------ --------- -------- ------ */
doublereal 
fnlpbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, logical *chng, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), doublereal *p0, doublereal *p1, doublecomplex *ev, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *fa, doublereal *fc, doublereal *dups, doublereal *tm, doublereal *dtm, doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{
  /* System generated locals */
  integer udotps_dim1, fa_dim1;
  doublereal ret_val;

    

    /* Local variables */
  integer ndim, ncol, nfpr, ifst, nllv, ntop;
  doublereal rdsz;
  integer ntot, ntst, i, j;

  integer iid, ibr, isp;

/* RETURNS A QUANTITY THAT CHANGES SIGN AT A LIMIT POINT (BVP) */


  fa_dim1 = *ndxloc;
  udotps_dim1 = *ndxloc;

  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  isp = iap->isp;
  iid = iap->iid;
  nfpr = iap->nfpr;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  /* Find the direction vector. */

  nllv = -1;
  ifst = 0;
  rdsz = 0.;

  solvbv(&ifst, iap, rap, par, icp, funi, bcni, icni, &rdsz, &nllv,
	 rlcur, rlold, rldot, ndxloc, ups, dups, uoldps, udotps, 
	 upoldp, dtm, fa, fc, 
	 p0, p1, thl, thu);

  for (i = 0; i < ndim; ++i) {
    ARRAY2D(udotps, ntst, i) = fc[i];
  }

  for (i = 0; i < nfpr; ++i) {
    rldot[i] = fc[ndim + i];
  }

  for (j = 0; j < ntst; ++j) {
    for (i = 0; i < ndim*ncol; ++i) {
      ARRAY2D(udotps, j, i) = ARRAY2D(fa, j, i);
    }
  }

  /* Scale the direction vector. */

  scaleb(iap, icp, ndxloc, udotps, rldot, dtm, 
	 thl, thu);
  if (iid >= 2 && iap->mynode == 0) {
    fprintf(fp9,"%4li%6li        Fold Function %14.6E\n",abs(ibr),ntop,rldot[0]);	
  }

  /* Set the quantity to be returned. */

  ret_val = rldot[0];
  *chng = TRUE_;
  rap->fldf = ret_val;


  return ret_val;
} /* fnlpbv_ */


/*     ------ --------- -------- ------ */
doublereal 
fnbpbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, logical *chng, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), doublereal *p0, doublereal *p1, doublecomplex *ev, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *fa, doublereal *fc, doublereal *dups, doublereal *tm, doublereal *dtm, doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{
  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  integer ndim, ntop, ntot;
  doublereal f;
  integer i;
  doublereal u;

  doublereal *pp;
  integer iid;
  doublereal det;
  integer ibr;
  doublereal det0;

  pp = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*(iap->ndim));


    
  ndim = iap->ndim;
  iid = iap->iid;

/* Save the determinant of the reduced system. */

  det = rap->det;
  det0 = det;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  /* Compute the determinant of P1. */

/* Computing 2nd power */
  for (i = 0; i < ndim * ndim; ++i) {
    pp[i] = p1[i];
  }
  ge(ndim, ndim, pp, 0, 1, &u, 1, &f, &det);
  rap->det = det;

  /* Set the determinant of the normalized reduced system. */

  if (det != 0.) {
    ret_val = det0 / det;
    *chng = TRUE_;
  } else {
    ret_val = 0.;
    *chng = FALSE_;
  }
  rap->biff = ret_val;

  if (iap->mynode > 0) {
    free(pp);
    return ret_val;
  }

  if (iid >= 2) {
    fprintf(fp9,"%4li%6li        BP   Function %14.6E\n",abs(ibr),ntop,ret_val);	

  }
  free(pp);
  return ret_val;
} /* fnbpbv_ */


/*     ------ --------- -------- ------ */
doublereal 
fnspbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, logical *chng, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), doublereal *p0, doublereal *p1, doublecomplex *ev, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *fa, doublereal *fc, doublereal *dups, doublereal *tm, doublereal *dtm, doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{
  /* System generated locals */
  doublereal ret_val;

  /* Local variables */
  doublereal amin;
  integer ndim, nins, ntop, ntot;
  doublecomplex ztmp;
  integer nins1=0;
  doublereal d;
  integer i, j;

  integer iid,ibr,loc=0,isp,isw;
  doublereal *wrk, azm1;

  wrk = (doublereal *)malloc(sizeof(doublereal)*(iap->ndim)*(iap->ndim));

  /* This function returns a quantity that changes sign when a complex */
  /* pair of eigenvalues of the linearized Poincare map moves in or out */
  /* of the unit circle or when a real eigenvalues passes through -1. */
  /* Local */



  /* Parameter adjustments */
  /*--ev;*/
  /*--p1;*/
  /*--p0;*/
    
    

  
  ndim = iap->ndim;
  isp = iap->isp;
  isw = iap->isw;
  iid = iap->iid;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  /* Initialize. */

  ret_val = 0.;
  rap->spbf = ret_val;
  d = 0.;
  *chng = FALSE_;

/*  Compute the Floquet multipliers */

  flowkm(&ndim, p0, p1, &iid, wrk, ev);
  free(wrk);
  /* Find the multiplier closest to z=1. */
  /* send_mult here! */
  send_mult(ibr,ntot+1,ndim,(doublecomplex *)&ev[0]);
  amin = RLARGE;
  for (j = 0; j < ndim; ++j) {
    doublecomplex tmp;
    tmp.r = ev[j].r - 1., tmp.i = ev[j].i;
    azm1 = z_abs(&tmp);
    if (azm1 <= amin) {
      amin = azm1;
      loc = j;
    }
  }
  if (loc != 0) {
    ztmp.r = ev[loc].r, ztmp.i = ev[loc].i;
    ev[loc].r = ev[0].r, ev[loc].i = ev[0].i;
    ev[0].r = ztmp.r, ev[0].i = ztmp.i;
  }

  /* Order the remaining Floquet multipliers by distance from |z|=1. */

  if (ndim >= 3) {
    for (i = 1; i < ndim - 1; ++i) {
      amin = RLARGE;
      for (j = i; j < ndim; ++j) {
	azm1 = z_abs(&ev[j]) - 1.;
	azm1 = fabs(azm1);
	if (azm1 <= amin) {
	  amin = azm1;
	  loc = j;
	}
      }
      if (loc != i) {
	ztmp.r = ev[loc].r, ztmp.i = ev[loc].i;
	ev[loc].r = ev[i].r, ev[loc].i = ev[i].i;
	ev[i].r = ztmp.r, ev[i].i = ztmp.i;
      }
    }
  }

  /* Print error message if the Floquet multiplier at z=1 is inaccurate. */
  /* (ISP is set to negative and detection of bifurations is discontinued) 
   */

  {
    doublecomplex tmp;
    tmp.r = ev[0].r - 1., tmp.i = ev[0].i;
    amin = z_abs(&tmp);
  }
  if (amin > (double).05 && isp == 2) {
    if (iap->mynode == 0) {
      if (iid >= 2) {
	fprintf(fp9,"%4li%6li NOTE:Multiplier inaccurate\n",abs(ibr),ntop);	

      }
      for (i = 0; i < ndim; ++i) {
	fprintf(fp9,"%4li%6li        Multiplier %3li %14.6E %14.6E\n",abs(ibr),ntop,i,ev[i].r,ev[i].i);	
      }
    }
    nins = 0;
    iap->nins = nins;
    if (iap->mynode == 0) {
      fprintf(fp9,"%4li%6li        Multipliers:   Stable: %3li\n",abs(ibr),ntop,nins);	
    }
    isp = -isp;
    iap->isp = isp;
    return ret_val;
  }

  /* Restart automatic detection if the Floquet multiplier at z=1 is */
  /* sufficiently accurate again. */

  if (isp < 0) {
    if (amin < (double).01) {
      if (iap->mynode == 0) {
	fprintf(fp9,"%4li%6li NOTE:Multiplier accurate again\n",abs(ibr),ntop);	
      }
      isp = -isp;
      iap->isp = isp;
    } else {
      if (iap->mynode == 0) {
	for (i = 0; i < ndim; ++i) {
	  fprintf(fp9,"%4li%6li        Multiplier %3li %14.6E %14.6E\n",abs(ibr),ntop,i,ev[i].r,ev[i].i);	
	}
      }
      return ret_val;
    }
  }

  /* Count the number of Floquet multipliers inside the unit circle. */

  if (ndim == 1) {
    d = 0.;
    ret_val = d;
    rap->spbf = ret_val;
  } else {
    nins1 = 1;
    for (i = 1; i < ndim; ++i) {
      if (z_abs(&ev[i]) <= 1.) {
	++nins1;
      }
    }
    if (isp == 2) {
      if (d_imag(&ev[1]) == 0. && ev[1].r > 0.) {
	/*            *Ignore if second multiplier is real positive */
	d = 0.;
      } else {
	d = z_abs(&ev[1]) - 1.;
      }
      if (isw == 2) {
	ret_val = 0.;
      } else {
	ret_val = d;
      }
      rap->spbf = ret_val;
      nins = iap->nins;
      if (nins1 != nins) {
	*chng = TRUE_;
      }
    }
  }

  nins = nins1;
  iap->nins = nins;
  if (iid >= 2 && (isp == 1 || isp == 2)) {
    if (iap->mynode == 0) {
      fprintf(fp9,"%4li%6li        SPB  Function %14.6E\n",abs(ibr),ntop,d);	

    }
  }

  /* Print the Floquet multipliers. */

  nins = iap->nins;
  if (iap->mynode == 0) {
    fprintf(fp9,"%4li%6li        Multipliers:   Stable: %3li\n",abs(ibr),ntop,nins);	

    for (i = 0; i < ndim; ++i) {
      fprintf(fp9,"%4li%6li        Multiplier %3li %14.6E %14.6E\n",abs(ibr),ntop,i,ev[i].r,ev[i].i);	

    }
  }


  return ret_val;
} /* fnspbv_ */


/*     ------ --------- -------- ------ */
doublereal 
fnuzbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, logical *chng, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), doublereal *p0, doublereal *p1, doublecomplex *ev, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *uoldps, doublereal *udotps, doublereal *upoldp, doublereal *fa, doublereal *fc, doublereal *dups, doublereal *tm, doublereal *dtm, doublereal *thl, doublereal *thu, integer *iuz, doublereal *vuz)
{
  /* System generated locals */
  doublereal ret_val;

    

    /* Local variables */
  integer ntop, ntot, iuzr, iid, ibr;

  
  iid = iap->iid;
  iuzr = iap->iuzr;
  ibr = iap->ibr;
  ntot = iap->ntot;
  ntop = (ntot + 1) % 10000;

  ret_val = par[abs(iuz[iuzr])] - vuz[iuzr];
  *chng = TRUE_;

  if (iid >= 3) {
    fprintf(fp9,"%4li%6li        User Func. %3li %14.6E\n",abs(ibr),ntop,iuzr,ret_val);	
  }

  return ret_val;
} /* fnuzbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
tpspbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublecomplex *ev)
{

    /* Local variables */
  doublereal amin;
  integer ndim;
  doublereal epss, d;
  integer i, itpst;
  doublereal ad;
  integer loc, itp, loc1;
  doublereal azm1;


  /* Determines type of secondary periodic bifurcation. */

  
  ndim = iap->ndim;

  epss = rap->epss;
  itpst = iap->itpst;

  /* Find the eigenvalue closest to z=1. */

  loc = 1;
  amin = RLARGE;
  for (i = 0; i < ndim; ++i) {
    doublecomplex tmp;
    tmp.r = ev[i].r - 1., tmp.i = ev[i].i;
    azm1 = z_abs(&tmp);
    if (azm1 <= amin) {
      amin = azm1;
      loc = i;
    }
  }

  /* Find the eigenvalue closest to the unit circle */
  /* (excluding the eigenvalue at z=1). */

  loc1 = 1;
  amin = RLARGE;
  for (i = 0; i < ndim; ++i) {
    if (i != loc) {
      d = z_abs(&ev[i]) - 1.;
      ad = fabs(d);
      if (ad <= amin) {
	amin = ad;
	loc1 = i;
      }
    }
  }

  if (fabs(d_imag(&ev[loc1])) > sqrt(epss)) {
    /*       ** torus bifurcation */
    itp = itpst * 10 + 8;
    iap->itp = itp;
    par[11] =fabs(atan2(d_imag(&ev[loc1]),ev[loc1].r));

  } else /* if(complicated condition) */ {
    if (ev[loc1].r < -.5) {
      /*       ** period doubling */
      itp = itpst * 10 + 7;
      iap->itp = itp;
    } else {
      /*       ** something else... */
      itp = 0;
      iap->itp = itp;
    }
  }

  return 0;
} /* tpspbv_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*                    Output (Boundary Value Problems) */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*     ---------- ------ */
/* Subroutine */ int 
stplbv(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *udotps, doublereal *tm, doublereal *dtm, doublereal *thl, doublereal *thu)
{
    /* Local variables */
  integer labw, ndim, ibrs, nins, iplt, itmp, jtmp, ntot;
  int iflag=0;
  integer i;
  double u_high[1000],u_low[1000],u_0[1000],u_bar[1000];
  doublereal a0, a1;
   /* used a fixed array here for maximum AUTO size 
     double u_high[NAUTO], u_low, u0, ubar   */

  integer istop, itpst, n2, ntots;

  doublereal rl0, rl1;

  integer iab, lab, ibr, ndm;
  doublereal amp;
  integer ips, itp, irs, npr, isw, nmx;
  doublereal umx[7];

  /* Writes the bifurcation diagram on unit 7 (Differential Equations) */
  /* (Also controls the writing of complete solutions on unit 8). */
  /* Every line written contains, in order, the following: */

/*  IBR    : The label of the branch. */
/*  NTOT   : The index of the point on the branch. */
/*           (Points are numbered consecutively along a branch). */
/*           If IPS=2 or 3, then the sign of NTOT indicates stability : */
/*            - = stable , + = unstable, or unknown. */
/*  ITP    : An integer indicating the type of point : */

/*             4  (  )  :   Output point (Every NPR steps along branch). 
*/
/*            -4  (UZ)  :   Output point (Zero of user function). */
/*             5  (LP)  :   Fold (fold). */
/*             6  (BP)  :   Branch point. */
/*             7  (PD)  :   Period doubling bifurcation. */
/*             8  (TR)  :   Bifurcation to an invariant torus. */
/*             9  (EP)  :   End point of branch, normal termination. */
/*            -9  (MX)  :   End point of branch, abnormal termination. */

/*  LAB        : The label of a special point. */
/*  PAR(ICP(1)): The principal parameter. */
/*  A          : The L2-norm of the solution vector, or other measure of 
*/
/*               the solution (see the user-supplied parameter IPLT). */
/*  MAX U(*)   : The maxima of the first few solution components. */
/*  PAR(ICP(*)): Further free parameters (if any). */

/* Local */

  
  ndim = iap->ndim;
  ips = iap->ips;
  irs = iap->irs;
  isw = iap->isw;
  iplt = iap->iplt;
  nmx = iap->nmx;
  npr = iap->npr;
  ndm = iap->ndm;
  itp = iap->itp;
  itpst = iap->itpst;
  ibr = iap->ibr;

  rl0 = rap->rl0;
  rl1 = rap->rl1;
  a0 = rap->a0;
  a1 = rap->a1;

  ntot = iap->ntot;
  ++ntot;
  iap->ntot = ntot;

/* ITP is set to 4 every NPR steps along a branch of solns and the entire 
*/
/* solution is written on unit 8. */

  if (npr != 0) {
    if (ntot % npr == 0 && itp % 10 == 0) {
      itp = itpst * 10 + 4;
    }
    iap->itp = itp;
  }

  /* Check whether limits of the bifurcation diagram have been reached : */

  iab = abs(iplt);
  if (iab == 0 || iab > ndm * 3) {
    amp = sqrt(rnrmsq(iap, &ndm, ndxloc, ups, dtm, thu));
  }
  if (iplt > 0 && iab <= ndm) {
    amp = rmxups(iap, ndxloc, &iab, ups);
  }
  if (iplt > ndm && iab <= (ndm * 2)) {
    amp = rintg(iap, ndxloc, iab - ndm, ups, dtm);
  }
  if (iplt > (ndm * 2) && iab <= ndm * 3) {
    integer tmp = iab - (ndm * 2);
    amp = rnrm2(iap, ndxloc, &tmp, ups, dtm);
  }
  if (iplt < 0 && iab <= ndm) {
    amp = rmnups(iap, ndxloc, &iab, ups);
  }

  rap->amp = amp;
  /* here is another place for byeauto
   call with iflag  */
  byeauto_(&iflag);
  istop = iap->istop;
  if (istop == 1) {
    /*        ** Maximum number of iterations reached somewhere. */
    itp = -9 - itpst * 10;
    iap->itp = itp;
  } else if (istop == -1) {
    /*        ** UZR endpoint */
    itp = itpst * 10 + 9;
    iap->itp = itp;
  } else {
    if (par[icp[0]] < rl0 || par[icp[0]] > rl1 || amp < a0 || amp > a1 || 
	ntot >= nmx || iflag==1 ) { /* more bye auto */
      istop = 1;
      iap->istop = istop;
      itp = itpst * 10 + 9;
      iap->itp = itp;
    }
  }

  /* All special points receive label: */

  labw = 0;
  if (itp % 10 != 0) {
    lab = iap->lab;
    ++lab;
    iap->lab = lab;
    labw = lab;
  }

  /* Compute maxima of solution components. */

  n2 = ndm;
  if (n2 > 7) {
    n2 = 7;
  }
  for (i = 0; i < n2; ++i) {
    itmp = i + 1;
    umx[i] = rmxups(iap, ndxloc, &itmp, ups);
  }

  /* compute max, min, etc of solutions   This is correct!! */
  for (i=0;i<NODE;i++){
    itmp=i+1;
    u_high[i]=rmxups(iap,ndxloc,&itmp,ups);
    u_low[i]=rmnups(iap,ndxloc,&itmp,ups);
    u_0[i]=ups[i*(*ndxloc)];
    u_bar[i]=rintg(iap, ndxloc, itmp, ups, dtm);
  }
    
  /* Determine stability, and write output on units 7 and 8. */

  ibrs = ibr;
  ntots = ntot;
  if (ips == 2 && abs(isw) != 2) {
    ibrs = -ibr;
    nins = iap->nins;
    if (nins == ndim) {
      ntots = -ntot;
    }
  }
  jtmp = NPARX;
  /* addbif max min  of variables & initial data 
    
   
  */
  addbif(iap,rap, ntots, ibrs, par,icp,labw,&amp, u_high, u_low, u_0, u_bar);  
  
  
  wrline(iap, rap, par, icp, &icp[jtmp], &ibrs, &ntots,
	 &labw, &amp, umx);

  /* Write plotting and restart data on unit 8. */

  if (itp % 10 != 0) {
    wrtbv8(iap, rap, par, icp, rldot, ndxloc, ups, 
	   udotps, tm, dtm);
  }

  return 0;
} /* stplbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
wrtbv8(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *rldot, integer *ndxloc, doublereal *ups, doublereal *udotps, doublereal *tm, doublereal *dtm)
{
  /* System generated locals */
  integer ups_dim1, udotps_dim1;

  /* Local variables */
  integer ndim, ncol, nfpr, ntpl, jtmp, mtot, ntot, ntst, i, j, k;
  doublereal t;
  integer k1, k2;
  doublereal rn;
  integer nrowpr, lab, ibr, nar, nrd, itp, isw;



  if(fp8_is_open==0) {
    fp8 = fopen(fort8,"w");
    fp8_is_open=1;
    if(fp8 == NULL) {
      fprintf(stderr,"Error:  Could not open fort.8\n");
      exit(1);
    }
  }

  /* Writes plotting and restart data on unit 8, viz.: */
  /* (1) data identifying the corresponding point on unit 7, */
  /* (2) the complete solution, */
  /* (3) the direction of the branch. */

  /* Specifically the following is written: */

/*  IBR   : The index of the branch. */
/*  NTOT  : The index of the point. */
/*  ITP   : The type of point (see STPLBV above). */
/*  LAB   : The label of the point. */
/*  NFPR : The number of free parameters used in the computation. */
/*  ISW   : The value of ISW used in the computation. */
/*  NTPL  : The number of points in the time interval [0,1] for which */
/*          solution values are wriiten. */
/*  NAR   : The number of values written per point. */
/*          (NAR=NDIM+1, since T and U(i), i=1,..,NDIM are written). */
/*  NROWPR: The number of lines printed following the identifying line */
/*          and before the next data set or the end of the file. */
/*          (Used for quickly skipping a data set when searching). */
/*  NTST  : The number of time intervals used in the discretization. */
/*  NCOL  : The number of collocation points used. */
/*  NPARX : The dimension of the par array (and the number of */
/*          number of values in the parameter block).*/

/*  Following the above described identifying line there are NTPL lines */
/* containing : */
/*     T , U-1(T) , U-2(T) , ... , U-NDIM(T), */
/* where NDIM is the dimension of the system of differential equations. */

/*  Following this is a line which lists the active parameters (from ICP) */

/* Following this is a line containing */
/*    RL-dot(i) , i=1,NFPR, */

/* and following this are NTPL lines each containing */
/*    U-dot-1(T), U-dot-2(T), ... , U-dot-NDIM(T). */

/* Finally the parameter values PAR(i) , i=1,NPARX, are written. */

/*  Above, RL-dot(.) and U-dot(.) specify the direction of the branch. */


  udotps_dim1 = *ndxloc;
  ups_dim1 = *ndxloc;

  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  isw = iap->isw;
  itp = iap->itp;
  nfpr = iap->nfpr;
  ibr = iap->ibr;
  ntot = iap->ntot;
  lab = iap->lab;

  /* Write information identifying the solution : */

  ntpl = ncol * ntst + 1;
  nar = ndim + 1;
  nrd = ndim / 7 + 2 + (ndim - 1) / 7;
  jtmp = NPARX;
  nrowpr = nrd * (ncol * ntst + 1) + (nfpr - 1) / 7 + 1 + (jtmp - 1) / 7 + 
    1 + (nfpr - 1) / 20 + 1;

  if (iap->mynode > 0) {
    return 0;
  }

  mtot = ntot % 10000;
  fprintf(fp8,"%5ld",ibr);
  fprintf(fp8,"%5ld",mtot);
  fprintf(fp8,"%5ld",itp);
  fprintf(fp8,"%5ld",lab);
  fprintf(fp8,"%5ld",nfpr);
  fprintf(fp8,"%5ld",isw);
  fprintf(fp8,"%5ld",ntpl);
  fprintf(fp8,"%5ld",nar);
  fprintf(fp8,"%5ld",nrowpr);
  fprintf(fp8,"%5ld",ntst);
  fprintf(fp8,"%5ld",ncol);
  fprintf(fp8,"%5d\n",NPARX);

/* Write the entire solution on unit 8 : */

  for (j = 0; j < ntst; ++j) {
    rn = 1. / ncol;
    for (i = 0; i < ncol; ++i) {
      k1 = i * ndim;
      k2 = (i + 1) * ndim;
      t = tm[j] + i * rn * dtm[j];
      fprintf(fp8,"    %19.10E",t);
      for (k = k1; k < k2; ++k) {
	if((k+1-k1)%7==0)
	  fprintf(fp8,"\n    ");
	fprintf(fp8,"%19.10E",ARRAY2D(ups, j, k));
      }
      fprintf(fp8,"\n");	 

    }
  }
  fprintf(fp8,"    %19.10E",tm[ntst]);
  for (i = 0; i < ndim; ++i) {
    if((i+1)%7==0)
      fprintf(fp8,"\n    ");
    fprintf(fp8,"%19.10E",ARRAY2D(ups, ntst, i));
  }
  fprintf(fp8,"\n");  

/* Write the free parameter indices: */
  for (i = 0; i < nfpr; ++i) {
    fprintf(fp8,"%5ld",icp[i]);
  }
  fprintf(fp8,"\n");  
    

/* Write the direction of the branch: */
  fprintf(fp8,"    ");
  for (i = 0; i < nfpr; ++i) {
    if((i>0)&&((i)%7==0))
      fprintf(fp8,"\n    ");
    fprintf(fp8,"%19.10E",rldot[i]);
  }
  fprintf(fp8,"\n");  

  for (j = 0; j < ntst; ++j) {
    for (i = 0; i < ncol; ++i) {
      k1 = i * ndim;
      k2 = (i + 1)* ndim;

      fprintf(fp8,"    ");
      for (k = k1; k < k2; ++k) {
	if((k!=k1)&&((k-k1)%7==0))
	  fprintf(fp8,"\n    ");
	fprintf(fp8,"%19.10E",ARRAY2D(udotps, j, k));
      }
      fprintf(fp8,"\n");
    }
  }
  fprintf(fp8,"    ");

  for (k = 0; k < ndim; ++k) {
    if((k!=0)&&(k%7==0))
      fprintf(fp8,"\n    ");
    fprintf(fp8,"%19.10E",ARRAY2D(udotps, ntst, k));
  }
  fprintf(fp8,"\n");

/* Write the parameter values. */

  fprintf(fp8,"    ");
  for (i = 0; i < NPARX; ++i) {
    if((i>0)&&(i%7==0))
      fprintf(fp8,"\n    ");
    fprintf(fp8,"%19.10E",par[i]);
  }
  fprintf(fp8,"\n");
  fflush(fp8);
  return 0;
} /* wrtbv8_ */


/*     ---------- ------ */
/* Subroutine */ int 
wrtbv9(iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *rlcur, integer *ndxloc, doublereal *ups, doublereal *tm, doublereal *dtm, doublereal *thl, doublereal *thu)
{
  
  /* System generated locals */
  integer ups_dim1;
  
  
  
  /* Local variables */
  integer ndim, ncol, nfpr, iplt, mtot, ntot, ntst, i, j, k;
  doublereal t;
  integer nfprp, k1, nitps, k2;
  doublereal ds, rn;
  
  integer iab, iid, ibr, ndm;
  doublereal amp=0.0;
  
  /* Writes additional output on unit 9. */

  ups_dim1 = *ndxloc;
  
  
  ndim = iap->ndim;
  ntst = iap->ntst;
  ncol = iap->ncol;
  iplt = iap->iplt;
  iid = iap->iid;
  ndm = iap->ndm;
  nfpr = iap->nfpr;
  ibr = iap->ibr;
  nitps = iap->nit;
  ntot = iap->ntot;
  ds = rap->ds;
  
  iab = abs(iplt);
  if (iab == 0 || iab > ndim) {
    amp = sqrt(rnrmsq(iap, &ndm, ndxloc, ups, dtm, thu));
  }
  if (iplt > 0 && iab <= ndim) {
    amp = rmxups(iap, ndxloc, &iab, ups);
  }
  if (iplt < 0 && iab <= ndim) {
    amp = rmnups(iap, ndxloc, &iab, ups);
  }
  rap->amp = amp;
  if (iid >= 2) {
    if (nfpr <= 5) {
      nfprp = nfpr;
    } else {
      nfprp = 5;
    }
    if (iap->mynode == 0) {
      if (nitps == 0 || iid >= 3) {
	fprintf(fp9,"========================================");
	fprintf(fp9,"========================================\n");
	fprintf(fp9,"  BR    PT  IT\n");	
	
      }
      mtot = (ntot + 1) % 10000;
      fprintf(fp9,"%4li%6li%4li    %14.6E%14.6E\n",ibr,mtot,nitps,rlcur[0],amp);	
      
    }
  }
  
  if (iid >= 5 && iap->mynode == 0) {
    fprintf(fp9," UPS :\n");	
    for (j = 0; j < ntst; ++j) {
      rn = 1. / ncol;
      for (i = 0; i < ncol; ++i) {
	t = tm[j] + i * rn * dtm[j];
	k1 = i * ndim;
	k2 = (i + 1) * ndim;
	fprintf(fp9," %14.6E",t);	
	for (k = k1; k < k2; ++k) {
	  if((k+1-k1)%7==0)
	    fprintf(fp9,"\n ");
	  fprintf(fp9," %14.6E",ARRAY2D(ups, j, k));	
	}
	fprintf(fp9,"\n");	
	
      }
    }
    fprintf(fp9," %14.6E",tm[ntst]);	
    for (i = 0; i < ndim; ++i) {
      if((i+1)%7==0)
	fprintf(fp9,"\n ");
      fprintf(fp9," %14.6E",ARRAY2D(ups, ntst, i));	
    }
    fprintf(fp9,"\n");	
  }
  return 0;
} /* wrtbv9_ */


/*     ---------- ------ */
/* Subroutine */ int 
pvlsae(iap_type *iap, rap_type *rap, doublereal *u, doublereal *par)
{

  integer ndm;

  
  setpae(iap, rap);
  ndm = iap->ndm;
  pvls(ndm, u, par);

  return 0;
} /* pvlsae_ */


/*     ---------- ------ */
/* Subroutine */ int 
pvlsbv(iap_type *iap, rap_type *rap, integer *icp, doublereal *dtm, integer *ndxloc, doublereal *ups, integer *ndim, doublereal *p0, doublereal *p1, doublereal *par)
{
  integer ndm;

  setpbv(iap, rap, dtm);
  ndm = iap->jac;
  pvls(ndm, ups, par);

  return 0;
} /* pvlsbv_ */


/*     ---------- ------ */
/* Subroutine */ int 
setpae(iap_type *iap, rap_type *rap)
{
  global_parameters.iav = iap;
  global_parameters.rav = rap;

  return 0;
} /* setpae_ */


/*     ---------- ------ */
/* Subroutine */ int 
setpbv(iap_type *iap, rap_type *rap, doublereal *dtm)
{
  global_parameters.iav = iap;
  global_parameters.rav = rap;

  global_parameters.dtv = dtm;
  return 0;
} /* setpbv_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*          System Dependent Subroutines for Timing AUTO */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

double time_start(void) {
  struct timeval time;
  double seconds,microseconds;
  gettimeofday(&time,NULL);
  seconds = (double)time.tv_sec;
  microseconds = (double)time.tv_usec;
  return seconds + microseconds/1e6;
}
double time_end(double start) {
  struct timeval time;
  double seconds,microseconds;
  gettimeofday(&time,NULL);
  seconds = (double)time.tv_sec;
  microseconds = (double)time.tv_usec;
  return (seconds + microseconds/1e6)-start;
}


doublereal 
getp(char *code, integer *ic, doublereal *ups, integer code_len)
{
  /* System generated locals */
  doublereal ret_val=0.0;

    /* Local variables */
  integer ntst;

  integer nxloc;

  integer ips;



  
  nxloc = global_parameters.iav->ntst + 1;

  ips = global_parameters.iav->ips;

  ntst = global_parameters.iav->ntst;


  if (abs(ips) <= 1 || ips == 5) {

    if (strcmp(code, "NRM") == 0 || strcmp(code, "nrm") == 0) {
      ret_val = fabs(ups[*ic - 1]);
    } else if (strcmp(code, "INT") == 0 || strcmp(code, "int") == 0) {
      ret_val = ups[*ic - 1];

    } else if (strcmp(code, "MAX") == 0 || strcmp(code, "max") == 0) {

      ret_val = ups[*ic - 1];

    } else if (strcmp(code, "MIN") == 0 || strcmp(code, "min") == 0) {

      ret_val = ups[*ic - 1];

    } else if (strcmp(code, "BV0") == 0 || strcmp(code, "bv0") == 0) {

      ret_val = ups[*ic - 1];

    } else if (strcmp(code, "BV1") == 0 || strcmp(code, "bv1") == 0) {

      ret_val = ups[*ic - 1];

    } else if (strcmp(code, "STP") == 0 || strcmp(code, "stp") == 0) {

      ret_val = global_parameters.rav->dsold;

    } else if (strcmp(code, "FLD") == 0 || strcmp(code, "fld") == 0) {

      ret_val = global_parameters.rav->fldf;

    } else if (strcmp(code, "HBF") == 0 || strcmp(code, "hbf") == 0) {

      ret_val = global_parameters.rav->hbff;

    } else if (strcmp(code, "BIF") == 0 || strcmp(code, "bif") == 0) {

      ret_val = global_parameters.rav->biff;

    } else if (strcmp(code, "SPB") == 0 || strcmp(code, "spb") == 0) {

      ret_val = (double)0.0;

    }

  } else {

    if (strcmp(code, "NRM") == 0 || strcmp(code, "nrm") == 0) {

      ret_val = rnrm2(global_parameters.iav, &nxloc, ic, ups, global_parameters.dtv);

    } else if (strcmp(code, "INT") == 0 || strcmp(code, "int") == 0) {

      ret_val = rintg(global_parameters.iav, &nxloc, *ic, ups, global_parameters.dtv);

    } else if (strcmp(code, "MAX") == 0 || strcmp(code, "max") == 0) {

      ret_val = rmxups(global_parameters.iav, &nxloc, ic, ups);

    } else if (strcmp(code, "MIN") == 0 || strcmp(code, "min") == 0) {

      ret_val = rmnups(global_parameters.iav, &nxloc, ic, ups);

    } else if (strcmp(code, "BV0") == 0 || strcmp(code, "bv0") == 0) {

      ret_val = ups[(*ic-1) * (global_parameters.iav->ntst + 1)];

    } else if (strcmp(code, "BV1") == 0 || strcmp(code, "bv1") == 0) {

      ret_val = ups[ntst + (*ic-1) * (global_parameters.iav->ntst + 1)];

    } else if (strcmp(code, "STP") == 0 || strcmp(code, "stp") == 0) {

      ret_val = global_parameters.rav->dsold;

    } else if (strcmp(code, "FLD") == 0 || strcmp(code, "fld") == 0) {

      ret_val = global_parameters.rav->fldf;

    } else if (strcmp(code, "HBF") == 0 || strcmp(code, "hbf") == 0) {

      ret_val = (double)0.;

    } else if (strcmp(code, "BIF") == 0 || strcmp(code, "bif") == 0) {

      ret_val = global_parameters.rav->biff;

    } else if (strcmp(code, "SPB") == 0 || strcmp(code, "spb") == 0) {

      ret_val = global_parameters.rav->spbf;

    }

  }


  return ret_val;

}

void allocate_global_memory(const iap_type iap) {
    free(global_scratch.dfu);
    free(global_scratch.dfp);
    free(global_scratch.uu1);
    free(global_scratch.uu2);
    free(global_scratch.ff1);
    free(global_scratch.ff2);
    global_scratch.dfu = (doublereal *)malloc(sizeof(doublereal)*(iap.ndim)*(iap.ndim));
    global_scratch.dfp = (doublereal *)malloc(sizeof(doublereal)*(iap.ndim)*NPARX);
    global_scratch.uu1 = (doublereal *)malloc(sizeof(doublereal)*(iap.ndim));
    global_scratch.uu2 = (doublereal *)malloc(sizeof(doublereal)*(iap.ndim));
    global_scratch.ff1 = (doublereal *)malloc(sizeof(doublereal)*(iap.ndim));
    global_scratch.ff2 = (doublereal *)malloc(sizeof(doublereal)*(iap.ndim));

    free(global_rotations.nrtn);
    global_rotations.nrtn = (integer *)malloc(sizeof(integer)*(iap.nbc));
}



