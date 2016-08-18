#include "do_fit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "browse.h"
#include "delay_handle.h"
#include "derived.h"
#include "form_ode.h"
#include "ggets.h"
#include "integrate.h"
#include "load_eqn.h"
#include "do_fit.h"
#include "numerics.h"
#include "parserslow.h"
#include "pop_list.h"
#include "solver/cv2.h"
#include "solver/dormpri.h"
#include "solver/gear.h"
#include "solver/rosenbrock.h"
#include "solver/stiff.h"
#include "storage.h"
#include "util/matrixalgebra.h"


/* --- Types --- */
typedef struct {
	char file[25];
	char varlist[25],collist[25];
	char parlist1[25],parlist2[25];
	int dim,npars,nvars,npts,maxiter;
	int icols[50],ipar[50],ivar[50];
	double tol,eps;
} FITINFO;


/* --- Forward declarations --- */
int (*solver)();
static void get_fit_info(double *y, double *a, double *t0, int *flag, double eps, double *yfit, double **yderv, int npts, int npars, int nvars, int *ivar, int *ipar);
static int get_fit_params(void);
static int marlevstep(double *t0, double *y0, double *y, double *sig, double *a, int npts, int nvars, int npars, int *ivar, int *ipar, double *covar, double *alpha, double *chisq, double *alambda, double *work, double **yderv, double *yfit, double *ochisq, int ictrl, double eps);
static int mrqcof(double *t0, double *y0, double *y, double *sig, double *a, int npts, int nvars, int npars, int *ivar, int *ipar, double *alpha, double *chisq, double *beta, double **yderv, double *yfit, double eps);
static void parse_collist(char *collist, int *icols, int *n);
static void parse_varlist(char *varlist, int *ivars, int *n);
static void parse_parlist(char *parlist, int *ipars, int *n);
static void print_fit_info(void);
static int run_fit(char *filename, int npts, int npars, int nvars, int maxiter, int ndim, double eps, double tol, int *ipar, int *ivar, int *icols, double *y0, double *a, double *yfit);


/* --- Data --- */
static FITINFO fin;


/* --- Functions --- */
void init_fit_info(void) {
	fin.tol=.001;
	fin.eps=1e-5;
	fin.dim=0;
	fin.npars=0;
	fin.nvars=0;
	fin.varlist[0]=0;
	fin.collist[0]=0;
	fin.parlist1[0]=0;
	fin.parlist2[0]=0;
	fin.npts=0;
	fin.maxiter=20;
	fin.file[0]=0;
}


void test_fit(void) {
	double *yfit,a[1000],y0[1000];
	int nvars = 0,npars,i,ok;
	char collist[30],parlist1[30],parlist2[30],varlist[30];
	fin.nvars=0;
	fin.npars=0;
	if(get_fit_params()==0)
		return;

	sprintf(collist,"%s",fin.collist);
	sprintf(varlist,"%s",fin.varlist);
	sprintf(parlist1,"%s",fin.parlist1);
	sprintf(parlist2,"%s",fin.parlist2);

	parse_collist(collist,fin.icols,&nvars);

	if(nvars<=0) {
		err_msg("No columns...");
		return;
	}
	fin.nvars=nvars;
	nvars=0;
	parse_varlist(varlist, fin.ivar, &nvars);

	if(fin.nvars!=nvars) {
		err_msg(" # columns != # fitted variables");
		return;
	}
	npars=0;
	parse_parlist(parlist1,fin.ipar,&npars);

	parse_parlist(parlist2,fin.ipar,&npars);

	if(npars<=0) {
		err_msg(" No parameters!");
		return;
	}
	fin.npars=npars;
	for(i=0;i<npars;i++) {
		if(fin.ipar[i]>=0) {
			if(fin.ipar[i]>=NODE) {
				err_msg(" Cant vary auxiliary/markov variables! ");
				return;
			}
		}
	}
	for(i=0;i<nvars;i++) {
		if(fin.icols[i]<2) {
			err_msg(" Illegal column must be >= 2");
			return;
		}
		if(fin.ivar[i]<0 || fin.ivar[i]>=NODE) {
			err_msg(" Fit only to variables! ");
			return;
		}
	}
	yfit=(double *)malloc(fin.npts*fin.nvars*sizeof(double));
	for(i=0;i<NODE;i++) {
		y0[i]=last_ic[i];
	}
	for(i=0;i<fin.npars;i++) {
		if(fin.ipar[i]<0) {
			a[i]=constants[-fin.ipar[i]];
		} else {
			a[i]=last_ic[fin.ipar[i]];
		}
	}

	print_fit_info();
	plintf(" Running the fit...\n");
	ok=run_fit(fin.file, fin.npts,fin.npars,fin.nvars,fin.maxiter,fin.dim,
			   fin.eps,fin.tol,
			   fin.ipar,fin.ivar,fin.icols,
			   y0,a,yfit);

	free(yfit);
	if(ok==0)
		return;

	/* get the latest par values ...  */

	for(i=0;i<npars;i++) {
		if(fin.ipar[i]<0) {
			constants[-fin.ipar[i]]=a[i];
		} else {
			last_ic[fin.ipar[i]]=a[i];
		}
	}
}


/* --- Static functions --- */
/*
  y     initial condition
  a     initial guesses for the parameters
  t0    vector of output times
  flag  1 for success  0 for failure
  eps   derivative step
  yfit  has y[i1](t0),...,y[im](t0), ..., y[i1](tn),...,y[im](tn)
		which are the values of the test functions at the npts
	times.  yfit is (npts)*nvars int
  yderv[npar][nvars*(npts)] is the derivative of yfit with rrspect
		to the parameter

  npts   is the number of times to be fitted
  npars    the number of parameters
  nvars    the number of variables
  ipar     the vector of parameters  negative are constants
		   positive are initial data
  ivar     the vector of variables

 */
static void get_fit_info(double *y, double *a, double *t0, int *flag, double eps,
						 double *yfit, double **yderv, int npts, int npars,
						 int nvars, int *ivar, int *ipar) {
	int i,iv,ip,istart=1,j,k,l,k0,ok;
	double yold[MAXODE],dp;
	double par;
	*flag=0;
	/* set up all initial data and parameter guesses  */
	for(l=0;l<npars;l++) {
		ip=ipar[l];
		if(ip<0) {
			constants[-ip]=a[l];
		} else {
			y[ip]=a[l];
		}
	}
	for(i=0;i<NODE;i++) {
		yold[i]=y[i];
		/*  plintf(" init y[%d]=%g \n",i,y[i]); */
	}
	if(DelayFlag) {
		/* restart initial data */
		if(do_init_delay(DELAY)==0)
			return;
	}
	evaluate_derived();
	/*   This gets the values at the desired points  */
	for(i=0;i<nvars;i++) {
		iv=ivar[i];
		yfit[i]=y[iv];
	}
	for(k=1;k<npts;k++) {
		k0=k*nvars;
		ok=one_step_int(y,t0[k-1],t0[k],&istart);
		if(ok==0) {
			for(i=0;i<NODE;i++) {
				y[i]=yold[i];
			}
			return;
		}

		for(i=0;i<nvars;i++) {
			iv=ivar[i];
			yfit[i+k0]=y[iv];
		}
	}
	if(METHOD==METHOD_CVODE)
		end_cvode();
	/*  Now we take the derivatives !!   */
	for(l=0;l<npars;l++) {
		istart=1;
		/* set up all the initial conditions   */
		for(j=0;j<nvars;j++) {
			/* no dependence on initial data ... */
			yderv[l][j]=0.0;
		}
		for(i=0;i<NODE;i++) {
			y[i]=yold[i];
		}
		ip=ipar[l];
		if(ip<0) {
			par=constants[-ip];
			dp=eps*MAX(eps,fabs(par));
			constants[-ip]=par+dp;
		} else {
			par=yold[ip];
			dp=eps*MAX(eps,fabs(par));
			y[ip]=par+dp;
			for(j=0;j<nvars;j++) {
				if(ip==ivar[j]) {
					/* ... except for those ICs that can vary */
					yderv[l][j]=1.0;
				}
			}
		}
		if(DelayFlag) {
			/* restart initial data */
			if(do_init_delay(DELAY)==0)
				return;
		}
		evaluate_derived();
		/* now loop through all the points */
		for(k=1;k<npts;k++) {
			k0=k*nvars;
			ok=one_step_int(y,t0[k-1],t0[k],&istart);
			if(ok==0) {
				for(i=0;i<NODE;i++) {
					y[i]=yold[i];
				}
				return;
			}
			for(i=0;i<nvars;i++) {
				iv=ivar[i];
				yderv[l][i+k0]=(y[iv]-yfit[i+k0])/dp;
			}
		}
		/* Now return the parameter to its old value */
		if(ip<0) {
			constants[-ip]=par;
		}
		evaluate_derived();
		if(METHOD==METHOD_CVODE)
			end_cvode();
	}
	*flag=1;
	for(i=0;i<NODE;i++) {
		y[i]=yold[i];
	}
}


static int get_fit_params(void) {
	static char *n[]={"File", "Fitvar","Params","Tolerance","Npts",
					  "NCols","To Col","Params","Epsilon","Max iter"};
	int status;
	char values[10][MAX_LEN_SBOX];
	sprintf(values[0],"%s",fin.file);
	sprintf(values[1],"%s",fin.varlist);
	sprintf(values[2],"%s",fin.parlist1);
	sprintf(values[3],"%g",fin.tol);
	sprintf(values[4],"%d",fin.npts);
	sprintf(values[5],"%d",fin.dim);
	sprintf(values[6],"%s",fin.collist);
	sprintf(values[7],"%s",fin.parlist2);
	sprintf(values[8],"%g",fin.eps);
	sprintf(values[9],"%d",fin.maxiter);
	status=do_string_box(10,5,2,"Fit",n,values,45);
	if(status!=0) {
		fin.tol=atof(values[3]);
		fin.npts=atoi(values[4]);
		fin.dim=atoi(values[5]);
		fin.eps=atof(values[8]);
		fin.maxiter=atoi(values[9]);
		sprintf(fin.file,"%s",values[0]);
		sprintf(fin.varlist,"%s",values[1]);
		sprintf(fin.parlist1,"%s",values[2]);
		sprintf(fin.collist,"%s",values[6]);
		sprintf(fin.parlist2,"%s",values[7]);
		return 1;
	}
	return 0;
}


/*   One step of Levenberg-Marquardt
nvars  the number of variables to fit
ivar   their indices
npars  the number of parameters to alter
ipar   their indices
npts   the number of times
ictrl  0 to start  1 to continue  2 to finish up

t0  the npts times
y0  the NODE initial data
y   the (npts)*nvars data points to fit
sig  the nvars  weights
a  the npars initial guesses of the things to be fit
chisq  the chisquare
alpha is work array and also the curvature matrix
covar is the covariance matrix (npars x npars)
alambda is control of step size; start negative 0 to get final value
work is an array of size npar*4+npar*npar
yderv is  npar x (nptts+1)*nvars
yfit  is  (npts)*nvars  on each completed step it has the fitted soln
eps   control numerical derivative
sigma  weights on nvars
*/
static int marlevstep(double *t0, double *y0, double *y, double *sig, double *a,
					  int npts, int nvars, int npars, int *ivar, int *ipar,
					  double *covar, double *alpha, double *chisq, double *alambda,
					  double *work, double **yderv, double *yfit, double *ochisq,
					  int ictrl, double eps) {
	int i,j,k,ierr,ipivot[1000];

	double *da,*atry,*beta,*oneda;
	da=work;
	atry=work+npars;
	beta=work+2*npars;
	oneda=work+3*npars;

	if(ictrl==0) {
		*alambda=.001;
		if(mrqcof(t0,y0,y,sig,a,npts,nvars,npars,
				  ivar,ipar,alpha,chisq,beta,
				  yderv,yfit,eps)==0) {
			return 0;
		}
		for(i=0;i<npars;i++) {
			atry[i]=a[i];
		}
		*ochisq=(*chisq);
	}
	for(j=0;j<npars;j++) {
		for(k=0;k<npars;k++) {
			covar[j+k*npars]=alpha[j+k*npars];
		}
		covar[j+j*npars]=alpha[j+j*npars]*(1+(*alambda));
		oneda[j]=beta[j];
	}
	sgefa(covar,npars,npars,ipivot,&ierr);
	if(ierr!=-1) {
		err_msg(" Singular matrix encountered...");
		return 0;
	}

	sgesl(covar,npars,npars,ipivot,oneda,0);
	for(j=0;j<npars;j++) {
		da[j]=oneda[j];
		/* plintf(" da[%d]=%g \n",j,da[j]); */
	}
	if(ictrl==2) {  /* all done invert alpha to get the covariance */
		for(j=0;j<(npars*npars);j++) {
			alpha[j]=covar[j];
		}
		for(j=0;j<npars;j++) {
			for(k=0;k<npars;k++) {
				oneda[k]=0.0;
			}
			oneda[j]=1.0;
			sgesl(alpha,npars,npars,ipivot,oneda,0);
			for(k=0;k<npars;k++) {
				covar[j+k*npars]=oneda[k];
			}
		}
		return 1;
	}
	for(j=0;j<npars;j++) {
		atry[j]=a[j]+da[j];
		/*    plintf(" aold[%d]=%g anew[%d]=%g \n",
	   j,a[j],j,atry[j]); */
	}
	if(mrqcof(t0,y0,y,sig,atry,npts,nvars,npars,
			  ivar,ipar,covar,chisq,da,
			  yderv,yfit,eps)==0)
		return 0;

	if(*chisq<*ochisq) {
		/* *ochisq=*chisq; */
		*alambda *= 0.1;
		for(j=0;j<npars;j++) {
			for(k=0;k<npars;k++) {
				alpha[j+k*npars]=covar[j+k*npars];
			}
			beta[j]=da[j];
			a[j]=atry[j];
		}
	} else {
		*alambda *= 10.0;
		/* *chisq=*ochisq; */
	}
	return 1;
}


static int mrqcof(double *t0, double *y0, double *y, double *sig, double *a,
				  int npts, int nvars, int npars, int *ivar, int *ipar,
				  double *alpha, double *chisq, double *beta, double **yderv,
				  double *yfit, double eps) {
	int flag,i,j,k,l,k0;
	double sig2i,dy,wt;

	get_fit_info(y0,a,t0,&flag,eps,yfit,yderv,npts,npars,nvars,ivar,ipar);
	if(flag==0) {
		err_msg(" Integration error ...\n");
		return 0;
	}
	for(i=0;i<npars;i++) {
		beta[i]=0.0;
		for(j=0;j<npars;j++) {
			alpha[i+j*npars]=0.0;
		}
	}
	*chisq=0.0;
	for(i=0;i<nvars;i++) {
		sig2i=1.0/(sig[i]*sig[i]);
		for(k=0;k<npts;k++) {
			k0=k*nvars+i;
			dy=y[k0]-yfit[k0];
			/*           plintf(" i=%d k=%d dy = %f \n",i,k,dy); */
			for(j=0;j<npars;j++) {
				wt=yderv[j][k0]*sig2i;
				for(l=0;l<npars;l++) {
					alpha[j+l*npars] += wt*yderv[l][k0];
				}
				beta[j] += dy*wt;
			}
			(*chisq) += dy*dy*sig2i;

			/* the last loop could be halved because of symmetry, but I am lazy
   and this is an insignificiant amount of the CPU time since
  the evaluation step is really where all the time is used
*/
		}
	}
	return 1;
}


/* gets a list of the data columns to use ... */
static void parse_collist(char *collist, int *icols, int *n) {
	char *item, *toksave;
	int v,i=0;

	item=strtok_r(collist, " ,", &toksave);

	if(item[0]==0) {
		return;
	}
	v=atoi(item);
	icols[i]=v;
	i++;
	while((item=strtok_r(NULL, " ,", &toksave))!=NULL) {
		v=atoi(item);
		icols[i]=v;
		i++;
	}
	*n=i;
}


static void parse_varlist(char *varlist, int *ivars, int *n) {
	char *item, *toksave;
	int v,i=0;

	item=strtok_r(varlist," ,", &toksave);
	if(item[0]==0) {
		return;
	}
	find_variable(item,&v);

	if(v<=0) {
		return;
	}
	ivars[i]=v-1;
	i++;
	while((item=strtok_r(NULL, " ,", &toksave))!=NULL) {
		find_variable(item,&v);
		if(v<=0) {
			return;
		}
		ivars[i]=v-1;
		i++;
	}
	*n=i;
}


static void parse_parlist(char *parlist, int *ipars, int *n) {
	char *item, *toksave;
	int v,i=0;
	int j;
	for(j=0;j<strlen(parlist);j++) {
		if(parlist[j]!=' ') {
			break;
		}
	}
	if(j==strlen(parlist)) {
		return;
	}
	if(strlen(parlist)==0) {
		return;
	}
	item=strtok_r(parlist," ,", &toksave);
	if(item[0]==0L) {
		return;
	}

	find_variable(item,&v);
	if(v>0) {
		ipars[i+*n]=v-1;
		i++;
	} else {
		v=get_param_index(item);
		if(v<=0) {
			return;
		}
		ipars[i+*n]=-v;
		i++;
	}
	while((item=strtok_r(NULL, " ,", &toksave))!=NULL) {
		find_variable(item,&v);
		if(v>0) {
			ipars[i+*n]=v-1;
			i++;
		} else {
			v=get_param_index(item);
			if(v<=0) {
				return;
			}
			ipars[i+*n]=-v;
			i++;
		}
	}
	*n=*n+i;
}


static void print_fit_info(void) {
	int i;
	plintf("dim=%d maxiter=%d npts=%d file=%s tol=%g eps=%g\n",
		   fin.dim,fin.maxiter,fin.npts,fin.file,fin.tol,fin.eps);

	for(i=0;i<fin.nvars;i++) {
		plintf(" variable %d to col %d \n",
			   fin.ivar[i],fin.icols[i]);
	}
	for(i=0;i<fin.npars;i++) {
		plintf(" P[%d]=%d \n",i,fin.ipar[i]);
	}
}

/*
   filename is where the data file is -- it is of the form:
   t1 y11 y12 .... y1m
   t2 ....
   ...
   tn yn1 ....     ynm
   icols gives the dependent variable columns -- we assume first col
   is the times
   ndim is the number of y-pts in the a row

*/
static int run_fit(char *filename, int npts, int npars, int nvars, int maxiter, int ndim,
				   double eps, double tol, int *ipar, int *ivar, int *icols, double *y0,
				   double *a, double *yfit) {
	double *t0,*y,sig[MAXODE],*covar,*alpha,chisq,ochisq,alambda,**yderv,*work;
	int i,j,k,ioff,ictrl=0,ok;
	FILE *fp;
	int niter=0,good_flag=0;
	double tol10=10*tol;
	double t,ytemp[MAXODE];

	if((fp=fopen(filename,"r"))==NULL) {
		err_msg("No such file...");
		return 0;
	}
	t0=(double *)malloc((npts+1)*sizeof(double));
	y=(double *)malloc((npts+1)*nvars*sizeof(double));
	/* load up the data to fit   */

	for(i=0;i<npts;i++) {
		if(fscanf(fp,"%lg ",&t)<1) {
			plintf("Error reading file %s", fp);
		}

		for(j=0;j<ndim-1;j++) {
			if(fscanf(fp,"%lg ",&ytemp[j])<1) {
				plintf("Error reading file %s", fp);
			}
		}
		t0[i]=t;

		ioff=nvars*i;
		for(k=0;k<nvars;k++) {
			y[ioff+k]=ytemp[icols[k]-2];
		}
	}
	plintf(" Data loaded ... %f %f ...  %f %f \n",
		   y[0],y[1],y[npts*nvars-2],y[npts*nvars-1]);

	work=(double *)malloc(sizeof(double)*(4*npars+npars*npars));
	yderv=(double **)malloc(npars*sizeof(double *));
	for(i=0;i<npars;i++) {
		yderv[i]=(double *)malloc((npts+1)*nvars*sizeof(double));
	}
	for(i=0;i<nvars;i++) {
		sig[i]=1.0;
	}
	covar=(double *)malloc(npars*npars*sizeof(double));
	alpha=(double *)malloc(npars*npars*sizeof(double));

	while(good_flag<3) {  /* take 3 good steps after convergence  */

		ok=marlevstep(t0,y0,y,sig,a,npts,nvars,npars,
					  ivar,ipar,covar,alpha,&chisq,&alambda,work,
					  yderv,yfit,&ochisq,ictrl,eps);
		niter++;
		plintf(" step %d is %d  -- lambda= %g  chisq= %g oldchi= %g\n",
			   niter,ok,alambda,chisq,ochisq);
		plintf(" params: ");
		for(i=0;i<npars;i++) {
			plintf(" %g ",a[i]);
		}
		plintf("\n");
		if((ok==0) || (niter>=maxiter))
			break;
		if(ochisq>chisq) {
			if(((ochisq-chisq)<tol10) || (((ochisq-chisq)/MAX(1.0,chisq))<tol)) {
				good_flag++;
				niter--;  /* compensate for good stuff ... */
			}
			ochisq=chisq;
		} else {
			chisq=ochisq;
		}
		ictrl=1;
	}

	if(ok==0) {
		err_msg("Error in step...");

		free(work);
		for(i=0;i<npars;i++) {
			free(yderv[i]);
		}
		free(yderv);
		free(alpha);
		free(covar);
		free(t0);
		free(y);
		return 0;
	}
	if(niter>=maxiter) {
		err_msg("Max iterations exceeded...");

		free(work);
		for(i=0;i<npars;i++) {
			free(yderv[i]);
		}
		free(yderv);
		free(alpha);
		free(covar);
		free(t0);
		free(y);

		return 1;
	}
	ictrl=2;
	marlevstep(t0,y0,y,sig,a,npts,nvars,npars,
			   ivar,ipar,covar,alpha,&chisq,&alambda,work,
			   yderv,yfit,&ochisq,ictrl,eps);
	err_msg(" Success! ");
	/* have the covariance matrix -- so what?   */
	plintf(" covariance: \n");
	for(i=0;i<npars;i++) {
		for(j=0;j<npars;j++) {
			plintf(" %g ",covar[i+npars*j]);
		}
		plintf("\n");
	}

	free(work);
	for(i=0;i<npars;i++) {
		free(yderv[i]);
	}
	free(yderv);
	free(alpha);
	free(covar);
	free(t0);
	free(y);

	return 1;
}
