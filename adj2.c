#include "adj2.h"
#include "phase_red.h"
#include "my_rhs.h"
#include "pop_list.h"
#include "browse.h"
#include "ggets.h"
#include "do_fit.h"
#include "lunch-new.h"
#include "gear.h"
#include "integrate.h"
#include "parserslow.h"


/* #include "stdfix.h"  */
#include "stdfix.h"
#include <string.h>
/*
  this has a bunch of numerical routines
  averaging
  adjoints
  transpose
  maximal liapunov exponent
  
*/

#include <stdio.h>
#include <math.h>

#include "xpplim.h"


#define MAX_LEN_SBOX 25
#define READEM 1

double evaluate();
double ndrand48();

extern double MyData[MAXODE];
extern int (*rhs)();
extern float **storage;
extern int storind,FOUR_HERE;
extern int NODE,INFLAG,NEQ,NJMP,FIX_VAR,NMarkov,nvec;
extern double TEND;
extern char uvar_names[MAXODE][12];
float **my_adj;
float **my_orb;           /* the orbit on the grid of the adjoint: time, the variables, and the auxiliary quantities */
float **my_I,**my_R,**my_B;   /* the isostable response I, the Floquet eigenfunction rho, and the correction B, in the same layout */
int ISO_HERE=0,B_HERE=0;
static char iso_note[160];
static PHASE ph_keep;        /* the double-precision results of the last New adj (the orbit, Z, I, rho, B) */
static int PH_KEEP=0;
float **my_h2;
int H2_HERE=0,h2_len=0;
void new_h2_fun();
void h2_back();
int adj_len;
int adj_ncol=0;           /* number of columns allocated for my_adj and my_orb */
int ADJ_GRID=800;         /* grid steps for the adjoint (multiple of 200) */
float **my_h;
extern double last_ic[MAXODE];
void orbit_back();
void iso_back(int which);
extern int METHOD,DelayFlag,NFlags,TORUS,itor[MAXODE];
void extra();



extern char *info_message;
struct {
  int here,col0,ncol,colskip;
  int row0,nrow,rowskip; 
  float **data;
  char firstcol[11];
} my_trans;
   
int TRANPOSE_HERE=0;
extern double NEWT_ERR;
double ADJ_EPS=1.e-8,ADJ_ERR=1.e-3;
int ADJ_MAXIT=20,ADJ_HERE=0,H_HERE=0,h_len,HODD_EV=0;
int AdjRange=0;
extern double DELTA_T,BOUND;
int *coup_fun[MAXODE];
char *coup_string[MAXODE];

extern int *my_ode[];
extern int NSYM,NSYM_START,NCON,NCON_START;
/* extern Window main_win; */
extern int DCURY;

void init_trans()
{
  my_trans.here=0;
  strcpy(my_trans.firstcol,uvar_names[0]);
  my_trans.ncol=2;
  my_trans.nrow=1;
  my_trans.rowskip=1;
  my_trans.colskip=1;
  my_trans.row0=1;
  my_trans.col0=2;
}

void dump_transpose_info(fp,f)
     FILE *fp;
     int f;
{
  char bob[256];
  if(f==READEM)
    fgets(bob,255,fp);
  else
    fprintf(fp,"# Transpose variables etc\n");
  io_string(my_trans.firstcol,11,fp,f);
  io_int(&my_trans.ncol,fp,f,"n columns");
  io_int(&my_trans.nrow,fp,f,"n rows");
  io_int(&my_trans.rowskip,fp,f,"row skip");
  io_int(&my_trans.colskip,fp,f,"col skip");
  io_int(&my_trans.row0,fp,f,"row 0");


}

int do_transpose()
{
 int i,status;
 static char *n[]={"*0Column 1","NCols","ColSkip","Row 1","NRows","RowSkip"};
 char values[6][MAX_LEN_SBOX];
 sprintf(values[0],"%s",my_trans.firstcol);
 sprintf(values[1],"%d",my_trans.ncol);
 sprintf(values[2],"%d",my_trans.colskip);
 sprintf(values[3],"%d",my_trans.row0);
 sprintf(values[4],"%d",my_trans.nrow);
 sprintf(values[5],"%d",my_trans.rowskip);
 if(my_trans.here){
   
   for(i=0;i<=my_trans.nrow;i++)free(my_trans.data[i]);
   free(my_trans.data);
   my_trans.here=0;
   data_back();
 }
 status=do_string_box(6,6,1,"Transpose Data",n,values,33);
 if(status!=0){
   find_variable(values[0],&i);
   if(i>-1)
     my_trans.col0=i+1;
   else
     {
       err_msg("No such columns");
       return 0;
     }
   strcpy(my_trans.firstcol,values[0]);
   i=atoi(values[4]);
   if(i>=NEQ)i=NEQ-1;
   my_trans.nrow=i;
   my_trans.ncol=atoi(values[1]);
   my_trans.colskip=atoi(values[2]);
   my_trans.row0=atoi(values[3]);
   my_trans.rowskip=atoi(values[5]);
   return (create_transpose());
 }
 return 0; 
     
 
}
 
int create_transpose()
{
  int i,j;
  int inrow,incol;
  my_trans.data=(float **)malloc(sizeof(float *)*(NEQ+1));
  for(i=0;i<=my_trans.nrow;i++)
    my_trans.data[i]=(float *)malloc(sizeof(float)*my_trans.ncol);
  for(i=my_trans.nrow+1;i<=NEQ;i++)my_trans.data[i]=storage[i];
  for(j=0;j<my_trans.ncol;j++)
    my_trans.data[0][j]=j+1;

  
  for(i=0;i<my_trans.ncol;i++){
    incol=my_trans.col0-1+i*my_trans.colskip;
    if(incol>NEQ)
      incol=NEQ;
    for(j=0;j<my_trans.nrow;j++){
      inrow=my_trans.row0+j*my_trans.rowskip;
      if(inrow>storind)
	inrow=storind;
      my_trans.data[j+1][i]=storage[incol][inrow];
    }
  }
  
  set_browser_data(my_trans.data,1);
  /*  my_browser.data=my_trans.data;
      my_browser.col0=1; */
   refresh_browser(my_trans.ncol);
   my_trans.here=1;
   return 1;
}
			 

void alloc_h_stuff()
{
  int i;
 for(i=0;i<NODE ;i++){
   coup_fun[i]=(int *)malloc(100*sizeof(int));
   coup_string[i]=(char *)malloc(80);
   strcpy(coup_string[i],"0");
 }
}
 

void data_back()
{
 FOUR_HERE=0;
 set_browser_data(storage,1); 
 /*  my_browser.data=storage;
     my_browser.col0=1; */
 refresh_browser(storind); 
}

void adj_back()
{
 if(ADJ_HERE){
   set_browser_data(my_adj,1);
   /*   my_browser.data=my_adj;

	my_browser.col0=1; */
   refresh_browser(adj_len);
 }
}

void h_back()
{
 if(H_HERE){
   set_browser_data(my_h,1);
   /*
   my_browser.data=my_h;

   my_browser.col0=1; */
   refresh_browser(h_len);
 }
}
/*  Here is how to do the range over adjoints and h functions
    unfortunately, h functions are always computed even if you dont want them 
    they will just be zeros

    Step 1. Compute a singel orbit, adjoint, and H function (to load the 
    program with the correct right-hand sides for H function. Or just load in
    set file where it was done
    Step 2.  Set transient to some reasonable number to assure convergence 
    onto the limit cycle as you change parameters and total to be at least
    2 periods beyond the transient
    Step 3. Se up Poincare map - period - stop on section. This lets you
    get the period
    Step 4. In numerics averaging - click on adjrange
    Step 5. Initconds range over the parameter. It should find the periodic
    orbit, adjoint, and H function and save. Files are of the form
    orbit.parname_parvalue.dat etc
*/
void make_adj_com(int com)
{
static char key[]="nmaohprifb   2";
 switch(key[com]){
 case 'n': 
   new_adjoint();
   break;
 case 'm':
   new_h_fun(0);
   break;
 case 'a':
   adj_back();
   break;
 case 'o':
   orbit_back();
   break;
 case 'h':
   h_back();
   break;
 case 'p':
   adjoint_parameters();
   break;
 case 'r':
   AdjRange=1;
   break;
 case '2':
   new_h2_fun();
   break;
 case 'i':
   iso_back(0);
   break;
 case 'f':
   iso_back(1);
   break;
 case 'b':
   iso_back(2);
   break;
 
 }
}


void adjoint_parameters()
{
  new_int("Grid steps for the adjoint (a multiple of 200, 200-6400) :",&ADJ_GRID);
  if(ADJ_GRID<200)ADJ_GRID=200;
  if(ADJ_GRID>6400)ADJ_GRID=6400;
}

void new_h_fun(silent)
     int silent;
{

 int i,n=2;
 if(!ADJ_HERE){
   err_msg("Must compute adjoint first!");
   return;
 }
 if(H_HERE){
   free(my_h[0]);
   free(my_h[1]);
   if(HODD_EV){
     free(my_h[2]);
     free(my_h[3]);
   }
   free(my_h);
   H_HERE=0;
   HODD_EV=0;
 }
   if(NEQ>2){
     HODD_EV=1;
     n=4;
   }
   h_len=adj_len;
   data_back();
   my_h=(float **)malloc(sizeof(float*)*(NEQ+1));
   for(i=0;i<n;i++)my_h[i]=(float *)malloc(sizeof(float)*h_len);
   for(i=n;i<=NEQ;i++)my_h[i]=my_orb[i];   /* the auxiliary columns are those of the orbit */
   if(make_h(my_orb,my_adj,my_h,h_len,0.0,NODE,silent )){
     double hmn=1.e300,hmx=-1.e300,hm=0.0;
     int jj;
     H_HERE=1;
     for(jj=0;jj<h_len-1;jj++){   /* the last point repeats the first */
       double v=my_h[1][jj];
       hm+=v; if(v<hmn)hmn=v; if(v>hmx)hmx=v;
     }
     plintf("H computed over %d points: H(0)=%.10g, mean %.10g, max %.10g, min %.10g\n",h_len-1,(double)my_h[1][0],hm/(double)(h_len-1),hmx,hmn);
     h_back();
   }
 ping();
  
}

void dump_h_stuff(fp,f)
     FILE *fp;
     int f;
{
  char bob[256];
  int i;
  if(f==READEM)
    fgets(bob,255,fp);
  else
    fprintf(fp,"# Coupling stuff for H funs\n");
 for(i=0;i<NODE ;i++)
   io_string(coup_string[i],79,fp,f);

}


int make_h(orb,adj,h,nt,dt,node,silent)
float **orb,**adj,**h;
double dt;
int node,nt,silent;
{

 int i,j,rval=0;
 double sum;
 double z;
 int n0=node+1+FIX_VAR,k2,k;
 char name[30];
 if(silent==0){
   for(i=0;i<NODE ;i++){
     sprintf(name,"Coupling for %s eqn:",uvar_names[i]);
     new_string(name,coup_string[i]);
     if(add_expr(coup_string[i],coup_fun[i],&j)){
       err_msg("Illegal formula");
       goto bye;
     }
   }
 }
 /*  formulae are fine .. lets do it ... */
 for(j=0;j<nt;j++)  /* j is phi variable  */
   {
     sum=0.0;
     
     for(k=0;k<nt-1;k++){      /* the last point repeats the first: nt-1 distinct points */
         k2=k+j;
       if(k2>=nt-1)k2=k2-nt+1;
       for(i=0;i<node;i++){
	 set_ivar(i+1,(double)orb[i+1][k]);
         set_ivar(i+n0+1,(double)orb[i+1][k2]);
       }
       z=0.0;
       update_based_on_current(); 

       for(i=0;i<node;i++){
	
	 z=evaluate(coup_fun[i]);
	
	 sum=sum+z*(double)adj[i+1][k];
       }
	
     }
     my_h[0][j]=orb[0][j];
     my_h[1][j]=(float)(sum/(double)(nt-1));
   }
   if(HODD_EV){
     for(k=0;k<nt;k++){
       k2=nt-k-1;
       my_h[2][k]=.5*(my_h[1][k]-my_h[1][k2]);
       my_h[3][k]=.5*(my_h[1][k]+my_h[1][k2]);
     }
   }
   rval=1;
      
 bye:
  NSYM=NSYM_START;
  NCON=NCON_START;
  return(rval);

}
 

 
  
 


/* free the adjoint, the orbit and the H function computed from them */
static void adj_free_all()
{
  int i;
  if(H_HERE){
    free(my_h[0]);
    free(my_h[1]);
    if(HODD_EV){
      free(my_h[2]);
      free(my_h[3]);
    }
    free(my_h);
    H_HERE=0;
    HODD_EV=0;
  }
  if(H2_HERE){for(i=0;i<=NEQ;i++)free(my_h2[i]);free(my_h2);H2_HERE=0;}
  if(PH_KEEP){phase_free(&ph_keep);PH_KEEP=0;}
  if(ADJ_HERE){
    for(i=0;i<adj_ncol;i++){
      free(my_adj[i]);
      free(my_orb[i]);
      if(ISO_HERE){free(my_I[i]);free(my_R[i]);}
      if(B_HERE)free(my_B[i]);
    }
    free(my_adj);
    free(my_orb);
    if(ISO_HERE){free(my_I);free(my_R);}
    if(B_HERE)free(my_B);
    ADJ_HERE=0;ISO_HERE=0;B_HERE=0;
  }
}

/* a data set with columns: time, then n columns, then zeros up to NEQ */
static float **iso_columns(const double *v,int n,int P,const PHASE *ph)
{
  float **d=(float **)malloc(sizeof(float *)*(NEQ+1));
  int i,k;
  if(d==NULL)return NULL;
  for(i=0;i<=NEQ;i++){
    d[i]=(float *)calloc((size_t)P,sizeof(float));
    if(d[i]==NULL){for(k=0;k<i;k++)free(d[k]);free(d);return NULL;}
  }
  for(k=0;k<P;k++){
    d[0][k]=(float)(ph->th[k]*ph->T);
    for(i=0;i<n;i++)d[i+1][k]=(float)v[(size_t)k*n+i];
  }
  return d;
}

void iso_back(int which)    /* 0: I, 1: rho, 2: B */
{
  if(!ADJ_HERE){err_msg("Must compute adjoint first!");return;}
  if(!ISO_HERE){
    char msg[260];
    sprintf(msg,"No isostable quantities: %s",iso_note);
    err_msg(msg);
    return;
  }
  if(which==2&&!B_HERE){err_msg(iso_note);return;}
  set_browser_data(which==0?my_I:(which==1?my_R:my_B),1);
  refresh_browser(adj_len);
}

/* the orbit that goes with the adjoint (or, with no adjoint yet, the stored solution) in the browser */
void orbit_back()
{
  if(ADJ_HERE){
    set_browser_data(my_orb,1);
    refresh_browser(adj_len);
  }
  else data_back();
}

/*  The adjoint Z of the stable limit cycle that starts at the initial conditions and has period Total (numerical
    parameters). That is just what (I)nitialconds (P)eriodic sets up. The orbit is computed again, in double precision,
    on a grid of ADJ_GRID steps (phase_red.c), so the stored solution is not used. Z is loaded into the browser;
    (O)rbit gives the orbit on the same grid.  */
void new_adjoint()
{
  int i,k,P,n=NODE;
  PHASE ph;
  double x0[MAXODE],yy[MAXODE],t;
  char msg[200];
  if(METHOD==0){err_msg("The adjoint is for differential equations, not maps");return;}
  if(METHOD==6||DelayFlag||NMarkov>0||NFlags>0){
    err_msg("The adjoint is not available with delays, Volterra equations, Markov variables or flags");
    return;
  }
  if(TORUS)for(i=0;i<n;i++)if(itor[i]){err_msg("The adjoint is not available with periodic (torus) variables");return;}
  if(n>PHASE_MAXN){
    sprintf(msg,"The adjoint is computed for up to %d equations",PHASE_MAXN);
    err_msg(msg);
    return;
  }
  if(!(TEND>0.0)){err_msg("Total must be the period of the cycle (Initialconds Periodic sets it)");return;}
  for(i=0;i<n;i++)x0[i]=last_ic[i];
  if(!phase_adjoint(n,x0,TEND,ADJ_GRID,&ph)){
    err_msg("No adjoint: start on a stable limit cycle, with Total equal to its period (Initialconds Periodic does this)");
    return;
  }
  if(ADJ_HERE){
    data_back();
    adj_free_all();
  }
  P=ph.P;
  adj_len=P;
  adj_ncol=NEQ+1;
  my_adj=(float **)malloc(sizeof(float *)*(NEQ+1));
  my_orb=(float **)malloc(sizeof(float *)*(NEQ+1));
  for(i=0;i<=NEQ;i++){
    my_adj[i]=(float *)calloc((size_t)P,sizeof(float));
    my_orb[i]=(float *)calloc((size_t)P,sizeof(float));
    if(my_adj[i]==NULL||my_orb[i]==NULL){
      err_msg("Insufficient storage");
      phase_free(&ph);
      adj_ncol=i+1;ADJ_HERE=1;      /* what was allocated is freed by adj_free_all */
      for(k=i+1;k<=NEQ;k++){my_adj[k]=NULL;my_orb[k]=NULL;}
      adj_ncol=NEQ+1;
      adj_free_all();
      return;
    }
  }
  for(k=0;k<P;k++){
    t=ph.th[k]*ph.T;
    my_adj[0][k]=(float)t;
    my_orb[0][k]=(float)t;
    for(i=0;i<n;i++){
      yy[i]=ph.U[(size_t)k*n+i];
      my_orb[i+1][k]=(float)yy[i];
      my_adj[i+1][k]=(float)ph.Z[(size_t)k*n+i];
    }
    if(NEQ>n){            /* the auxiliary quantities along the orbit */
      for(i=n;i<NEQ;i++)yy[i]=0.0;
      extra(yy,t,n,NEQ);
      for(i=n;i<NEQ;i++)my_orb[i+1][k]=(float)yy[i];
    }
  }
  ADJ_HERE=1;
  ISO_HERE=0;B_HERE=0;
  iso_note[0]=0;
  if(ph.iso){
    my_I=iso_columns(ph.I,n,P,&ph);
    my_R=iso_columns(ph.R,n,P,&ph);
    ISO_HERE=(my_I!=NULL&&my_R!=NULL);
    if(!ISO_HERE){free(my_I);free(my_R);}
    if(ph.hasB){
      my_B=iso_columns(ph.B,n,P,&ph);
      B_HERE=(my_B!=NULL);
    }
    else sprintf(iso_note,"%s",phase_b_note(ph.bNote));
  }
  else sprintf(iso_note,"%s",phase_iso_note(ph.isoSkip));
  plintf("Adjoint of the cycle: period %.12g (grid of %d steps), Z.f = 1 to %.1e, Z periodic to %.1e\n",
	 ph.T,ph.Ng,ph.errNorm,ph.errPer);
  plintf("  Z at phase 0:");
  for(i=0;i<n&&i<8;i++)plintf(" %.8g",ph.Z[i]);
  plintf("\n");
  if(ph.iso){
    plintf("  isostable: kappa = %.10g, nu = exp(kappa T) = %.10g, next multiplier / nu = %.3g; I.rho = 1 to %.1e, Z.rho = 0 to %.1e\n",
	   ph.kappa,ph.nu,ph.nuRatio,ph.errDot,ph.errZR);
    if(ph.hasB)plintf("  B found, periodic to %.1e\n",ph.errPerB);
    else plintf("  B not computed: %s\n",iso_note);
  }
  else plintf("  no isostable quantities: %s\n",iso_note);
  if(ph.closure1>=0.0)plintf("  (the grid orbit was closed by Newton steps: closure error %.2e -> %.2e)\n",ph.closure0,ph.closure1);
  if(ph.multiple>1)plintf("  Total is %d times the period: the adjoint is for the period %.12g\n",ph.multiple,ph.T);
  if(fabs(ph.T-TEND)>1.e-3*TEND){
    sprintf(msg,"Total=%g is not the period: the cycle has period %.8g",TEND,ph.T);
    err_msg(msg);
  }
  else if(ph.errNorm>1.e-4||ph.errPer>1.e-4){
    err_msg("The adjoint may be inaccurate: is the cycle stable? Try more grid steps in (P)arameters");
  }
  ph_keep=ph;PH_KEEP=1;           /* kept for H2 */
  adj_back();
  ping();
}

/* this computes the periodic orbit and stores it in 
   the usual place  given initial data and period */
  


void test_test()
{
  double x[2];
  x[0]=.35249;
  x[1]=.2536;
  compute_one_orbit(x,14.6);
 
}

void compute_one_orbit(double *ic,double per)
{
  double oldtotal=TEND;
  TEND=per;
  /*   plintf(" %g %g \n",ic[0],ic[1]); */
  run_from_x(ic);
  new_adjoint();
  TEND=oldtotal;
}
  



/*    ADJOINT ROUTINE
 *
 *
      This assumes that you have already computed the periodic orbit
	and have stored in in an array **orbit
      including time in the first column
     
     
      The righthand sides of the equations are
	rhs(t,y,yp,n)
      and the coupling function for ``H'' functions is
	couple(y,yhat,f,n)

	where yhat is presynaptic and y is postynaptic
     variable.  f returns the coupling vector.

    adjoint is the same size as orbit and when returned has 
    t in the first column.
  */

int adjoint(orbit,adjnt,nt,dt,eps,minerr,maxit,node)
 float **orbit,**adjnt;
 double dt,eps,minerr;
 int nt,node,maxit;
 {
  double **jac,*yold,ytemp,*fold,*fdev;
  double *yprime,*work;
  double t,prod,del;
  int i,j,k,l,k2,rval=0;
  int n2=node*node;
  double error;
   
   work = (double *)malloc((n2+4*node)*sizeof(double));
   yprime = (double *)malloc(node*sizeof(double));
   yold=(double *)malloc(node*sizeof(double));
   fold=(double *)malloc(node*sizeof(double));
   fdev=(double *)malloc(node*sizeof(double));
  jac = (double **)malloc(n2*sizeof(double *));
  
  for(i=0;i<n2;i++)
  {
   jac[i]=(double *)malloc(nt*sizeof(double));
       if(jac[i]==NULL){
       err_msg("Insufficient storage");
	return(0);
	}
   }
  
  /*  Now we compute the
	transpose time reversed jacobian  --  this is complex !! */
  for(k=0;k<nt;k++){
	l=nt-1-k;  /* reverse the limit cycle  */
	for(i=0;i<node;i++)yold[i]=(double)orbit[i+1][l];
        rhs(0.0,yold,fold,node);
	for(j=0;j<node;j++){
		ytemp=yold[j];
		del=eps*fabs(ytemp);
		if(del<eps)del=eps;
		
		yold[j]+=del;
		rhs(0.0,yold,fdev,node);
		yold[j]=ytemp;
		for(i=0;i<node;i++)
			jac[i+node*j][k]=(fdev[i]-fold[i])/del;
		
	      }
 
      }
 
  /* now we iterate to get a good adjoint using implicit Euler's method */
  ytemp=0.0;
 for(i=0;i<node;i++){
	yold[i]=1.+.01*(ndrand48()-.5); /* random initial data */
	
	ytemp+=fabs(yold[i]);
	}
  for(i=0;i<node;i++){
    yold[i]=yold[i]/ytemp;
    fdev[i]=yold[i];
  }
	
  plintf("%f %f \n",yold[0],yold[1]);
      
 for(l=0;l<maxit;l++){
	for(k=0;k<nt-1;k++){
		k2=k+1;
		if(k2>=nt)k2=k2-nt;
		if(step_eul(jac,k,k2,yold,work,node,dt)==0){
		  rval=0;
		  goto bye;
		} 
		/* rk_interp(jac,k,k2,yold,work,node,dt,5); */
	      }
	ytemp=0.0;
	error=0.0;
	
         for(i=0;i<node;i++){
	  if(fabs(yold[i])>BOUND){
	    
	     rval=0;
	  err_msg("Out of bounds");
	  goto bye;
	  }
	error+=fabs(yold[i]-fdev[i]);
	ytemp+=fabs(yold[i]);
	
	}
	
	 
	
        for(i=0;i<node;i++){ yold[i]=yold[i]/ytemp;
			     fdev[i]=yold[i];
			   }
	printf("%f %f \n",yold[0],yold[1]);
        plintf("err=%f \n",error);
	if(error<minerr)break; /*  exit if error small   */
 }
 /*  onelast time to compute the adjoint  */
 prod=0.0;   /* for normalization   */
 t=0.0;
 for(k=0;k<nt;k++){
        l=nt-k-1;
	t+=dt; 
        for(i=0;i<node;i++)fdev[i]=(double)orbit[i+1][l];
	rhs(0.0,fdev,yprime,node);
	for(j=0;j<node;j++){
	adjnt[j+1][l]=(float)yold[j];
	prod+=yold[j]*yprime[j]*dt;
      }
	k2=k+1;
	if(k2>=nt)k2-=nt;
	 if(step_eul(jac,k,k2,yold,work,node,dt)==0){
	  rval=0;
	  goto bye;
	    }	 
	/* rk_interp(jac,k,k2,yold,work,node,dt,5); */
	

      }

        	 
	prod=prod/t;
  plintf(" Multiplying the adjoint by 1/%g to normalize\n",prod);
  for(k=0;k<nt;k++){
     for(j=0;j<node;j++)adjnt[j+1][k]=adjnt[j+1][k]/(float)prod;
     adjnt[0][k]=orbit[0][k];
   }
  rval=1;

 bye: 
   free(work);  
   free(yprime);
   free(yold);
   free(fold);
   free(fdev);
   for(i=0;i<n2;i++)
   free(jac[i]);
   free(jac); 
   return(rval);
 }
 
 

void eval_rhs(jac,k1,k2,t,y,yp,node)
     double t,**jac,*y,*yp;
     int node,k1,k2;
{
  int i;
  int j;
  for(j=0;j<node;j++){
    yp[j]=0.0;
    for(i=0;i<node;i++)
 yp[j]=yp[j]+(jac[i+j*node][k1]*(1.0-t)+jac[i+j*node][k2]*t)*y[i];
  }
}

int rk_interp(jac,k1,k2,y,work,neq,del,nstep)
double *y,del,*work,**jac;
int neq,k1,k2,nstep;
{
 int i,j;
 double *yval[3],dt=del/nstep;
 double t=0.0,t1,t2;
 yval[0]=work;
 yval[1]=work+neq;
 yval[2]=work+neq+neq;
 for(j=0;j<nstep;j++)
 {
  eval_rhs(jac,k1,k2,t/del,y,yval[1],neq);
  for(i=0;i<neq;i++)
  {
   yval[0][i]=y[i]+dt*yval[1][i]/6.00;
   yval[2][i]=y[i]+dt*yval[1][i]*0.5;
  }
  t1=t+.5*dt;
  eval_rhs(jac,k1,k2,t1/del,yval[2],yval[1],neq);
  for(i=0;i<neq;i++)
  {
   yval[0][i]=yval[0][i]+dt*yval[1][i]/3.00;
   yval[2][i]=y[i]+.5*dt*yval[1][i];
  }
  eval_rhs(jac,k1,k2,t1/del,yval[2],yval[1],neq);
  for(i=0;i<neq;i++)
  {
   yval[0][i]=yval[0][i]+dt*yval[1][i]/3.000;
   yval[2][i]=y[i]+dt*yval[1][i];
  }
  t2=t+dt;
  eval_rhs(jac,k1,k2,t2/del,yval[2],yval[1],neq);
  for(i=0;i<neq;i++)y[i]=yval[0][i]+dt*yval[1][i]/6.00;
  t=t2;
 }
return(1);  
}
 
	 
int step_eul(jac,k,k2,yold,work,node,dt)
double *work,*yold,**jac,dt;
int k,k2,node;
{

int j,i,n2=node*node,info;
int ipvt[MAXODE];
double *mat,*fold;
fold=work;
mat=work+node;

  for(j=0;j<node;j++){
    fold[j]=0.0;
    for(i=0;i<node;i++)
      fold[j]=fold[j]+jac[i+j*node][k]*yold[i];
  }
  for(j=0;j<node;j++)yold[j]=yold[j]+.5*dt*fold[j];
  for(i=0;i<n2;i++)mat[i]=-jac[i][k2]*dt*.5;
  for(i=0;i<node;i++)mat[i+i*node]=1.+mat[i+i*node];
  sgefa(mat,node,node,ipvt,&info);
if(info!=-1){
  
  err_msg("Univertible Jacobian");
  return(0);
}
sgesl(mat,node,node,ipvt,yold,0);
return(1);
}
		


		
void norm_vec(v,mu,n) /* returns the length of the vector and the unit vector */
     double *v,*mu;
     int n;
{
  int i;
  double sum=0.0;
  for(i=0;i<n;i++)
    sum+=(v[i]*v[i]);
  sum=sqrt(sum);
  if(sum>0)
    for(i=0;i<n;i++)
      v[i]=v[i]/sum;
  *mu=sum;
  return;
}


 






/*  H2: the second-order interaction function of two identical cells coupled through the formulas of (M)ake H
    (phase_red.c, after Wilson and Ermentrout, Phys. Rev. Lett. 123, 164101 (2019)). Needs I, rho and B.  */
static int h2_active[MAXODE];

static void xpp_couple(const double *own,const double *partner,double *out,void *ctx)
{
  int i,n0=NODE+1+FIX_VAR;
  (void)ctx;
  for(i=0;i<NODE;i++){
    set_ivar(i+1,own[i]);
    set_ivar(i+n0+1,partner[i]);
  }
  update_based_on_current();
  for(i=0;i<NODE;i++)out[i]=h2_active[i]?evaluate(coup_fun[i]):0.0;
}

void new_h2_fun()
{
  int i,j,nc=0,idx[MAXODE],Mp;
  PHASEH h;
  char msg[260];
  if(!ADJ_HERE||!PH_KEEP){err_msg("Must compute adjoint first!");return;}
  if(!ph_keep.iso){sprintf(msg,"H2 needs I, rho and B: %s",iso_note);err_msg(msg);return;}
  if(!ph_keep.hasB){sprintf(msg,"H2 needs B: %s",iso_note);err_msg(msg);return;}
  if((ph_keep.P-1)%PHASE_HM!=0){err_msg("H2 needs a grid with a multiple of 200 steps");return;}
  if(!H_HERE){                     /* the coupling formulas are asked for in Make H */
    new_h_fun(0);
    if(!H_HERE)return;
  }
  for(i=0;i<NODE;i++){             /* make_h forgets the compiled formulas when it ends */
    if(add_expr(coup_string[i],coup_fun[i],&j)){
      err_msg("Illegal formula");
      NSYM=NSYM_START;NCON=NCON_START;
      return;
    }
    h2_active[i]=(strcmp(coup_string[i],"0")!=0);
    if(h2_active[i])idx[nc++]=i;
  }
  if(nc==0){err_msg("No coupling is defined: give the coupling in (M)ake H");NSYM=NSYM_START;NCON=NCON_START;return;}
  if(!phase_hfun(&ph_keep,idx,nc,xpp_couple,NULL,&h)){
    err_msg("Could not compute H2");
    NSYM=NSYM_START;NCON=NCON_START;
    return;
  }
  NSYM=NSYM_START;NCON=NCON_START;
  Mp=h.Mp;
  if(H2_HERE){for(i=0;i<=NEQ;i++)free(my_h2[i]);free(my_h2);H2_HERE=0;}
  h2_len=Mp+1;
  my_h2=(float **)malloc(sizeof(float *)*(NEQ+1));
  for(i=0;i<=NEQ;i++)my_h2[i]=(float *)calloc((size_t)h2_len,sizeof(float));
  for(j=0;j<=Mp;j++){
    int jm=(Mp-j)%Mp;
    my_h2[0][j]=(float)(h.T*j/Mp);
    my_h2[1][j]=(float)h.H2[j];
    if(NEQ>=2)my_h2[2][j]=(float)(0.5*(h.H2[j]-h.H2[jm]));      /* odd part */
    if(NEQ>=3)my_h2[3][j]=(float)(0.5*(h.H2[j]+h.H2[jm]));      /* even part */
  }
  H2_HERE=1;
  {
    double mx=0.0;
    for(j=0;j<=Mp;j++)if(fabs(h.H2[j])>mx)mx=fabs(h.H2[j]);
    plintf("H2 computed (two coupled cells, second order): H2(0)=%.10g, max |H2|=%.8g\n",h.H2[0],mx);
    plintf("  g2+g3 = -2 odd(H2); gamma = (g2'(0)+g3'(0))/4 = %.8g (from A: %.8g, from Bf: %.8g)\n",h.gamma,h.gammaA,h.gammaB);
    if(h.notes&1)plintf("  note: the slowest multiplier is close to 1, where the expansion in epsilon breaks down\n");
    if(h.notes&2)plintf("  note: the next multiplier is not much smaller (ratio %.3g), so one isostable is not enough\n",ph_keep.nuRatio);
  }
  phase_h_free(&h);
  h2_back();
  ping();
}

void h2_back()
{
  if(H2_HERE){
    set_browser_data(my_h2,1);
    refresh_browser(h2_len);
  }
}
