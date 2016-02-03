#include "load_eqn.h"
#include "parserslow.h"

#include "read_dir.h"

#include "main.h"
#include "ggets.h"
#include "dae_fun.h"
#include "derived.h"
#include "extra.h"
#include "init_conds.h"
#include "browse.h"
#include "txtread.h"
#include "numerics.h"
#include "integrate.h"
#include "adj2.h"
#include "arrayplot.h"
#include "lunch-new.h"
#include "graphics.h"

/*#include "macdirent.h"
*/

#include <dirent.h>
#include "userbut.h"
#include "volterra2.h"
#include "storage.h"
#include "tabular.h"

#include <stdlib.h> 
#include <string.h>
#include <stdio.h>
#include "xpplim.h"

#define PARAM 1
#define IC 2


#define DFNORMAL 1
#define MAXOPT 1000
#define READEM 1

extern OptionsSet notAlreadySet;

typedef struct {
  int nbins,nbins2,type,col,col2,fftc;
  double xlo,xhi;
  double ylo,yhi;
  char cond[80];
} HIST_INFO;

extern HIST_INFO hist_inf;
extern int spec_col,spec_wid,spec_win,spec_col2,post_process;



int nsrand48(int seed);

char *interopt[MAXOPT];
int Nopts=0;
int RunImmediately=0;
extern char dll_lib[256];
extern char dll_fun[256];
extern int dll_flag;

extern char UserBlack[8];
extern char UserWhite[8];
extern char UserMainWinColor[8];
extern char UserDrawWinColor[8];
/*extern char UserBGBitmap[100];*/
extern char UserBGBitmap[XPP_MAX_NAME];

extern int UserGradients;
extern int UserMinWidth;
extern int UserMinHeight;
extern int XPPVERBOSE;
extern FILE *logfile;

extern int OVERRIDE_QUIET;
extern int OVERRIDE_LOGFILE;

extern int SLIDER1;
extern int SLIDER2;
extern int SLIDER3;
extern char SLIDER1VAR[20];
extern char SLIDER2VAR[20];
extern char SLIDER3VAR[20];
extern double SLIDER1LO;
extern double SLIDER2LO;
extern double SLIDER3LO;
extern double SLIDER1HI;
extern double SLIDER2HI;
extern double SLIDER3HI;
extern double SLIDER1INIT;
extern double SLIDER2INIT;
extern double SLIDER3INIT;

extern int NCBatch,DFBatch;
extern int DF_GRID;

typedef struct {
  char *name;
  char *does;
  unsigned int use;
} INTERN_SET;

extern int XNullColor,YNullColor,StableManifoldColor,UnstableManifoldColor;
int IX_PLT[10],IY_PLT[10],IZ_PLT[10],NPltV;
int MultiWin=0;
extern int SimulPlotFlag;
double X_LO[10],Y_LO[10],X_HI[10],Y_HI[10];
int START_LINE_TYPE=1;
INTERN_SET intern_set[MAX_INTERN_SET];
int Nintern_set=0;

extern int STOCH_FLAG;
extern char uvar_names[MAXODE][12]; 
extern struct {
         char item[30],item2[30];
	 int steps,steps2,reset,oldic,index,index2,cycle,type,type2,movie;
	 double plow,phigh,plow2,phigh2;
         int rtype;
       } range;

extern int custom_color;
extern int del_stab_flag;
extern int MaxPoints;
extern double THETA0,PHI0;
extern int tfBell;
extern int DoTutorial;
/*void set_option(char *s1,char *s2);
*/

double atof();
char *get_first();
char *get_next();
/*   this file has all of the phaseplane parameters defined   
     and created.  All other files should use external stuff
    to use them. (Except eqn forming stuff)
 */

extern char batchout[256];
extern int batch_range; 
 double last_ic[MAXODE];
extern char PlotFormat[100];
 extern char big_font_name[100],small_font_name[100];
 extern int PaperWhite;

extern int PSColorFlag,PS_FONTSIZE,PS_Color;
extern char PS_FONT[100];
extern double PS_LW;

extern int SEc,UEc,SPc,UPc;
int (*solver)();

 int rung_kut();
 char delay_string[MAXODE][80];
 int itor[MAXODE];
 /*char this_file[100];
 */
 char this_file[XPP_MAX_NAME];
 char this_internset[XPP_MAX_NAME];
 float oldhp_x,oldhp_y,my_pl_wid,my_pl_ht;
 int mov_ind;
 int  storind,STORFLAG,INFLAG,MAXSTOR;
 double x_3d[2],y_3d[2],z_3d[2];
 int IXPLT,IYPLT,IZPLT;
 int AXES,TIMPLOT,PLOT_3D;
 double MY_XLO,MY_YLO,MY_XHI,MY_YHI;
 double TOR_PERIOD=6.2831853071795864770;
 int TORUS=0;
 int NEQ;
 char options[100];  

/*   Numerical stuff ....   */

 double DELTA_T,TEND,T0,TRANS,
	NULL_ERR,EVEC_ERR,NEWT_ERR;
 double BOUND,DELAY,TOLER,ATOLER,HMIN,HMAX;
 double BVP_EPS,BVP_TOL;

 double POIPLN;

 extern int RandSeed;
 int MaxEulIter;
double EulTol;
extern int cv_bandflag,cv_bandupper,cv_bandlower;
 int NMESH,NJMP,METHOD,color_flag,NC_ITER;
 int EVEC_ITER;
 int BVP_MAXIT,BVP_FLAG;

 int POIMAP,POIVAR,POISGN,SOS;
   int FFT,NULL_HERE,POIEXT;
  int HIST,HVAR,hist_ind,FOREVER;

 /*  control of range stuff  */

 int PAUSER,ENDSING,SHOOT,PAR_FOL;


/*  custon color stuff  */

extern char ColorVia[15];
extern double ColorViaLo,ColorViaHi;
extern int ColorizeFlag;


/* AUTO STUFF  */
extern int auto_ntst,auto_nmx,auto_npr,auto_ncol;
extern double auto_ds,  auto_dsmax,  auto_dsmin;
extern double auto_rl0,auto_rl1,auto_a0,auto_a1;
extern double auto_epss,auto_epsl,auto_epsu;
extern int auto_var;
extern double auto_xmin,auto_xmax,auto_ymin,auto_ymax;

 extern int PltFmtFlag;

 int xorfix,silent,got_file;

/*Logical negate OR on options set. Result overwrites the first OptionsSet
in the argument list.*/
void notBothOptions(OptionsSet nasA,OptionsSet nasB)
{
   nasA.BIG_FONT_NAME = (nasA.BIG_FONT_NAME & nasB.BIG_FONT_NAME);
   nasA.SMALL_FONT_NAME = (nasA.SMALL_FONT_NAME & nasB.SMALL_FONT_NAME);
   nasA.BACKGROUND = (nasA.BACKGROUND & nasB.BACKGROUND);
   nasA.IXPLT = (nasA.IXPLT & nasB.IXPLT);
   nasA.IYPLT = (nasA.IYPLT & nasB.IYPLT);
   nasA.IZPLT = (nasA.IZPLT & nasB.IZPLT);
   nasA.AXES = (nasA.AXES & nasB.AXES);
   nasA.NMESH = (nasA.NMESH & nasB.NMESH);
   nasA.METHOD = (nasA.METHOD & nasB.METHOD);
   nasA.TIMEPLOT = (nasA.TIMEPLOT & nasB.TIMEPLOT);
   nasA.MAXSTOR = (nasA.MAXSTOR & nasB.MAXSTOR);
   nasA.TEND = (nasA.TEND & nasB.TEND);
   nasA.DT = (nasA.DT & nasB.DT);
   nasA.T0 = (nasA.T0 & nasB.T0);
   nasA.TRANS = (nasA.TRANS & nasB.TRANS);
   nasA.BOUND = (nasA.BOUND & nasB.BOUND);
   nasA.TOLER = (nasA.TOLER & nasB.TOLER);
   nasA.DELAY = (nasA.DELAY & nasB.DELAY);
   nasA.XLO = (nasA.XLO & nasB.XLO);
   nasA.XHI = (nasA.XHI & nasB.XHI);
   nasA.YLO = (nasA.YLO & nasB.YLO);
   nasA.YHI = (nasA.YHI & nasB.YHI);  
   nasA.UserBlack = (nasA.UserBlack & nasB.UserBlack);
   nasA.UserWhite = (nasA.UserWhite & nasB.UserWhite);
   nasA.UserMainWinColor = (nasA.UserMainWinColor & nasB.UserMainWinColor);
   nasA.UserDrawWinColor = (nasA.UserDrawWinColor & nasB.UserDrawWinColor);
   nasA.UserGradients = (nasA.UserGradients & nasB.UserGradients);
   nasA.UserBGBitmap = (nasA.UserBGBitmap & nasB.UserBGBitmap);
   nasA.UserMinWidth = (nasA.UserMinWidth & nasB.UserMinWidth);
   nasA.UserMinHeight = (nasA.UserMinHeight & nasB.UserMinHeight);
   nasA.YNullColor = (nasA.YNullColor & nasB.YNullColor);
   nasA.XNullColor = (nasA.XNullColor & nasB.XNullColor);
   nasA.StableManifoldColor = (nasA.StableManifoldColor & nasB.StableManifoldColor);
   nasA.UnstableManifoldColor = (nasA.UnstableManifoldColor & nasB.UnstableManifoldColor);
   nasA.START_LINE_TYPE = (nasA.START_LINE_TYPE & nasB.START_LINE_TYPE);
   nasA.RandSeed = (nasA.RandSeed & nasB.SMALL_FONT_NAME);
   nasA.PaperWhite = (nasA.PaperWhite & nasB.PaperWhite);
   nasA.COLORMAP = (nasA.COLORMAP & nasB.COLORMAP);
   nasA.NPLOT = (nasA.NPLOT & nasB.NPLOT);
   nasA.DLL_LIB = (nasA.DLL_LIB & nasB.DLL_LIB);
   nasA.DLL_FUN = (nasA.DLL_FUN & nasB.DLL_FUN);
   nasA.XP = (nasA.XP & nasB.XP);
   nasA.YP = (nasA.YP & nasB.YP);
   nasA.ZP = (nasA.ZP & nasB.ZP);
   nasA.NOUT = (nasA.NOUT & nasB.NOUT);
   nasA.VMAXPTS = (nasA.VMAXPTS & nasB.VMAXPTS);
   nasA.TOR_PER = (nasA.TOR_PER & nasB.TOR_PER);
   nasA.JAC_EPS = (nasA.JAC_EPS & nasB.JAC_EPS);
   nasA.NEWT_TOL = (nasA.NEWT_TOL & nasB.NEWT_TOL);
   nasA.NEWT_ITER = (nasA.NEWT_ITER & nasB.NEWT_ITER);
   nasA.FOLD = (nasA.FOLD & nasB.FOLD);
   nasA.DTMIN = (nasA.DTMIN & nasB.DTMIN);
   nasA.DTMAX = (nasA.DTMAX & nasB.DTMAX);
   nasA.ATOL = (nasA.ATOL & nasB.ATOL);
   nasA.TOL = (nasA.TOL & nasB.TOL);
   nasA.BANDUP = (nasA.BANDUP & nasB.BANDUP);
   nasA.BANDLO = (nasA.BANDLO & nasB.BANDLO);
   nasA.PHI = (nasA.PHI & nasB.PHI);
   nasA.THETA = (nasA.THETA & nasB.THETA);
   nasA.XMIN = (nasA.XMIN & nasB.XMIN);
   nasA.XMAX = (nasA.XMAX & nasB.XMAX);
   nasA.YMIN = (nasA.YMIN & nasB.YMIN);
   nasA.YMAX = (nasA.YMAX & nasB.YMAX);
   nasA.ZMIN = (nasA.ZMIN & nasB.ZMIN);
   nasA.ZMAX = (nasA.ZMAX & nasB.ZMAX);
   nasA.POIVAR = (nasA.POIVAR & nasB.POIVAR);
   nasA.OUTPUT = (nasA.OUTPUT & nasB.OUTPUT);
   nasA.POISGN = (nasA.POISGN & nasB.POISGN);  
   nasA.POIEXT = (nasA.POIEXT & nasB.POIEXT);
   nasA.POISTOP = (nasA.POISTOP & nasB.POISTOP);
   nasA.STOCH = (nasA.STOCH & nasB.STOCH);
   nasA.POIPLN = (nasA.POIPLN & nasB.POIPLN);
   nasA.POIMAP = (nasA.POIMAP & nasB.POIMAP);
   nasA.RANGEOVER = (nasA.RANGEOVER & nasB.RANGEOVER);
   nasA.RANGESTEP = (nasA.RANGESTEP & nasB.RANGESTEP);
   nasA.RANGELOW = (nasA.RANGELOW & nasB.RANGELOW);
   nasA.RANGEHIGH = (nasA.RANGEHIGH & nasB.RANGEHIGH);
   nasA.RANGERESET = (nasA.RANGERESET & nasB.RANGERESET);
   nasA.RANGEOLDIC = (nasA.RANGEOLDIC & nasB.RANGEOLDIC);
   nasA.RANGE = (nasA.RANGE & nasB.RANGE);
   nasA.NTST = (nasA.NTST & nasB.NTST);
   nasA.NMAX = (nasA.NMAX & nasB.NMAX);
   nasA.NPR = (nasA.NPR & nasB.NPR);
   nasA.NCOL = (nasA.NCOL & nasB.NCOL);
   nasA.DSMIN = (nasA.DSMIN & nasB.DSMIN);
   nasA.DSMAX = (nasA.DSMAX & nasB.DSMAX);
   nasA.DS = (nasA.DS & nasB.DS);
   nasA.PARMAX = (nasA.PARMAX & nasB.PARMAX);
   nasA.NORMMIN = (nasA.NORMMIN & nasB.NORMMIN);
   nasA.NORMMAX = (nasA.NORMMAX & nasB.NORMMAX);
   nasA.EPSL = (nasA.EPSL & nasB.EPSL);
   nasA.EPSU = (nasA.EPSU & nasB.EPSU);
   nasA.EPSS = (nasA.EPSS & nasB.EPSS);
   nasA.RUNNOW = (nasA.RUNNOW & nasB.RUNNOW);
   nasA.SEC = (nasA.SEC & nasB.SEC);
   nasA.UEC = (nasA.UEC & nasB.UEC);
   nasA.SPC = (nasA.SPC & nasB.SPC);
   nasA.UPC = (nasA.UPC & nasB.UPC);
   nasA.AUTOEVAL = (nasA.AUTOEVAL & nasB.AUTOEVAL);
   nasA.AUTOXMAX = (nasA.AUTOXMAX & nasB.AUTOXMAX);
   nasA.AUTOYMAX = (nasA.AUTOYMAX & nasB.AUTOYMAX);
   nasA.AUTOXMIN = (nasA.AUTOXMIN & nasB.AUTOXMIN);
   nasA.AUTOYMIN = (nasA.AUTOYMIN & nasB.AUTOYMIN);
   nasA.AUTOVAR = (nasA.AUTOVAR & nasB.AUTOVAR);
   nasA.PS_FONT = (nasA.PS_FONT & nasB.PS_FONT);
   nasA.PS_LW = (nasA.PS_LW  & nasB.PS_LW );   
   nasA.PS_FSIZE = (nasA.PS_FSIZE & nasB.PS_FSIZE);
   nasA.PS_COLOR = (nasA.PS_COLOR & nasB.PS_COLOR);
   nasA.FOREVER = (nasA.FOREVER & nasB.FOREVER); 
   nasA.BVP_TOL = (nasA.BVP_TOL & nasB.BVP_TOL); 
   nasA.BVP_EPS = (nasA.BVP_EPS & nasB.BVP_EPS); 
   nasA.BVP_MAXIT = (nasA.BVP_MAXIT & nasB.BVP_MAXIT);
   nasA.BVP_FLAG = (nasA.BVP_FLAG & nasB.BVP_FLAG);
   nasA.SOS = (nasA.SOS & nasB.SOS);
   nasA.FFT = (nasA.FFT & nasB.FFT);
   nasA.HIST = (nasA.HIST & nasB.HIST);
   nasA.PltFmtFlag = (nasA.PltFmtFlag & nasB.PltFmtFlag);
   nasA.ATOLER = (nasA.ATOLER & nasB.ATOLER);
   nasA.MaxEulIter = (nasA.MaxEulIter & nasB.MaxEulIter);
   nasA.EulTol    = (nasA.EulTol    & nasB.EulTol);
   nasA.EVEC_ITER = (nasA.EVEC_ITER & nasB.EVEC_ITER);
   nasA.EVEC_ERR  = (nasA.EVEC_ERR  & nasB.EVEC_ERR);
   nasA.NULL_ERR  = (nasA.NULL_ERR  & nasB.NULL_ERR);
   nasA.NEWT_ERR  = (nasA.NEWT_ERR  & nasB.NEWT_ERR);
   nasA.NULL_HERE = (nasA.NULL_HERE & nasB.NULL_HERE);
}

   

void dump_torus(fp,f)
     FILE *fp;
     int f;
{
  int i;
  char bob[256];
  if(f==READEM)
    fgets(bob,255,fp);
  else
    fprintf(fp,"# Torus information \n");
  io_int(&TORUS,fp,f," Torus flag 1=ON");
  io_double(&TOR_PERIOD,fp,f,"Torus period");
  if(TORUS){
    for(i=0;i<NEQ;i++)
      io_int(&itor[i],fp,f,uvar_names[i]);
  }
}


void load_eqn()
{
 int no_eqn=1,okay=0;
 int i;
 int std=0;
 FILE *fptr;
  init_ar_ic();
  for(i=0;i<MAXODE;i++)
 {
  itor[i]=0;
/*  last_ic[i]=0.0; */
  strcpy(delay_string[i],"0.0");
 }
/* Moved to main
 do_comline(argc,argv); */
 if(strcmp(this_file,"/dev/stdin")==0)std=1;
 struct dirent *dp;
 if (got_file==1&&(std==0) &&(dp=(struct dirent*)opendir(this_file))!=NULL)
  {

  	no_eqn = 1;
	okay=0; 	
	       change_directory(this_file);
	       okay=make_eqn();
	       return;
	
  }
  else
  {
  if(got_file==1&&(fptr=fopen(this_file,"r"))!=NULL)
  {
    if(std==1)sprintf(this_file,"console");
   okay=get_eqn(fptr);
   if(std==0)
     fclose(fptr);
  
   if(okay==1)no_eqn=0;
  }
 }
 if(no_eqn)
   {
     while(okay==0)
       {
       	 struct dirent *dp;
	 char odeclassrm[256];
	 if (getenv("XPPSTART")!=NULL)
	 {
	      
	 	sprintf(odeclassrm,"%s",getenv("XPPSTART"));
	
		if ((dp=(struct dirent*)opendir(odeclassrm))!=NULL)
		{
	 	       change_directory(odeclassrm);
		}
	 }
	 
	 okay=make_eqn();
       } 
   }   
}

/*
load_eqn()
{
 int no_eqn=1,okay=0;
 int i;
 int std=0;
 FILE *fptr;

 init_ar_ic();
  for(i=0;i<MAXODE;i++)
 {
  itor[i]=0;
  strcpy(delay_string[i],"0.0");
 }

 if(strcmp(this_file,"/dev/stdin")==0)std=1;

  struct dirent *dp;
  if(got_file==1)
  {
  	if ((dp=opendir(this_file))!=NULL)
	{
		
		change_directory(this_file);
		no_eqn=1;
	}
	else
	{
		if (fptr=fopen(this_file,"r")!=NULL)
  		{ 
		plintf("Here we are 2\n");
			if(std==1)sprintf(this_file,"console");
   		plintf("Here we are 4\n");
   			okay=get_eqn(fptr);
		plintf("Here we are 3\n");
			if(std==0)
     			  fclose(fptr); 
			  
		        if(okay==1)no_eqn=0;
		}
	}
  }
  if(no_eqn)
   {
     while(okay==0)
       {
	 okay=make_eqn();
       }
     
   }
   
}

*/


void set_X_vals()
{
	/*
	Set up the default look here.
	*/
	
	tfBell=1;
 	/*PaperWhite=0;*/
	/*
	No gradients tends to look cleaner but some
	may prefer gradients improved contrast/readability.
	*/
	/*UserGradients=1;
	*/
	/*fixed is the new X11 default fixed font. 9x15 is dead and gone.
	*/
	if (strlen(big_font_name)==0)
 	{
		strcpy(big_font_name,"fixed");
	}
	
	if (strlen(small_font_name)==0)
 	{
 		strcpy(small_font_name,"6x13");
	}
	
	if (strlen(UserBlack)==0)
 	{
        	sprintf(UserBlack,"#%s","000000");
	}
	
	if (strlen(UserWhite)==0)
 	{
		sprintf(UserWhite,"#%s","EDE9E3");
	}
	
	if (strlen(UserMainWinColor)==0)
 	{
		sprintf(UserMainWinColor,"#%s","808080");
	}
	
	if (strlen(UserDrawWinColor)==0)
 	{
		sprintf(UserDrawWinColor,"#%s","FFFFFF");
	}
	
	if (UserGradients<0)
 	{
		UserGradients=1;
	}
}



void set_all_vals()
{
 int i;
 
 FILE *fp;
 
 if (notAlreadySet.TIMEPLOT){TIMPLOT=1;notAlreadySet.TIMEPLOT=0;};
 if (notAlreadySet.FOREVER){FOREVER=0;notAlreadySet.FOREVER=0;};
 if (notAlreadySet.BVP_TOL){BVP_TOL=1.e-5;notAlreadySet.BVP_TOL=0;};
 if (notAlreadySet.BVP_EPS){BVP_EPS=1.e-5;notAlreadySet.BVP_EPS=0;};
 if (notAlreadySet.BVP_MAXIT){BVP_MAXIT=20;notAlreadySet.BVP_MAXIT=0;};
 if (notAlreadySet.BVP_FLAG){BVP_FLAG=0;notAlreadySet.BVP_FLAG=0;};
 if (notAlreadySet.NMESH){NMESH=40;notAlreadySet.NMESH=0;};
 if (notAlreadySet.NOUT){NJMP=1;notAlreadySet.NOUT=0;};
 if (notAlreadySet.SOS){SOS=0;notAlreadySet.SOS=0;};
 if (notAlreadySet.FFT){FFT=0;notAlreadySet.FFT=0;};
 if (notAlreadySet.HIST){HIST=0;notAlreadySet.HIST=0;};
 if (notAlreadySet.PltFmtFlag){PltFmtFlag=0;notAlreadySet.PltFmtFlag=0;};
 if (notAlreadySet.AXES){AXES=0;notAlreadySet.AXES=0;};
 if (notAlreadySet.TOLER){TOLER=0.001;notAlreadySet.TOLER=0;};
 if (notAlreadySet.ATOLER){ATOLER=0.001;notAlreadySet.ATOLER=0;};
 if (notAlreadySet.MaxEulIter){MaxEulIter=10;notAlreadySet.MaxEulIter=0;}
 if (notAlreadySet.EulTol){EulTol=1.e-7;notAlreadySet.EulTol=0;};
 if (notAlreadySet.DELAY){DELAY=0.0;notAlreadySet.DELAY=0;};
 if (notAlreadySet.DTMIN){HMIN=1e-12;notAlreadySet.DTMIN=0;};
 if (notAlreadySet.EVEC_ITER){EVEC_ITER=100;notAlreadySet.EVEC_ITER=0;};
 if (notAlreadySet.EVEC_ERR){EVEC_ERR=.001;notAlreadySet.EVEC_ERR=0;};
 if (notAlreadySet.NULL_ERR){NULL_ERR=.001;notAlreadySet.NULL_ERR=0;};
 if (notAlreadySet.NEWT_ERR){NEWT_ERR=.001;notAlreadySet.NEWT_ERR=0;};
 if (notAlreadySet.NULL_HERE){NULL_HERE=0;notAlreadySet.NULL_HERE=0;};
 del_stab_flag=DFNORMAL;
 if (notAlreadySet.DTMAX){HMAX=1.000;notAlreadySet.DTMAX=0;};
 if (notAlreadySet.POIMAP){POIMAP=0;notAlreadySet.POIMAP=0;};
 if (notAlreadySet.POIVAR){POIVAR=1;notAlreadySet.POIVAR=0;};
 if (notAlreadySet.POIEXT){POIEXT=0;notAlreadySet.POIEXT=0;};
 if (notAlreadySet.POISGN){POISGN=1;notAlreadySet.POISGN=0;};
 if (notAlreadySet.POIPLN){POIPLN=0.0;notAlreadySet.POIPLN=0;};

 storind=0;
 mov_ind=0;


 STORFLAG=0;


 INFLAG=0;
 oldhp_x=-100000.0 ;
 oldhp_y=-100000.0;
 solver=rung_kut;
 PLOT_3D=0;
 if (notAlreadySet.METHOD){METHOD=3;notAlreadySet.METHOD=0;};
 if (notAlreadySet.XLO){MY_XLO=0.0;x_3d[0]=MY_XLO;notAlreadySet.XLO=0;notAlreadySet.XMIN=0;};
 if (notAlreadySet.XHI){MY_XHI=20.0;x_3d[1]=MY_XHI;notAlreadySet.XHI=0;notAlreadySet.XMAX=0;};
 if (notAlreadySet.YLO){MY_YLO=-1;y_3d[0]=MY_YLO;notAlreadySet.YLO=0;notAlreadySet.YMIN=0;};
 if (notAlreadySet.YHI){MY_YHI=1;y_3d[0]=MY_YHI;notAlreadySet.YHI=0;notAlreadySet.YMAX=0;};
 
 if (notAlreadySet.BOUND){BOUND=100;notAlreadySet.BOUND=0;};
 if (notAlreadySet.MAXSTOR){MAXSTOR=5000;notAlreadySet.MAXSTOR=0;};
 my_pl_wid=10000. ;
 my_pl_ht=7000.  ;

 /* TORUS=0; */ 
 if (notAlreadySet.T0){T0=0.0;notAlreadySet.T0=0;};
 if (notAlreadySet.TRANS){TRANS=0.0;notAlreadySet.TRANS=0;};
 if (notAlreadySet.DT){DELTA_T=.05;notAlreadySet.DT=0;};
 
 if (notAlreadySet.XMIN){x_3d[0]=-12;notAlreadySet.XMIN=0;notAlreadySet.XLO=0;};
 if (notAlreadySet.XMAX){x_3d[1]=12;notAlreadySet.XMAX=0;notAlreadySet.XHI=0;};
 if (notAlreadySet.YMIN){y_3d[0]=-12;notAlreadySet.YMIN=0;notAlreadySet.YLO=0;};
 if (notAlreadySet.YMAX){y_3d[1]=12;notAlreadySet.YMAX=0;notAlreadySet.YHI=0;};
 if (notAlreadySet.ZMIN){z_3d[0]=-12;notAlreadySet.ZMIN=0;};
 if (notAlreadySet.ZMAX){z_3d[1]=12;notAlreadySet.ZMAX=0;};
 
 if (notAlreadySet.TEND){TEND=20.00;notAlreadySet.TEND=0;};
 TOR_PERIOD=6.2831853071795864770;
 if (notAlreadySet.IXPLT){IXPLT=0;notAlreadySet.IXPLT=0;}
 if (notAlreadySet.IYPLT){IYPLT=1;notAlreadySet.IYPLT=0;}
 if (notAlreadySet.IZPLT){IZPLT=1;notAlreadySet.IZPLT=0;}
 
 if (notAlreadySet.NPLOT){
   if (NEQ>2){if(notAlreadySet.IZPLT){IZPLT=2;}}
 NPltV=1;
 for(i=0;i<10;i++){
   IX_PLT[i]=IXPLT;
   IY_PLT[i]=IYPLT;
   IZ_PLT[i]=IZPLT;
   X_LO[i]=0;
   Y_LO[i]=-1;
   X_HI[i]=20;
   Y_HI[i]=1;
 }
 notAlreadySet.NPLOT=0;
 }
 /* internal options go here  */
 set_internopts(NULL);
 

 if((fp=fopen(options,"r"))!=NULL)
 {
  read_defaults(fp);
  fclose(fp);
 } 


 init_range();
 init_trans();
 init_my_aplot();
 init_txtview();

  chk_volterra();  

/*                           */

 if(IZPLT>NEQ)IZPLT=NEQ;
 if(IYPLT>NEQ)IYPLT=NEQ;
 if(IXPLT==0||IYPLT==0)
   TIMPLOT=1;
 else 
   TIMPLOT=0;
 if(x_3d[0]>=x_3d[1]){
   x_3d[0]=-1;
   x_3d[1]=1;
 }
if(y_3d[0]>=y_3d[1]){
   y_3d[0]=-1;
   y_3d[1]=1;
 }
if(z_3d[0]>=z_3d[1]){
   z_3d[0]=-1;
   z_3d[1]=1;
 }
 if(MY_XLO>=MY_XHI){
   MY_XLO=-2.0;
   MY_XHI=2.0;
 }
if(MY_YLO>=MY_YHI){
   MY_YLO=-2.0;
   MY_YHI=2.0;
 }
 if(AXES<5){
   x_3d[0]=MY_XLO;
   y_3d[0]=MY_YLO;
   x_3d[1]=MY_XHI;
   y_3d[1]=MY_YHI;
 } 
 init_stor(MAXSTOR,NEQ+1);
 if(AXES>=5)PLOT_3D=1;
 chk_delay(); /* check for delay allocation */
 alloc_h_stuff();

 alloc_v_memory();  /* allocate stuff for volterra equations */
 alloc_meth();
 arr_ic_start(); /* take care of all predefined array ics */
 

}


void read_defaults(fp)
 FILE *fp;
 {
 char bob[100];
 char *ptr;
 fgets(bob,80,fp);
 ptr=get_first(bob," ");
 if (notAlreadySet.BIG_FONT_NAME)
 {
 	strcpy(big_font_name,ptr);
	notAlreadySet.BIG_FONT_NAME=0;
 }
 
 fgets(bob,80,fp);
 ptr=get_first(bob," ");
 if (notAlreadySet.SMALL_FONT_NAME)
 {
        
 	strcpy(small_font_name,ptr);
	notAlreadySet.SMALL_FONT_NAME=0;
 }
 
 if (notAlreadySet.PaperWhite){fil_int(fp,&PaperWhite);notAlreadySet.PaperWhite=0;};
 if (notAlreadySet.IXPLT){fil_int(fp,&IXPLT);notAlreadySet.IXPLT=0;};
 if (notAlreadySet.IYPLT){fil_int(fp,&IYPLT);notAlreadySet.IYPLT=0;};
 if (notAlreadySet.IZPLT){fil_int(fp,&IZPLT);notAlreadySet.IZPLT=0;};
 if (notAlreadySet.AXES){fil_int(fp,&AXES);notAlreadySet.PaperWhite=0;};
 if (notAlreadySet.NOUT){fil_int(fp,&NJMP);notAlreadySet.NOUT=0;};
 if (notAlreadySet.NMESH){fil_int(fp,&NMESH);notAlreadySet.NMESH=0;};
 if (notAlreadySet.METHOD){fil_int(fp,&METHOD);notAlreadySet.METHOD=0;};

 if (notAlreadySet.TIMEPLOT){fil_int(fp,&TIMPLOT);notAlreadySet.TIMEPLOT=0;};
 if (notAlreadySet.MAXSTOR){fil_int(fp,&MAXSTOR);notAlreadySet.MAXSTOR=0;};
 if (notAlreadySet.TEND){fil_flt(fp,&TEND);notAlreadySet.TEND=0;};
 if (notAlreadySet.DT){fil_flt(fp,&DELTA_T);notAlreadySet.DT=0;};
 if (notAlreadySet.T0){fil_flt(fp,&T0);notAlreadySet.T0=0;};
 if (notAlreadySet.TRANS){fil_flt(fp,&TRANS);notAlreadySet.TRANS=0;};
 if (notAlreadySet.BOUND){fil_flt(fp,&BOUND);notAlreadySet.BOUND=0;};
 if (notAlreadySet.DTMIN){fil_flt(fp,&HMIN);notAlreadySet.DTMIN=0;};
 if (notAlreadySet.DTMAX){fil_flt(fp,&HMAX);notAlreadySet.DTMIN=0;};
 if (notAlreadySet.TOLER){fil_flt(fp,&TOLER);notAlreadySet.TOLER=0;};
 if (notAlreadySet.DELAY){fil_flt(fp,&DELAY);notAlreadySet.DELAY=0;};
 if (notAlreadySet.XLO){fil_flt(fp,&MY_XLO);notAlreadySet.XLO=0;};
 if (notAlreadySet.XHI){fil_flt(fp,&MY_XHI);notAlreadySet.XHI=0;};
 if (notAlreadySet.YLO){fil_flt(fp,&MY_YLO);notAlreadySet.YLO=0;};
 if (notAlreadySet.YHI){fil_flt(fp,&MY_YHI);notAlreadySet.YHI=0;};

 
}

void fil_flt(fpt,val)
FILE *fpt;
double *val;
{
 char bob[80];
 fgets(bob,80,fpt);
 *val=atof(bob);
}

void fil_int(fpt,val)
int *val;
FILE *fpt;
{
 char bob[80];
 fgets(bob,80,fpt);
 *val=atoi(bob);
}



/* here is some new code for internal set files:
   format of the file is a long string of the form:
   { x=y, z=w, q=p , .... }
*/



void add_intern_set(name,does)
     char *name,*does;
{
  char bob[256],ch;
  int i,n,j=Nintern_set,k=0;
  if(Nintern_set>=MAX_INTERN_SET){
   plintf(" %s not added -- too many must be less than %d \n",
	   name,MAX_INTERN_SET);
    return;
  }
  intern_set[j].use=1;
  n=strlen(name);
  intern_set[j].name=(char *)malloc((n+1));
  strcpy(intern_set[j].name,name);
  n=strlen(does);
  bob[0]='$';
  bob[1]=' ';
  k=2;
  for(i=0;i<n;i++){
    ch=does[i];
    if(ch==','){
      bob[k]=' ';
      k++;
    }
    if(ch=='}'||ch=='{')
      continue;
    if(ch!=','){
      bob[k]=ch;
      k++;
    }
  }
  bob[k]=0;
  intern_set[j].does=(char *)malloc(n+3);
  strcpy(intern_set[j].does,bob);
 plintf(" added %s doing %s \n",
	 intern_set[j].name,intern_set[j].does);
  Nintern_set++;
}
      

void extract_action(char *ptr)
{
  char name[256],value[256];
 char tmp[2048];
  char *junk,*mystring;
  /* plintf("ptr=%s \n",ptr);  */
  strcpy(tmp,ptr);
  junk=get_first(tmp," ");
  if (junk == NULL)
  {
  	/*No more tokens--should this throw an error?*/
  }
  
  while((mystring=get_next(" ,;\n"))!=NULL){
   split_apart(mystring,name,value);
      if(strlen(name)>0&&strlen(value)>0)
       do_intern_set(name,value);
    } 
}

void extract_internset(j)
     int j;
{
  extract_action(intern_set[j].does);
}

void do_intern_set(name1,value)
     char *name1,*value;
{
  int i;
  char name[20];
  convert(name1,name);

  i=find_user_name(IC,name);
  if(i>-1){
    last_ic[i]=atof(value);
  }
  else {
    i=find_user_name(PARAM,name);
    if(i>-1){
      set_val(name,atof(value));
    }
    else {
      /*     set_option(name,value,0,NULL); */
      set_option(name,value,1,NULL);
   }
  }
 alloc_meth();
 do_meth();
}
/*  ODE options stuff  here !!   */

int msc(s1,s2)
     char *s1,*s2;
{

 int n=strlen(s1),i;
 if(strlen(s2)<n)return(0);
 for(i=0;i<n;i++)
   if(s1[i]!=s2[i])return(0);
 return(1);
}  
  
void set_internopts(OptionsSet *mask)
{
  int i;
  char *ptr,name[20],value[80],*junk,*mystring;
  if(Nopts==0)return;
 /*  parsem here   */
  for(i=0;i<Nopts;i++){
    ptr=interopt[i];
    junk=get_first(ptr," ,");
    if (junk == NULL)
    {
    	/*No more tokens.  Should this throw an error?*/
    }	
    while((mystring=get_next(" ,\n\r"))!=NULL)
    {
      split_apart(mystring,name,value);
      if(strlen(name)>0&&strlen(value)>0)
      {
        /*
	if (strcmp("mwcolor",name)==0)
	{
		if (strlen(UserMainWinColor)!=0)
		{
		 	continue;
		}
	}
	
	if (strcmp("dwcolor",name)==0)
	{
		if (strlen(UserDrawWinColor)!=0)
		{
		 	continue;
		}
	}
	
	if (strcmp("forecolor",name)==0)
	{
		if (strlen(UserWhite)!=0)
		{
		 	continue;
		}
	}
	
	if (strcmp("backcolor",name)==0)
	{
		if (strlen(UserBlack)!=0)
		{
		 	continue;
		}
	}
	
	if (strcmp("backimage",name)==0)
	{
		if (strlen(UserBGBitmap)!=0)
		{
		 	continue;
		}
	}
	
	if (strcmp("smallfont",name)==0)
	{
		if (strlen(small_font_name)!=0)
		{
		 	continue;
		}
	}
	
	if (strcmp("bigfont",name)==0)
	{
		if (strlen(big_font_name)!=0)
		{
		 	continue;
		}
	}
	*/
	set_option(name,value,0,mask);
      }
    }
  }
 
  for(i=0;i<Nopts;i++)
  {
    free(interopt[i]);
  }  
  Nopts = 0;  
}

void set_internopts_xpprc_and_comline()
{
  int i;
  char *ptr,name[20],value[80],*junk,*mystring;
  if(Nopts==0)return;
 /*  parsem here   */
 /*Check for QUIET and LOGFILE options first...*/
  char intrnoptcpy[255]; /*Must use copy to avoid side effects of strtok used in get_first below*/
  for(i=0;i<Nopts;i++){
    strcpy(intrnoptcpy,interopt[i]);
    ptr=intrnoptcpy;
    junk=get_first(ptr," ,");  
    if (junk == NULL)
    {
    	/*No more tokens.  Should this throw an error?*/
    }	
    while((mystring=get_next(" ,\n\r"))!=NULL)
    {
      split_apart(mystring,name,value);
      strupr(name); 
      
      if (strlen(name)==5)
      {
      	      strupr(name);
	      if(strcmp(name,"QUIET")==0)
	      {
		set_option(name,value,0,NULL);
	      }
      }
      else if (strlen(name)==7)
      {	
	      strupr(name);
	      
	      if(strcmp(name,"LOGFILE")==0)
	      {
	      	      set_option(name,value,0,NULL);
	      }
    }
  }
  }
  
  /*We make a BOOLEAN MASK using the current OptionsSet*/
  /*This allows options to be overwritten multiple times within .xpprc
  but prevents overwriting across comline, .xpprc etc.
  */ 
  OptionsSet *tempNAS = (OptionsSet*)malloc(sizeof(OptionsSet));
  *tempNAS = notAlreadySet;
  
  for(i=0;i<Nopts;i++){
    ptr=interopt[i];
    junk=get_first(ptr," ,");
    while((mystring=get_next(" ,\n\r"))!=NULL)
    {
      split_apart(mystring,name,value);
      if(strlen(name)>0&&strlen(value)>0)
      {
        set_option(name,value,0,tempNAS);
      }	
    }
  }
  free(tempNAS);
   
  /*
  We leave a fresh start for options specified in the ODE file.
  */
  for(i=0;i<Nopts;i++)
  {
  	free(interopt[i]);
  } 
  
  Nopts=0;
}


void split_apart(bob, name,value)
char *bob,*name,*value;
{
 int k,i,l;


 l=strlen(bob);
 k=strcspn(bob,"=");
 if(k==l)
 {
  value[0]=0;
  strcpy(name,bob);
  }
  else
  {
  strncpy(name,bob,k);
  name[k]='\0';
  for(i=k+1;i<l;i++)value[i-k-1]=bob[i];
  value[l-k-1]='\0';
    }

}




void check_for_xpprc()
{
  FILE *fp;
  char rc[256];
  char bob[256];
  sprintf(rc,"%s/.xpprc",getenv("HOME"));
  fp=fopen(rc,"r");
  if(fp==NULL){
    /*   plintf("Didnt find rc \n"); */
    return;
  }
  while(!feof(fp)){
    bob[0]='\0';
    fgets(bob,255,fp);
    if(bob[0]=='@'){
      stor_internopts(bob);

    }
  }
  fclose(fp);
}


void stor_internopts(s1)
     char *s1;
{
  int n=strlen(s1);
  if(Nopts>MAXOPT){
   plintf("WARNING -- to many options set %s ignored\n",s1);
    return;
  }
  interopt[Nopts]=(char *)malloc(n+1);
  sprintf(interopt[Nopts],"%s",s1);
  Nopts++;

}
  


void set_option(s1,s2,force,mask)
     char *s1,*s2;
     int force;
     OptionsSet *mask;
{
  int i,j,f;
 char xx[4],yy[4],zz[4];
 char xxl[6],xxh[6],yyl[6],yyh[6];
 static char mkey[]="demragvbqsc582y";
 static char Mkey[]="DEMRAGVBQSC582Y";
 strupr(s1);
 if(msc("QUIET",s1)){
   if(!(msc(s2,"0")||msc(s2,"1")))
   {
   	plintf("QUIET option must be 0 or 1.\n");
	exit(-1);
   }
   if (OVERRIDE_QUIET==0)/*Will be 1 if -quiet was specified on the command line.*/
   {
   	XPPVERBOSE=(atoi(s2)==0);
   }
   return;
 }
 if(msc("LOGFILE",s1)){
   if (OVERRIDE_LOGFILE==0) /*Will be 1 if -logfile was specified on the command line.*/
   {
      if (logfile != NULL)       
      { 		         
     	  fclose(logfile);       
      } 		         
      logfile=fopen(s2,"w");     
   }
   return;
 }
 if(msc("BELL",s1)){
   if(!(msc(s2,"0")||msc(s2,"1")))
   {
   	plintf("BELL option must be 0 or 1.\n");
	exit(-1);
   }
   tfBell=atoi(s2);
   return;
 }
 if(msc("BUT",s1)){
    add_user_button(s2);
    return;
  }
 if((msc("BIGFONT",s1))||(msc("BIG",s1))){
    if ((notAlreadySet.BIG_FONT_NAME||force) || ((mask!=NULL)&&(mask->BIG_FONT_NAME==1)))
    {
    	strcpy(big_font_name,s2);
	notAlreadySet.BIG_FONT_NAME=0;
    }
    return;
  }
  if((msc("SMALLFONT",s1))||(msc("SMALL",s1))){;
    if ((notAlreadySet.SMALL_FONT_NAME||force) || ((mask!=NULL)&&(mask->SMALL_FONT_NAME==1)))
    {
    	strcpy(small_font_name,s2);
	notAlreadySet.SMALL_FONT_NAME=0;
    }
    return;
  }
  if(msc("FORECOLOR",s1)){
    if ((notAlreadySet.UserBlack||force) || ((mask!=NULL)&&(mask->UserBlack==1)))
    {
    	sprintf(UserBlack,"#%s",s2);
	notAlreadySet.UserBlack=0;
    }
    return;
  }
  if(msc("BACKCOLOR",s1)){
    if ((notAlreadySet.UserWhite||force) || ((mask!=NULL)&&(mask->UserWhite==1)))
    {
    	sprintf(UserWhite,"#%s",s2);
	notAlreadySet.UserWhite=0;
    }
    return;
  }
  if(msc("MWCOLOR",s1)){
    if ((notAlreadySet.UserMainWinColor||force) || ((mask!=NULL)&&(mask->UserMainWinColor==1)))
    {
      /* printf("Setting MWCOLOR=%s\n",s2); */
        sprintf(UserMainWinColor,"#%s",s2);
	notAlreadySet.UserMainWinColor=0;
    }
    return;
  }
  if(msc("DWCOLOR",s1)){
    if ((notAlreadySet.UserDrawWinColor||force) || ((mask!=NULL)&&(mask->UserDrawWinColor==1)))
    {
    	sprintf(UserDrawWinColor,"#%s",s2);
	notAlreadySet.UserDrawWinColor=0;
    }
    return;
  }
  if(msc("GRADS",s1)){
    if ((notAlreadySet.UserGradients||force) || ((mask!=NULL)&&(mask->UserGradients==1)))
    {
	    if(!(msc(s2,"0")||msc(s2,"1")))
	    {
   		 plintf("GRADS option must be 0 or 1.\n");
		 exit(-1);
	    }
	    UserGradients=atoi(s2);
	    notAlreadySet.UserGradients=0;
    }
    return;
  }



  if(msc("PLOTFMT",s1)){
    if ((notAlreadySet.PLOTFORMAT||force) || ((mask!=NULL)&&(mask->PLOTFORMAT==1)))
    {
    	strcpy(PlotFormat,s2);
	notAlreadySet.PLOTFORMAT=0;
    }
    return;
  }

  

  if(msc("BACKIMAGE",s1)){
    if ((notAlreadySet.UserBGBitmap||force) || ((mask!=NULL)&&(mask->UserBGBitmap==1)))
    {
    	strcpy(UserBGBitmap,s2);
	notAlreadySet.UserBGBitmap=0;
    }
    return;
  }
  if(msc("WIDTH",s1)){
    if ((notAlreadySet.UserMinWidth||force)|| ((mask!=NULL)&&(mask->UserMinWidth==1)))
    {
       UserMinWidth=atoi(s2);
       notAlreadySet.UserMinWidth=0;
    }
    return;
  }
  if(msc("HEIGHT",s1)){
    if ((notAlreadySet.UserMinHeight||force) || ((mask!=NULL)&&(mask->UserMinHeight==1)))
    {
         UserMinHeight=atoi(s2);
	 notAlreadySet.UserMinHeight=0;
    }
    return;
  }
  if(msc("YNC",s1)){
    if ((notAlreadySet.YNullColor||force) || ((mask!=NULL)&&(mask->YNullColor==1)))
    {
	  i=atoi(s2);
	  if(i>-1&&i<11)
	  {
	   YNullColor=i;
	  }
	   notAlreadySet.YNullColor=0;
    }
  return;
  }
if(msc("XNC",s1)){
    if ((notAlreadySet.XNullColor||force) || ((mask!=NULL)&&(mask->XNullColor==1)))
    {
	    i=atoi(s2);
	  if(i>-1&&i<11)
	  {
	   XNullColor=i; 
	   notAlreadySet.XNullColor=0;
	  }
	  
    }
  return;
  }

if(msc("SMC",s1)){

    if ((notAlreadySet.StableManifoldColor||force) || ((mask!=NULL)&&(mask->StableManifoldColor==1)))
    {
  
	  i=atoi(s2);
	  if(i>-1&&i<11)
	  {
	   StableManifoldColor=i;
	   notAlreadySet.StableManifoldColor=0;
	  }
    }
  return;
  }
if(msc("UMC",s1)){
    if ((notAlreadySet.UnstableManifoldColor||force) || ((mask!=NULL)&&(mask->UnstableManifoldColor==1)))
    {
	    i=atoi(s2);
	    if(i>-1&&i<11)
	    {
	     UnstableManifoldColor=i;
	     notAlreadySet.UnstableManifoldColor=0;
	    }
    }
   return;
  }

  if(msc("LT",s1)){
     if ((notAlreadySet.START_LINE_TYPE||force) || ((mask!=NULL)&&(mask->START_LINE_TYPE==1)))
     {
     	
	    i=atoi(s2);
	    if(i<2&&i>-6)
	    {  
	      START_LINE_TYPE=i; 
	      reset_all_line_type();
	      notAlreadySet.START_LINE_TYPE=0;
	      }
     }
     return;
  }
  if(msc("SEED",s1)){ 
     if ((notAlreadySet.RandSeed||force) || ((mask!=NULL)&&(mask->RandSeed==1)))
     {
	    i=atoi(s2);
	    if(i>=0){
	      RandSeed=i;
	      nsrand48(RandSeed);  
	      notAlreadySet.RandSeed=0;
	    }
     }
    return;
  }
 if(msc("BACK",s1)){
   if ((notAlreadySet.PaperWhite||force) || ((mask!=NULL)&&(mask->PaperWhite==1)))
   {
	   if(s2[0]=='w'||s2[0]=='W')
	   {
	   	PaperWhite=1;
	   }
	   else 
	   {  
	   	PaperWhite=0;
	   }
	   notAlreadySet.PaperWhite=0;
   }
    return;
  }
 if(msc("COLORMAP",s1)){
     if ((notAlreadySet.COLORMAP||force) || ((mask!=NULL)&&(mask->COLORMAP==1)))
     {
   		i=atoi(s2);
   		if(i<7)custom_color=i;
		notAlreadySet.COLORMAP=0;

     }
   return;
 }
   if(msc("NPLOT",s1)){
     if ((notAlreadySet.NPLOT||force) || ((mask!=NULL)&&(mask->NPLOT==1)))
     {
    	NPltV=atoi(s2);
	notAlreadySet.NPLOT=0;
     }
    return;
  }

   if(msc("DLL_LIB",s1)){
      if ((notAlreadySet.DLL_LIB||force) || ((mask!=NULL)&&(mask->DLL_LIB==1)))
     {
     sprintf(dll_lib,"%s",s2);
     dll_flag+=1;
     notAlreadySet.DLL_LIB=0;
     }
     return;
   }
   if(msc("DLL_FUN",s1)){
     if ((notAlreadySet.DLL_FUN||force) || ((mask!=NULL)&&(mask->DLL_FUN==1)))
     {
     	sprintf(dll_fun,"%s",s2);
     	dll_flag+=2;
     	notAlreadySet.DLL_FUN=0;
     }
     return;
   }
   /* can now initialize several plots */
   if(msc("SIMPLOT",s1)){
     SimulPlotFlag=1;
     return;
   }
   if(msc("MULTIWIN",s1)){
     MultiWin=1;
     return;
   }
 for(j=2;j<=8;j++){
      sprintf(xx,"XP%d",j);
      sprintf(yy,"YP%d",j);
      sprintf(zz,"ZP%d",j);
      sprintf(xxh,"XHI%d",j);    
      sprintf(xxl,"XLO%d",j);    
      sprintf(yyh,"YHI%d",j);    
      sprintf(yyl,"YLO%d",j);    
    if(msc(xx,s1)){
    find_variable(s2,&i);
    if(i>-1)IX_PLT[j]=i;
    return;
  }
   if(msc(yy,s1)){
     find_variable(s2,&i);
    if(i>-1)IY_PLT[j]=i;
    return;
  }
   if(msc(zz,s1)){
     find_variable(s2,&i);
    if(i>-1)IZ_PLT[j]=i;
    return;
  }
   if(msc(xxh,s1)){
     X_HI[j]=atof(s2);
     return;
   }
   if(msc(xxl,s1)){
     X_LO[j]=atof(s2);
     return;
   }
if(msc(yyh,s1)){
     Y_HI[j]=atof(s2);
     return;
   }
if(msc(yyl,s1)){
     Y_LO[j]=atof(s2);
     return;
   }
 }
   if(msc("XP",s1)){
     if ((notAlreadySet.XP||force) || ((mask!=NULL)&&(mask->XP==1)))
     {
    	find_variable(s2,&i);
    	if(i>-1)IXPLT=i;
	notAlreadySet.XP=0;
	notAlreadySet.IXPLT=0;
     }
    return;
  }
   if(msc("YP",s1)){
     if ((notAlreadySet.YP||force) || ((mask!=NULL)&&(mask->YP==1)))
     {
     	find_variable(s2,&i);
    	if(i>-1)IYPLT=i;
	notAlreadySet.YP=0;
	notAlreadySet.IYPLT=0;
     }
    return;
  }
   if(msc("ZP",s1)){
     if ((notAlreadySet.ZP||force) || ((mask!=NULL)&&(mask->ZP==1)))
     {
     	find_variable(s2,&i);
    	if(i>-1)IZPLT=i;

     	notAlreadySet.ZP=0;
	notAlreadySet.IZPLT=0;
     }
    return;
  }
   if(msc("AXES",s1)){
     if ((notAlreadySet.AXES||force) || ((mask!=NULL)&&(mask->AXES==1)))
     {
	 if(s2[0]=='3')
	 {
	   AXES=5;
	 }
	 else 
	 {
	   AXES=0;
	 } 
        
	 notAlreadySet.AXES=0;
     }
    return;
  }

   if(msc("NJMP",s1)){
     if ((notAlreadySet.NOUT||force) || ((mask!=NULL)&&(mask->NOUT==1)))
     {
    	NJMP=atoi(s2);
        notAlreadySet.NOUT=0;
     }
    return;
  }
  if(msc("NOUT",s1)){
     if ((notAlreadySet.NOUT||force) || ((mask!=NULL)&&(mask->NOUT==1)))
     {
      NJMP=atoi(s2);
      notAlreadySet.NOUT=0;
     }
    return;
  }
   if(msc("NMESH",s1)){
     if ((notAlreadySet.NMESH||force) || ((mask!=NULL)&&(mask->NMESH==1)))
     {
    	NMESH=atoi(s2);
	notAlreadySet.NMESH=0;
     }
    return;
  }
   if(msc("METH",s1)){
     if ((notAlreadySet.METHOD||force) || ((mask!=NULL)&&(mask->METHOD==1)))
     {
    for(i=0;i<15;i++)
      if(s2[0]==mkey[i]||s2[0]==Mkey[i])
	METHOD=i;
      
       notAlreadySet.METHOD=0;
     }
    return;
  }
   if(msc("VMAXPTS",s1)){
     if ((notAlreadySet.VMAXPTS||force) || ((mask!=NULL)&&(mask->VMAXPTS==1)))
     {
     	MaxPoints=atoi(s2);
	notAlreadySet.VMAXPTS=0;
     
     }
     return;
   }
   if(msc("MAXSTOR",s1)){ 
     if ((notAlreadySet.MAXSTOR||force) || ((mask!=NULL)&&(mask->MAXSTOR==1)))
     {
    	MAXSTOR=atoi(s2);
        notAlreadySet.MAXSTOR=0;
     } 
    return;
  }
   if(msc("TOR_PER",s1)){
     if ((notAlreadySet.TOR_PER||force) || ((mask!=NULL)&&(mask->TOR_PER==1)))
     {
     	TOR_PERIOD=atof(s2);
     	TORUS=1;
	notAlreadySet.TOR_PER=0;
     }
     return;
   }
   if(msc("JAC_EPS",s1)){
     if ((notAlreadySet.JAC_EPS||force) || ((mask!=NULL)&&(mask->JAC_EPS==1)))
     {
     	NEWT_ERR=atof(s2);
        notAlreadySet.JAC_EPS=0;
     }
     return;
   }
   if(msc("NEWT_TOL",s1)){
     if ((notAlreadySet.NEWT_TOL||force) || ((mask!=NULL)&&(mask->NEWT_TOL==1)))
     {
     	EVEC_ERR=atof(s2);
	notAlreadySet.NEWT_TOL=0;
     
     }
     return;
   }
   if(msc("NEWT_ITER",s1)){
     if ((notAlreadySet.NEWT_ITER||force) || ((mask!=NULL)&&(mask->NEWT_ITER==1)))
     {
     	EVEC_ITER=atoi(s2);
	notAlreadySet.NEWT_ITER=0;
     }
     return;
   }
  if(msc("FOLD",s1)){
     if ((notAlreadySet.FOLD||force) || ((mask!=NULL)&&(mask->FOLD==1)))
     {
     find_variable(s2,&i);
     if(i>0){
       itor[i-1]=1;
      TORUS=1;
     }
     
     }
     return;
   }
   if(msc("TOTAL",s1)){
    if ((notAlreadySet.TEND||force) || ((mask!=NULL)&&(mask->TEND==1)))
     {
    	TEND=atof(s2);
	notAlreadySet.TEND=0;
    }
    return;
  }
  if(msc("DTMIN",s1)){
     if ((notAlreadySet.DTMIN||force) || ((mask!=NULL)&&(mask->DTMIN==1)))
     {
    	HMIN=atof(s2);
         notAlreadySet.DTMIN=0;
     }
    return;
  }
  if(msc("DTMAX",s1)){
     if ((notAlreadySet.DTMAX||force) || ((mask!=NULL)&&(mask->DTMAX==1)))
     {
    	HMAX=atof(s2);
	notAlreadySet.DTMAX=0;
      }
    return;
  }
   if(msc("DT",s1)){
     if ((notAlreadySet.DT||force) || ((mask!=NULL)&&(mask->DT==1)))
     {
    	DELTA_T=atof(s2);
	notAlreadySet.DT=0;
     }
    return;
  }
   if(msc("T0",s1)){
     if ((notAlreadySet.T0||force) || ((mask!=NULL)&&(mask->T0==1)))
     { 
    	T0=atof(s2);
        notAlreadySet.T0=0;
     }
    return;
  }
   if(msc("TRANS",s1)){
     if ((notAlreadySet.TRANS||force) || ((mask!=NULL)&&(mask->TRANS==1)))
     {
     	TRANS=atof(s2);
        notAlreadySet.TRANS=0;
     }
    return;
  }
   if(msc("BOUND",s1)){
     if ((notAlreadySet.BOUND||force) || ((mask!=NULL)&&(mask->BOUND==1)))
     {
       BOUND=atof(s2);
       notAlreadySet.BOUND=0;
     }
    return;
  }
   if(msc("ATOL",s1)){
     if ((notAlreadySet.ATOLER||force) || ((mask!=NULL)&&(mask->ATOLER==1)))
     {
     	ATOLER=atof(s2);
        notAlreadySet.ATOLER=0;
     }
     return;
   }
   if(msc("TOL",s1)){
     if ((notAlreadySet.TOLER||force) || ((mask!=NULL)&&(mask->TOLER==1)))
     {
	TOLER=atof(s2);
	notAlreadySet.TOLER=0;
     }
    return;
  }
    
   if(msc("DELAY",s1)){
     if ((notAlreadySet.DELAY||force) || ((mask!=NULL)&&(mask->DELAY==1)))
     {
    	DELAY=atof(s2);
	notAlreadySet.DELAY=0;
     }
    return;
  }
   if(msc("BANDUP",s1)){
     if ((notAlreadySet.BANDUP||force) || ((mask!=NULL)&&(mask->BANDUP==1)))
     {
     	cv_bandflag=1;
     	cv_bandupper=atoi(s2);
     	notAlreadySet.BANDUP=0;
     }
     return;
   }
  if(msc("BANDLO",s1)){
     if ((notAlreadySet.BANDLO||force) || ((mask!=NULL)&&(mask->BANDLO==1)))
     {
     	cv_bandflag=1;
     	cv_bandlower=atoi(s2);
     	notAlreadySet.BANDLO=0;
     }
     return;
   }
  
  if(msc("PHI",s1)){
     if ((notAlreadySet.PHI||force) || ((mask!=NULL)&&(mask->PHI==1)))
     {
    	PHI0=atof(s2);
	notAlreadySet.PHI=0;
     }
    return;
  }
   if(msc("THETA",s1)){
     if ((notAlreadySet.THETA||force) || ((mask!=NULL)&&(mask->THETA==1)))
     {
    	THETA0=atof(s2);
	notAlreadySet.THETA=0;
     }
    return;
  }
   if(msc("XLO",s1)){
     if ((notAlreadySet.XLO||force) || ((mask!=NULL)&&(mask->XLO==1)))
     {
    	MY_XLO=atof(s2);
	notAlreadySet.XLO=0;
     }
    return;
  }
   if(msc("YLO",s1)){
    if ((notAlreadySet.YLO||force) || ((mask!=NULL)&&(mask->YLO==1)))
    {
    	MY_YLO=atof(s2);
	notAlreadySet.YLO=0;
    }
    return;
  }
  
   if(msc("XHI",s1)){
    if ((notAlreadySet.XHI||force) || ((mask!=NULL)&&(mask->XHI==1)))
    {
    	MY_XHI=atof(s2);
        notAlreadySet.XHI=0;
    }
    return;
  }
   if(msc("YHI",s1)){
     if ((notAlreadySet.YHI||force) || ((mask!=NULL)&&(mask->YHI==1)))
     {
    	MY_YHI=atof(s2);
        notAlreadySet.YHI=0;
     }
    return;
  }
   if(msc("XMAX",s1)){
     if ((notAlreadySet.XMAX||force) || ((mask!=NULL)&&(mask->XMAX==1)))
     {
    	x_3d[1]=atof(s2);
	notAlreadySet.XMAX=0;
     
     }
    return;
  }
   if(msc("YMAX",s1)){
     if ((notAlreadySet.YMAX||force) || ((mask!=NULL)&&(mask->YMAX==1)))
     {
        y_3d[1]=atof(s2);
	notAlreadySet.YMAX=0;
     }
    return;
  }
   if(msc("ZMAX",s1)){
     if ((notAlreadySet.ZMAX||force) || ((mask!=NULL)&&(mask->ZMAX==1)))
     {
        z_3d[1]=atof(s2);
	notAlreadySet.ZMAX=0;
     }
    return;
  }
   if(msc("XMIN",s1)){
     /*  printf("Trying to set XMIN %d =%s\n",notAlreadySet.XMIN,s2); */
     if ((notAlreadySet.XMIN||force) || ((mask!=NULL)&&(mask->XMIN==1)))
     {
        x_3d[0]=atof(s2);
	notAlreadySet.XMIN=0; 
	if ((notAlreadySet.XLO||force) || ((mask!=NULL)&&(mask->XLO==1)))
	{
    	   MY_XLO=atof(s2);
	   notAlreadySet.XLO=0;
	}
     }
    return;
  }
   if(msc("YMIN",s1)){
     if ((notAlreadySet.YMIN||force) || ((mask!=NULL)&&(mask->YMIN==1)))
     {
    	y_3d[0]=atof(s2);
	notAlreadySet.YMIN=0;
	if ((notAlreadySet.YLO||force) || ((mask!=NULL)&&(mask->YLO==1)))
	{
    	   MY_YLO=atof(s2);
	   notAlreadySet.YLO=0;
	}
     }
    return;
  }
 if(msc("ZMIN",s1)){
     if ((notAlreadySet.ZMIN||force) || ((mask!=NULL)&&(mask->ZMIN==1)))
     {
    	z_3d[0]=atof(s2);
	notAlreadySet.ZMIN=0;
     }
    return;
  }

 if(msc("POIMAP",s1)){
     if ((notAlreadySet.POIMAP||force) || ((mask!=NULL)&&(mask->POIMAP==1)))
     {
   	if(s2[0]=='m'||s2[0]=='M')POIMAP=2;
   	if(s2[0]=='s'||s2[0]=='S')POIMAP=1;
   	if(s2[0]=='p'||s2[0]=='P')POIMAP=3;
   	notAlreadySet.POIMAP=0;
   }
   return;
 }

 if(msc("POIVAR",s1)){
     if ((notAlreadySet.POIVAR||force) || ((mask!=NULL)&&(mask->POIVAR==1)))
     {
    	find_variable(s2,&i);
    	if(i>-1)POIVAR=i;
	
	notAlreadySet.POIVAR=0;
     
     }
    return;
  }
 if(msc("OUTPUT",s1)){
     if ((notAlreadySet.OUTPUT||force) || ((mask!=NULL)&&(mask->OUTPUT==1)))
     {
   	strcpy(batchout,s2);
	notAlreadySet.OUTPUT=0;
     }
   return;
 }
  
 if(msc("POISGN",s1)){
     if ((notAlreadySet.POISGN||force) || ((mask!=NULL)&&(mask->POISGN==1)))
     {
   	POISGN=atoi(s2);
	notAlreadySet.POISGN=0;
     }
   return;
 }
 
 if(msc("POISTOP",s1)){
     if ((notAlreadySet.POISTOP||force) || ((mask!=NULL)&&(mask->POISTOP==1)))
     {
   	SOS=atoi(s2);
	notAlreadySet.POISTOP=0;
     }
   return;
 }
 if(msc("STOCH",s1)){
     if ((notAlreadySet.STOCH||force)|| ((mask!=NULL)&&(mask->STOCH==1)))
     {
   	STOCH_FLAG=atoi(s2);
	notAlreadySet.STOCH=0;
     }
   return;
 }
 if(msc("POIPLN",s1)){
     if ((notAlreadySet.POIPLN||force)|| ((mask!=NULL)&&(mask->POIPLN==1)))
     {
   	POIPLN=atof(s2);
	notAlreadySet.POIPLN=0;
     }
   return;
 }
  
 

 if(msc("RANGEOVER",s1)){
     if ((notAlreadySet.RANGEOVER||force)|| ((mask!=NULL)&&(mask->RANGEOVER==1)))
     {
    	strcpy(range.item,s2);
	notAlreadySet.RANGEOVER=0;
     }

    return;
  }
 if(msc("RANGESTEP",s1)){
     if ((notAlreadySet.RANGESTEP||force)|| ((mask!=NULL)&&(mask->RANGESTEP==1)))
     {
        
   	range.steps=atoi(s2);
	notAlreadySet.RANGESTEP=0;
     }
   return;
 }
  
 if(msc("RANGELOW",s1)){
     if ((notAlreadySet.RANGELOW||force)|| ((mask!=NULL)&&(mask->RANGELOW==1)))
     {
   	range.plow=atof(s2);
   	notAlreadySet.RANGELOW=0;
     }

   return;
 }

 if(msc("RANGEHIGH",s1)){
     if ((notAlreadySet.RANGEHIGH||force)|| ((mask!=NULL)&&(mask->RANGEHIGH==1)))
     {
   	range.phigh=atof(s2);
	notAlreadySet.RANGEHIGH=0;
     }
   return;
 }
 
 if(msc("RANGERESET",s1)){
     if ((notAlreadySet.RANGERESET||force)|| ((mask!=NULL)&&(mask->RANGERESET==1)))
     {
	 if(s2[0]=='y'||s2[0]=='Y')
	 {
	  range.reset=1;
	 }
	 else
	 {
	  range.reset=0;
	 } 
	  notAlreadySet.RANGERESET=0;
     }
  	return;
   }

 if(msc("RANGEOLDIC",s1)){
     if ((notAlreadySet.RANGEOLDIC||force)|| ((mask!=NULL)&&(mask->RANGEOLDIC==1)))
     {
  	if(s2[0]=='y'||s2[0]=='Y')
	{
   		range.oldic=1;
   	}
	else
	{ 
   		range.oldic=0;
	}
	
   	notAlreadySet.RANGEOLDIC=0;
     }
      return;
 }
 
   
 if(msc("RANGE",s1)){
     if ((notAlreadySet.RANGE||force)|| ((mask!=NULL)&&(mask->RANGE==1)))
     {
   	batch_range=atoi(s2);
	notAlreadySet.RANGE=0;
     }
   return;
 }
 
 if(msc("NTST",s1)){
     if ((notAlreadySet.NTST||force)|| ((mask!=NULL)&&(mask->NTST==1)))
     {
   	auto_ntst=atoi(s2);
	notAlreadySet.NTST=0;
     }
   return;
 }
if(msc("NMAX",s1)){
   if ((notAlreadySet.NMAX||force)|| ((mask!=NULL)&&(mask->NMAX==1)))
   {
   	auto_nmx=atoi(s2);
	notAlreadySet.NMAX=0;
   }
   return;
 }
if(msc("NPR",s1)){
   if ((notAlreadySet.NPR||force)|| ((mask!=NULL)&&(mask->NPR==1)))
   {
   	auto_npr=atoi(s2);
	notAlreadySet.NPR=0;
   }
   return;
 }
 if(msc("NCOL",s1)){
   if ((notAlreadySet.NCOL||force)|| ((mask!=NULL)&&(mask->NCOL==1)))
   {
   	auto_ncol=atoi(s2);
   	notAlreadySet.NCOL=0;
   }
   return;
 }


if(msc("DSMIN",s1)){
   if ((notAlreadySet.DSMIN||force)|| ((mask!=NULL)&&(mask->DSMIN==1)))
   {
   	auto_dsmin=atof(s2);
	notAlreadySet.DSMIN=0;
   }
   return;
 }
if(msc("DSMAX",s1)){
   if ((notAlreadySet.DSMAX||force)|| ((mask!=NULL)&&(mask->DSMAX==1)))
   {
   	auto_dsmax=atof(s2);
   	notAlreadySet.DSMAX=0;
   }
   return;
 }
if(msc("DS",s1)){
    if ((notAlreadySet.DS||force)|| ((mask!=NULL)&&(mask->DS==1)))
    {
   	auto_ds=atof(s2);
	notAlreadySet.DS=0;
    }
 
   return;
 }
if(msc("PARMIN",s1)){
   if ((notAlreadySet.XMAX||force)|| ((mask!=NULL)&&(mask->XMAX==1)))
   {
   	auto_rl0=atof(s2);
	notAlreadySet.XMAX=0;
   }
   return;
 }
if(msc("PARMAX",s1)){
    if ((notAlreadySet.PARMAX||force)|| ((mask!=NULL)&&(mask->PARMAX==1)))
    {
   	auto_rl1=atof(s2);
	notAlreadySet.PARMAX=0;
    }
   return;
 }
if(msc("NORMMIN",s1)){
     if ((notAlreadySet.NORMMIN||force)|| ((mask!=NULL)&&(mask->NORMMIN==1)))
     {
   	auto_a0=atof(s2);
	notAlreadySet.NORMMIN=0;
     }
   return;
 }
if(msc("NORMMAX",s1)){
     if ((notAlreadySet.NORMMAX||force)|| ((mask!=NULL)&&(mask->NORMMAX==1)))
     {
   	auto_a1=atof(s2);
   	notAlreadySet.NORMMAX=0;
     }
   return;
 }
 if(msc("EPSL",s1)){
     if ((notAlreadySet.EPSL||force)|| ((mask!=NULL)&&(mask->EPSL==1)))
     {
   	auto_epsl=atof(s2);
	notAlreadySet.EPSL=0;
     }
   return;
 }

if(msc("EPSU",s1)){
     if ((notAlreadySet.EPSU||force)|| ((mask!=NULL)&&(mask->EPSU==1)))
     {
   	auto_epsu=atof(s2);
	notAlreadySet.EPSU=0;
     }
   return;
 }
if(msc("EPSS",s1)){
     if ((notAlreadySet.EPSS||force)|| ((mask!=NULL)&&(mask->EPSS==1)))
     {
   	auto_epss=atof(s2);
	notAlreadySet.EPSS=0;
     }
   return;
 }
 if(msc("RUNNOW",s1)){
     if ((notAlreadySet.RUNNOW||force)|| ((mask!=NULL)&&(mask->RUNNOW==1)))
     {
   	RunImmediately=atoi(s2);
	notAlreadySet.RUNNOW=0;
     }
   return;
 }

 if(msc("SEC",s1)){
     if ((notAlreadySet.SEC||force)|| ((mask!=NULL)&&(mask->SEC==1)))
     {
   	SEc=atoi(s2);
	notAlreadySet.SEC=0;
     }
   return;
 }
 if(msc("UEC",s1)){
     if ((notAlreadySet.UEC||force)|| ((mask!=NULL)&&(mask->UEC==1)))
     {
   	UEc=atoi(s2);
	notAlreadySet.UEC=0;
     }
   return;
 }
 if(msc("SPC",s1)){
     if ((notAlreadySet.SPC||force)|| ((mask!=NULL)&&(mask->SPC==1)))
     {
   	SPc=atoi(s2);
	notAlreadySet.SPC=0;
     }
   return;
 }
 if(msc("UPC",s1)){
     if ((notAlreadySet.UPC||force)|| ((mask!=NULL)&&(mask->UPC==1)))
     {
   	UPc=atoi(s2);
	notAlreadySet.UPC=0;
     }
   return;
 }

 if(msc("AUTOEVAL",s1)){
     if ((notAlreadySet.AUTOEVAL||force)|| ((mask!=NULL)&&(mask->AUTOEVAL==1)))
     {
   	f=atoi(s2);
   	set_auto_eval_flags(f);
	notAlreadySet.AUTOEVAL=0;
    }
   return;
 }
if(msc("AUTOXMAX",s1)){
     if ((notAlreadySet.AUTOXMAX||force)|| ((mask!=NULL)&&(mask->AUTOXMAX==1)))
     {
 	auto_xmax=atof(s2);
	notAlreadySet.AUTOXMAX=0;
     }
 return;
}
if(msc("AUTOYMAX",s1)){
     if ((notAlreadySet.AUTOYMAX||force)|| ((mask!=NULL)&&(mask->AUTOYMAX==1)))
     {
 		auto_ymax=atof(s2);
		notAlreadySet.AUTOYMAX=0;
     }
 return;
}
if(msc("AUTOXMIN",s1)){
     if ((notAlreadySet.AUTOXMIN||force)|| ((mask!=NULL)&&(mask->AUTOXMIN==1)))
     {
 	auto_xmin=atof(s2);
	notAlreadySet.AUTOXMIN=0;
     }
 return;
}
if(msc("AUTOYMIN",s1)){
     if ((notAlreadySet.AUTOYMIN||force)|| ((mask!=NULL)&&(mask->AUTOYMIN==1)))
     {
 	auto_ymin=atof(s2);
	notAlreadySet.AUTOYMIN=0;
     }
 return;
}
if(msc("AUTOVAR",s1)){
     if ((notAlreadySet.AUTOVAR||force)|| ((mask!=NULL)&&(mask->AUTOVAR==1)))
     {
     	find_variable(s2,&i);
    	if(i>0)auto_var=i-1;
	notAlreadySet.AUTOVAR=0;
    }
    return;
  }

/* postscript options */

 if(msc("PS_FONT",s1)){
     if ((notAlreadySet.PS_FONT||force)|| ((mask!=NULL)&&(mask->PS_FONT==1)))
     {
   	strcpy(PS_FONT,s2);
	notAlreadySet.PS_FONT=0;
     }
   return;
 }

if(msc("PS_LW",s1)){
   if ((notAlreadySet.PS_LW||force)|| ((mask!=NULL)&&(mask->PS_LW==1)))
   {
  	PS_LW=atof(s2);
	notAlreadySet.PS_LW=0;
   }
   return;
 }

if(msc("PS_FSIZE",s1)){
     if ((notAlreadySet.PS_FSIZE||force)|| ((mask!=NULL)&&(mask->PS_FSIZE==1)))
     {
  	PS_FONTSIZE=atoi(s2);
	notAlreadySet.PS_FSIZE=0;
     }
   return;
 }

if(msc("PS_COLOR",s1)){
     if ((notAlreadySet.PS_COLOR||force)|| ((mask!=NULL)&&(mask->PS_COLOR==1)))
     {
  	PSColorFlag=atoi(s2);
  	PS_Color=PSColorFlag;
	notAlreadySet.PS_COLOR=0;
     }
   return;
 }
if(msc("TUTORIAL",s1)){
   if(!(msc(s2,"0")||msc(s2,"1")))
   {
   	plintf("TUTORIAL option must be 0 or 1.\n");
	exit(-1);
   }
   if ((notAlreadySet.TUTORIAL||force) || ((mask!=NULL)&&(mask->TUTORIAL==1)))
   {
   	DoTutorial=atoi(s2);
	notAlreadySet.TUTORIAL=0;
   }
   return;
 }
 if(msc("S1",s1)){
     if ((notAlreadySet.SLIDER1||force) || ((mask!=NULL)&&(mask->SLIDER1==1)))
     {
	strncpy(SLIDER1VAR,s2,20);
	SLIDER1VAR[19]= '\0';
	notAlreadySet.SLIDER1=0;
     }
    return;
  }

if(msc("S2",s1)){
     if ((notAlreadySet.SLIDER2||force) || ((mask!=NULL)&&(mask->SLIDER2==1)))
     {
    	strncpy(SLIDER2VAR,s2,20);
	SLIDER2VAR[19]= '\0';
	notAlreadySet.SLIDER2=0;
     }
    return;
   }
 if(msc("S3",s1)){
     if ((notAlreadySet.SLIDER3||force) || ((mask!=NULL)&&(mask->SLIDER3==1)))
     {	
     	strncpy(SLIDER3VAR,s2,20);
	SLIDER3VAR[19]= '\0';
	notAlreadySet.SLIDER3=0;
     }
    return;
  }
  if(msc("SLO1",s1)){
     if ((notAlreadySet.SLIDER1LO||force) || ((mask!=NULL)&&(mask->SLIDER1LO==1)))
     {
    	SLIDER1LO=atof(s2);
	notAlreadySet.SLIDER1LO=0;
     }
    return;
  }

if(msc("SLO2",s1)){
     if ((notAlreadySet.SLIDER2LO||force) || ((mask!=NULL)&&(mask->SLIDER2LO==1)))
     {
    	SLIDER2LO=atof(s2);
	notAlreadySet.SLIDER2LO=0;
     }
    return;
   }
 if(msc("SLO3",s1)){
     if ((notAlreadySet.SLIDER3LO||force) || ((mask!=NULL)&&(mask->SLIDER3LO==1)))
     {
    	SLIDER3LO=atof(s2);
	notAlreadySet.SLIDER3LO=0;
     }
    return;
  }
 if(msc("SHI1",s1)){
     if ((notAlreadySet.SLIDER1HI||force) || ((mask!=NULL)&&(mask->SLIDER1HI==1)))
     {
    	SLIDER1HI=atof(s2);
	notAlreadySet.SLIDER1HI=0;
     }
    return;
  }
 if(msc("SHI2",s1)){
     if ((notAlreadySet.SLIDER2HI||force) || ((mask!=NULL)&&(mask->SLIDER2HI==1)))
     {
    	SLIDER2HI=atof(s2);
	notAlreadySet.SLIDER2HI=0;
     }
    return;
   }
 if(msc("SHI3",s1)){
     if ((notAlreadySet.SLIDER3HI||force) || ((mask!=NULL)&&(mask->SLIDER3HI==1)))
     {
    	SLIDER3HI=atof(s2);
	notAlreadySet.SLIDER3HI=0;
     }
    return;
 }

 /* postprocessing options
    This is rally only relevant for batch jobs as it 
    writes files then
 */

 if(msc("POSTPROCESS",s1)){
     if ((notAlreadySet.POSTPROCESS||force) || ((mask!=NULL)&&(mask->POSTPROCESS==1)))
     {
    	post_process=atoi(s2);
	notAlreadySet.POSTPROCESS=0;
     }
    return;
   }
   
 if(msc("HISTLO",s1)){
     if ((notAlreadySet.HISTLO||force) || ((mask!=NULL)&&(mask->HISTLO==1)))
     {
    	hist_inf.xlo=atof(s2);
	notAlreadySet.HISTLO=0;
     }
    return;
  }

 if(msc("HISTHI",s1)){
     if ((notAlreadySet.HISTHI||force) || ((mask!=NULL)&&(mask->HISTHI==1)))
     {
    	hist_inf.xhi=atof(s2);
	notAlreadySet.HISTHI=0;
     }
    return;
  }

 if(msc("HISTBINS",s1)){
     if ((notAlreadySet.HISTBINS||force) || ((mask!=NULL)&&(mask->HISTBINS==1)))
     {
    	hist_inf.nbins=atoi(s2);
	notAlreadySet.HISTBINS=0;
     }
    return;
  }

 if(msc("HISTCOL",s1)){
     if ((notAlreadySet.HISTCOL||force) || ((mask!=NULL)&&(mask->HISTCOL==1)))
     {
       find_variable(s2,&i);
       if(i>(-1)) hist_inf.col=i;
	notAlreadySet.HISTCOL=0;
     }
    return;
  }

 if(msc("SPECCOL",s1)){
     if ((notAlreadySet.SPECCOL||force) || ((mask!=NULL)&&(mask->SPECCOL==1)))
     {
       find_variable(s2,&i);
       if(i>(-1)) spec_col=i;
	notAlreadySet.SPECCOL=0;
     }
    return;
  }

 if(msc("SPECCOL2",s1)){
     if ((notAlreadySet.SPECCOL2||force) || ((mask!=NULL)&&(mask->SPECCOL2==1)))
     {
       find_variable(s2,&i);
       if(i>(-1)) spec_col2=i;
	notAlreadySet.SPECCOL2=0;
     }
    return;
  }

 if(msc("SPECWIDTH",s1)){
     if ((notAlreadySet.SPECWIDTH||force) || ((mask!=NULL)&&(mask->SPECWIDTH==1)))
     {
       spec_wid=atoi(s2);
	notAlreadySet.SPECWIDTH=0;
     }
    return;
  }

 if(msc("SPECWIN",s1)){
     if ((notAlreadySet.SPECWIN||force) || ((mask!=NULL)&&(mask->SPECWIN==1)))
     {
       spec_win=atoi(s2);
	notAlreadySet.SPECWIN=0;
     }
    return;
  }


  if(msc("DFGRID",s1)){
     if ((notAlreadySet.DFGRID||force)|| ((mask!=NULL)&&(mask->DFGRID==1)))
     { 
     	DF_GRID=atoi(s2);
	notAlreadySet.DFGRID=0;
     }
   return;
 }
  if(msc("DFDRAW",s1)){ 
     if ((notAlreadySet.DFBATCH||force)|| ((mask!=NULL)&&(mask->DFBATCH==1)))
     { 
     	DFBatch=atoi(s2);
	notAlreadySet.DFBATCH=0;
     }
   return;
 }
   if(msc("NCDRAW",s1)){
     if ((notAlreadySet.NCBATCH||force)|| ((mask!=NULL)&&(mask->NCBATCH==1)))
     { 
     	NCBatch=atoi(s2);
	notAlreadySet.NCBATCH=0;
     }
   return;
   }

   /* colorize customizing !! */
   if(msc("COLORVIA",s1))
     {
       if ((notAlreadySet.COLORVIA||force)|| ((mask!=NULL)&&(mask->COLORVIA==1)))
       strcpy(ColorVia,s2);
       	notAlreadySet.COLORVIA=0;
       return;
     }
   if(msc("COLORIZE",s1))
     {
          if ((notAlreadySet.COLORIZE||force)|| ((mask!=NULL)&&(mask->COLORIZE==1)))
       ColorizeFlag=atoi(s2);
          	notAlreadySet.COLORIZE=0;
          return;
     }
   if(msc("COLORLO",s1))
     {
              if ((notAlreadySet.COLORLO||force)|| ((mask!=NULL)&&(mask->COLORLO==1)))
       ColorViaLo=atof(s2);
	             	notAlreadySet.COLORLO=0;
          return;
     }
   if(msc("COLORHI",s1))
     {
              if ((notAlreadySet.COLORHI||force)|| ((mask!=NULL)&&(mask->COLORHI==1)))
       ColorViaHi=atof(s2);
              	notAlreadySet.COLORHI=0;
       return;
     }

plintf("!! Option %s not recognized\n",s1); 
  
}





