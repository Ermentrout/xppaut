#include "load_eqn.h"

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "adj2.h"
#include "arrayplot.h"
#include "auto_nox.h"
#include "browse.h"
#include "color.h"
#include "dae_fun.h"
#include "delay_handle.h"
#include "derived.h"
#include "extra.h"
#include "solver/gear.h"
#include "ggets.h"
#include "graf_par.h"
#include "graphics.h"
#include "histogram.h"
#include "init_conds.h"
#include "integrate.h"
#include "lunch-new.h"
#include "main.h"
#include "many_pops.h"
#include "markov.h"
#include "my_ps.h"
#include "nullcline.h"
#include "numerics.h"
#include "parserslow.h"
#include "read_dir.h"
#include "storage.h"
#include "util/strutil.h"
#include "tabular.h"
#include "txtread.h"
#include "userbut.h"
#include "solver/runge_kutta.h"
#include "solver/volterra2.h"
#include "xpplim.h"

/* --- Macros --- */
#define DFNORMAL 1
#define MAXOPT 1000


/* --- Forward declarations --- */
static void check_volterra(void);
static void do_intern_set(char *name1, char *value);
static void init_X_vals(void);
static void split_apart(char *bob, char *name, char *value);
int (*solver)();


/* --- Data --- */
static float oldhp_x,oldhp_y,my_pl_wid,my_pl_ht;

char delay_string[MAXODE][80];
char this_file[XPP_MAX_NAME];
char this_internset[XPP_MAX_NAME];
double last_ic[MAXODE];
double MY_XLO,MY_YLO,MY_XHI,MY_YHI;
double TOR_PERIOD=6.2831853071795864770;
double x_3d[2],y_3d[2],z_3d[2];
double X_LO[10],Y_LO[10],X_HI[10],Y_HI[10];

int AXES,TIMPLOT,PLOT_3D;
int itor[MAXODE];
int IX_PLT[10],IY_PLT[10],IZ_PLT[10],NPltV;
int IXPLT,IYPLT,IZPLT;
int mov_ind;
int MultiWin=0;
int NEQ;
int Nintern_set=0;
int RunImmediately=0;
int START_LINE_TYPE=1;
int  storind,STORFLAG,INFLAG,MAXSTOR;
int TORUS=0;
INTERN_SET intern_set[MAX_INTERN_SET];


/*   Numerical stuff ....   */
double DELTA_T,TEND,T0,TRANS,NULL_ERR,EVEC_ERR,NEWT_ERR;
double BOUND,DELAY,TOLER,ATOLER,HMIN,HMAX;
double BVP_EPS,BVP_TOL;
double POIPLN;
int MaxEulIter;
double EulTol;
int NMESH,NJMP,color_flag,NC_ITER;
int EVEC_ITER;
int BVP_MAXIT,BVP_FLAG;
int POIMAP,POIVAR,POISGN,SOS;
int FFT,NULL_HERE,POIEXT;
int HIST,HVAR,hist_ind,FOREVER;
/*  control of range stuff  */
int PAUSER,ENDSING,SHOOT,PAR_FOL;
int xorfix,silent;


/* --- Functions --- */
/* here is some new code for internal set files:
 * format of the file is a long string of the form:
 * { x=y, z=w, q=p , .... }
 */
void add_intern_set(char *name, char *does) {
	char bob[DEFAULT_STRING_LENGTH],ch;
	int i,n,j=Nintern_set,k=0;
	if(Nintern_set>=MAX_INTERN_SET) {
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
	for(i=0;i<n;i++) {
		ch=does[i];
		if(ch==',') {
			bob[k]=' ';
			k++;
		}
		if(ch=='}' || ch=='{') {
			continue;
		}
		if(ch!=',') {
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


/*
 * Load an ~/.xpprc if it exists.
 *
 * The RC file has the following syntax:
 *
 *   file   <- line*
 *   line   <- '@' <SP> assign ((<SP>|',') assign)*
 *   assign <- char+ '=' char+
 *
 * Lines not starting with '@' are ignored.
 */
void loadeqn_load_xpprc(void) {
	char rc[DEFAULT_STRING_LENGTH];

	sprintf(rc,"%s/.xpprc",getenv("HOME"));
	FILE *fp=fopen(rc, "r");
	if(fp==NULL) {
		return;
	}
	while(!feof(fp)) {
		char bob[DEFAULT_STRING_LENGTH];
		bob[0]='\0';
		if(fgets(bob, sizeof(bob) - 1, fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
		if (bob[0] != '@') {
			continue;
		}
		loadeqn_set_internopt(bob);
	}
	fclose(fp);
}


void extract_action(char *ptr) {
	char name[XPP_MAX_NAME],value[XPP_MAX_NAME];
	char tmp[2048];
	char *mystring, *toksave;
	strcpy(tmp,ptr);
	strtok_r(tmp," ", &toksave);

	while((mystring=strtok_r(NULL, " ,;\n", &toksave))!=NULL) {
		split_apart(mystring,name,value);
		if(strlen(name)>0 && strlen(value)>0) {
			do_intern_set(name,value);
		}
	}
}


void extract_internset(int j) {
	extract_action(intern_set[j].does);
}


void load_eqn(void) {
	int no_eqn=1,okay=0;
	int i;
	int std=0;
	FILE *fptr;
	init_ar_ic();
	for(i=0;i<MAXODE;i++) {
		itor[i]=0;
		strcpy(delay_string[i],"0.0");
	}
	if(strcmp(this_file,"/dev/stdin")==0) {
		std=1;
	}
	struct dirent *dp;
	if(this_file[0] && (std==0)  &&
	   (dp=(struct dirent*)opendir(this_file))!=NULL) {
		no_eqn = 1;
		okay=0;
		change_directory(this_file);
		okay=make_eqn();
		return;
	} else {
		if(this_file[0] && (fptr=fopen(this_file,"r"))!=NULL) {
			if(std==1) {
				sprintf(this_file,"console");
			}
			okay=get_eqn(fptr);
			if(std==0) {
				fclose(fptr);
			}
			if(okay==1) {
				no_eqn=0;
			}
		}
	}
	if(no_eqn) {
		while(okay==0) {
			struct dirent *dp;
			char odeclassrm[XPP_MAX_NAME];
			if(getenv("XPPSTART")!=NULL) {
				sprintf(odeclassrm,"%s",getenv("XPPSTART"));
				if((dp=(struct dirent*)opendir(odeclassrm))!=NULL)	{
					change_directory(odeclassrm);
				}
			}
			okay=make_eqn();
		}
	}
}


void loadeqn_init_options(void) {
	TIMPLOT=1;
	FOREVER=0;
	BVP_TOL=1.e-5;
	BVP_EPS=1.e-5;
	BVP_MAXIT=20;
	BVP_FLAG=0;
	NMESH=40;
	NJMP=1;
	SOS=0;
	FFT=0;
	HIST=0;
	PltFmtFlag=0;
	AXES=0;
	TOLER=0.001;
	ATOLER=0.001;
	MaxEulIter=10;
	EulTol=1.e-7;
	DELAY=0.0;
	HMIN=1e-12;
	EVEC_ITER=100;
	EVEC_ERR=.001;
	NULL_ERR=.001;
	NEWT_ERR=.001;
	NULL_HERE=0;
	del_stab_flag=DFNORMAL;
	HMAX=1.000;
	POIMAP=0;
	POIVAR=1;
	POIEXT=0;
	POISGN=1;
	POIPLN=0.0;
	storind=0;
	mov_ind=0;
	STORFLAG=0;
	INFLAG=0;

	oldhp_x=-100000.0 ;
	oldhp_y=-100000.0;
	solver=runge_kutta;
	PLOT_3D=0;

	METHOD=METHOD_RK4;
	MY_XLO=0.0;
	MY_XHI=20.0;
	MY_YLO=-1;
	MY_YHI=1;
	x_3d[0]=-12;
	x_3d[1]=12;
	y_3d[0]=-12;
	y_3d[1]=12;
	z_3d[0]=-12;
	z_3d[1]=12;

	BOUND=100;
	MAXSTOR=5000;
	my_pl_wid=10000. ;
	my_pl_ht=7000.  ;

	/* TORUS=0; */
	T0=0.0;
	TRANS=0.0;
	DELTA_T=.05;
	TEND=20.00;
	TOR_PERIOD=6.2831853071795864770;
	IXPLT=0;
	IYPLT=1;
	IZPLT=1;
	if(NEQ>2) {
		IZPLT=2;
	}
	NPltV=1;
	for(int i=0;i<10;i++) {
		IX_PLT[i]=IXPLT;
		IY_PLT[i]=IYPLT;
		IZ_PLT[i]=IZPLT;
		X_LO[i]=0;
		Y_LO[i]=-1;
		X_HI[i]=20;
		Y_HI[i]=1;
	}

	integrate_init_range();
	adj2_init_trans();
	init_my_aplot();
	init_txtview();
	init_X_vals();
}


void loadeqn_set_internopt(char *line) {
	char name[20], value[80], *mystring;
	char *toksave;

	strtok_r(line, " ,", &toksave);
	while ((mystring = strtok_r(NULL, " ,\n\r", &toksave)) != NULL) {
		split_apart(mystring, name, value);
		if (strlen(name) > 0 && strlen(value) > 0) {
			set_option(name, value);
		}
	}
}


void loadeqn_setup_all(void) {
	if(IZPLT>NEQ) {
		IZPLT=NEQ;
	}
	if(IYPLT>NEQ) {
		IYPLT=NEQ;
	}
	if(IXPLT==0 || IYPLT==0) {
		TIMPLOT=1;
	} else {
		TIMPLOT=0;
	}
	if(x_3d[0]>=x_3d[1]) {
		x_3d[0]=-1;
		x_3d[1]=1;
	}
	if(y_3d[0]>=y_3d[1]) {
		y_3d[0]=-1;
		y_3d[1]=1;
	}
	if(z_3d[0]>=z_3d[1]) {
		z_3d[0]=-1;
		z_3d[1]=1;
	}
	if(MY_XLO>=MY_XHI) {
		MY_XLO=-2.0;
		MY_XHI=2.0;
	}
	if(MY_YLO>=MY_YHI) {
		MY_YLO=-2.0;
		MY_YHI=2.0;
	}
	if(AXES<5) {
		x_3d[0]=MY_XLO;
		y_3d[0]=MY_YLO;
		x_3d[1]=MY_XHI;
		y_3d[1]=MY_YHI;
	}
	init_stor(MAXSTOR,NEQ+1);
	if(AXES>=5) {
		PLOT_3D=1;
	}
	check_volterra();

	integrate_setup_range();
	adj2_setup_trans();
	init_stor(MAXSTOR, NEQ + 1);
	check_delay(); /* check for delay allocation */
	alloc_h_stuff();
	alloc_v_memory();  /* allocate stuff for volterra equations */
	alloc_meth();
	arr_ic_start(); /* take care of all predefined array ics */
}


void set_option(char *s1, char *s2) {
	int i,j;
	char xx[4],yy[4],zz[4];
	char xxl[6],xxh[6],yyl[6],yyh[6];
	/* Must match enum Method. */
	static char mkey[]="demragvbqsc582y";
	static char Mkey[]="DEMRAGVBQSC582Y";
	strupr(s1);
	if(strprefix("QUIET",s1)) {
		if(!(strprefix(s2,"0") || strprefix(s2,"1")))	{
			plintf("QUIET option must be 0 or 1.\n");
			exit(-1);
		}
		XPPVERBOSE=(atoi(s2)==0);
		return;
	}
	if(strprefix("LOGFILE",s1)) {
		/*Will be 1 if -logfile was specified on the command line.*/
		if(logfile != NULL) {
			fclose(logfile);
		}
		logfile=fopen(s2,"w");
		return;
	}
	if(strprefix("BELL",s1)) {
		if(!(strprefix(s2,"0") || strprefix(s2,"1")))	{
			plintf("BELL option must be 0 or 1.\n");
			exit(-1);
		}
		tfBell=atoi(s2);
		return;
	}
	if(strprefix("BUT",s1)) {
		add_user_button(s2);
		return;
	}
	if((strprefix("BIGFONT",s1)) || (strprefix("BIG",s1))) {
		strcpy(big_font_name,s2);
		return;
	}
	if((strprefix("SMALLFONT",s1)) || (strprefix("SMALL",s1))) {
		strcpy(small_font_name,s2);
		return;
	}
	if(strprefix("FORECOLOR",s1)) {
		sprintf(UserBlack,"#%s",s2);
		return;
	}
	if(strprefix("BACKCOLOR",s1)) {
		sprintf(UserWhite,"#%s",s2);
		return;
	}
	if(strprefix("MWCOLOR",s1)) {
		sprintf(UserMainWinColor,"#%s",s2);
		return;
	}
	if(strprefix("DWCOLOR",s1)) {
		sprintf(UserDrawWinColor,"#%s",s2);
		return;
	}
	if(strprefix("GRADS",s1)) {
		if(!(strprefix(s2,"0") || strprefix(s2,"1"))) {
			plintf("GRADS option must be 0 or 1.\n");
			exit(-1);
		}
		UserGradients=atoi(s2);
		return;
	}
	if(strprefix("PLOTFMT",s1)) {
		strcpy(PlotFormat,s2);
		return;
	}
	if(strprefix("BACKIMAGE",s1)) {
		strcpy(UserBGBitmap,s2);
		return;
	}
	if(strprefix("WIDTH",s1)) {
		UserMinWidth=atoi(s2);
		return;
	}
	if(strprefix("HEIGHT",s1)) {
		UserMinHeight=atoi(s2);
		return;
	}
	if(strprefix("YNC",s1)) {
		i=atoi(s2);
		if(i>-1 && i<11) {
			YNullColor=i;
		}
		return;
	}
	if(strprefix("XNC",s1)) {
		i=atoi(s2);
		if(i>-1 && i<11) {
			XNullColor=i;
		}
		return;
	}
	if(strprefix("SMC",s1)) {
		i=atoi(s2);
		if(i>-1 && i<11) {
			StableManifoldColor=i;
		}
		return;
	}
	if(strprefix("UMC",s1)) {
		i=atoi(s2);
		if(i>-1 && i<11) {
			UnstableManifoldColor=i;
		}
		return;
	}
	if(strprefix("LT",s1)) {
		i=atoi(s2);
		if(i<2 && i>-6) {
			START_LINE_TYPE=i;
			reset_all_line_type();
		}
		return;
	}
	if(strprefix("SEED",s1)) {
		i=atoi(s2);
		if(i>=0) {
			RandSeed=i;
			nsrand48(RandSeed);
		}
		return;
	}
	if(strprefix("BACK",s1)) {
		if(s2[0]=='w' || s2[0]=='W') {
			PaperWhite=1;
		} else {
			PaperWhite=0;
		}
		return;
	}
	if(strprefix("COLORMAP",s1)) {
		i=atoi(s2);
		if(i<7) {
			custom_color=i;
		}
		return;
	}
	if(strprefix("NPLOT",s1)) {
		NPltV=atoi(s2);
		return;
	}
	if(strprefix("DLL_LIB",s1)) {
		sprintf(dll_lib,"%s",s2);
		dll_flag+=1;
		return;
	}
	if(strprefix("DLL_FUN",s1)) {
		sprintf(dll_fun,"%s",s2);
		dll_flag+=2;
		return;
	}
	/* can now initialize several plots */
	if(strprefix("SIMPLOT",s1)) {
		SimulPlotFlag=1;
		return;
	}
	if(strprefix("MULTIWIN",s1)) {
		MultiWin=1;
		return;
	}
	for(j=2;j<=8;j++) {
		sprintf(xx,"XP%d",j);
		sprintf(yy,"YP%d",j);
		sprintf(zz,"ZP%d",j);
		sprintf(xxh,"XHI%d",j);
		sprintf(xxl,"XLO%d",j);
		sprintf(yyh,"YHI%d",j);
		sprintf(yyl,"YLO%d",j);
		if(strprefix(xx,s1)) {
			find_variable(s2,&i);
			if(i>-1) {
				IX_PLT[j]=i;
			}
			return;
		}
		if(strprefix(yy,s1)) {
			find_variable(s2,&i);
			if(i>-1) {
				IY_PLT[j]=i;
			}
			return;
		}
		if(strprefix(zz,s1)) {
			find_variable(s2,&i);
			if(i>-1) {
				IZ_PLT[j]=i;
			}
			return;
		}
		if(strprefix(xxh,s1)) {
			X_HI[j]=atof(s2);
			return;
		}
		if(strprefix(xxl,s1)) {
			X_LO[j]=atof(s2);
			return;
		}
		if(strprefix(yyh,s1)) {
			Y_HI[j]=atof(s2);
			return;
		}
		if(strprefix(yyl,s1)) {
			Y_LO[j]=atof(s2);
			return;
		}
	}
	if(strprefix("XP",s1)) {
		find_variable(s2,&i);
		if(i>-1) {
			IXPLT=i;
		}
		return;
	}
	if(strprefix("YP",s1)) {
		find_variable(s2,&i);
		if(i>-1) {
			IYPLT=i;
		}
		return;
	}
	if(strprefix("ZP",s1)) {
		find_variable(s2,&i);
		if(i>-1) {
			IZPLT=i;
		}
		return;
	}
	if(strprefix("AXES",s1)) {
		if(s2[0]=='3') {
			AXES=5;
		} else {
			AXES=0;
		}
		return;
	}
	if(strprefix("NJMP",s1)) {
		NJMP=atoi(s2);
		return;
	}
	if(strprefix("NOUT",s1)) {
		NJMP=atoi(s2);
		return;
	}
	if(strprefix("NMESH",s1)) {
		NMESH=atoi(s2);
		return;
	}
	if(strprefix("METH",s1)) {
		for(i=0;i<15;i++) {
			if(s2[0]==mkey[i] || s2[0]==Mkey[i]) {
				METHOD=i;
			}
		}
		return;
	}
	if(strprefix("VMAXPTS",s1)) {
		MaxPoints=atoi(s2);
		return;
	}
	if(strprefix("MAXSTOR",s1)) {
		MAXSTOR=atoi(s2);
		return;
	}
	if(strprefix("TOR_PER",s1)) {
		TOR_PERIOD=atof(s2);
		TORUS=1;
		return;
	}
	if(strprefix("JAC_EPS",s1)) {
		NEWT_ERR=atof(s2);
		return;
	}
	if(strprefix("NEWT_TOL",s1)) {
		EVEC_ERR=atof(s2);
		return;
	}
	if(strprefix("NEWT_ITER",s1)) {
		EVEC_ITER=atoi(s2);
		return;
	}
	if(strprefix("FOLD",s1)) {
		find_variable(s2,&i);
		if(i>0) {
			itor[i-1]=1;
			TORUS=1;
		}
		return;
	}
	if(strprefix("TOTAL",s1)) {
		TEND=atof(s2);
		return;
	}
	if(strprefix("DTMIN",s1)) {
		HMIN=atof(s2);
		return;
	}
	if(strprefix("DTMAX",s1)) {
		HMAX=atof(s2);
		return;
	}
	if(strprefix("DT",s1)) {
		DELTA_T=atof(s2);
		return;
	}
	if(strprefix("T0",s1)) {
		T0=atof(s2);
		return;
	}
	if(strprefix("TRANS",s1)) {
		TRANS=atof(s2);
		return;
	}
	if(strprefix("BOUND",s1)) {
		BOUND=atof(s2);
		return;
	}
	if(strprefix("ATOL",s1)) {
		ATOLER=atof(s2);
		return;
	}
	if(strprefix("TOL",s1)) {
		TOLER=atof(s2);
		return;
	}
	if(strprefix("DELAY",s1)) {
		DELAY=atof(s2);
		return;
	}
	if(strprefix("BANDUP",s1)) {
		cv_bandflag=1;
		cv_bandupper=atoi(s2);
		return;
	}
	if(strprefix("BANDLO",s1)) {
		cv_bandflag=1;
		cv_bandlower=atoi(s2);
		return;
	}
	if(strprefix("PHI",s1)) {
		PHI0=atof(s2);
		return;
	}
	if(strprefix("THETA",s1)) {
		THETA0=atof(s2);
		return;
	}
	if(strprefix("XLO",s1)) {
		MY_XLO=atof(s2);
		return;
	}
	if(strprefix("YLO",s1)) {
		MY_YLO=atof(s2);
		return;
	}
	if(strprefix("XHI",s1)) {
		MY_XHI=atof(s2);
		return;
	}
	if(strprefix("YHI",s1)) {
		MY_YHI=atof(s2);
		return;
	}
	if(strprefix("XMAX",s1)) {
		x_3d[1]=atof(s2);
		return;
	}
	if(strprefix("YMAX",s1)) {
		y_3d[1]=atof(s2);
		return;
	}
	if(strprefix("ZMAX",s1)) {
		z_3d[1]=atof(s2);
		return;
	}
	if(strprefix("XMIN",s1)) {
		x_3d[0]=atof(s2);
		MY_XLO=atof(s2);
		return;
	}
	if(strprefix("YMIN",s1)) {
		y_3d[0]=atof(s2);
		MY_YLO=atof(s2);
		return;
	}
	if(strprefix("ZMIN",s1)) {
		z_3d[0]=atof(s2);
		return;
	}
	if(strprefix("POIMAP",s1)) {
		if(s2[0]=='m' || s2[0]=='M') {
			POIMAP=2;
		}
		if(s2[0]=='s' || s2[0]=='S') {
			POIMAP=1;
		}
		if(s2[0]=='p' || s2[0]=='P') {
			POIMAP=3;
		}
		return;
	}
	if(strprefix("POIVAR",s1)) {
		find_variable(s2,&i);
		if(i>-1) {
			POIVAR=i;
		}
		return;
	}
	if(strprefix("OUTPUT",s1)) {
		strcpy(batchout,s2);
		return;
	}
	if(strprefix("POISGN",s1)) {
		POISGN=atoi(s2);
		return;
	}
	if(strprefix("POISTOP",s1)) {
		SOS=atoi(s2);
		return;
	}
	if(strprefix("STOCH",s1)) {
		STOCH_FLAG=atoi(s2);
		return;
	}
	if(strprefix("POIPLN",s1)) {
		POIPLN=atof(s2);
		return;
	}
	if(strprefix("RANGEOVER",s1)) {
		strcpy(range.item,s2);
		return;
	}
	if(strprefix("RANGESTEP",s1)) {
		range.steps=atoi(s2);
		return;
	}
	if(strprefix("RANGELOW",s1)) {
		range.plow=atof(s2);
		return;
	}
	if(strprefix("RANGEHIGH",s1)) {
		range.phigh=atof(s2);
		return;
	}
	if(strprefix("RANGERESET",s1)) {
		if(s2[0]=='y' || s2[0]=='Y') {
			range.reset=1;
		} else {
			range.reset=0;
		}
		return;
	}
	if(strprefix("RANGEOLDIC",s1)) {
		if(s2[0]=='y' || s2[0]=='Y') {
			range.oldic=1;
		} else {
			range.oldic=0;
		}
		return;
	}
	if(strprefix("RANGE",s1)) {
		batch_range=atoi(s2);
		return;
	}
	if(strprefix("NTST",s1)) {
		auto_ntst=atoi(s2);
		return;
	}
	if(strprefix("NMAX",s1)) {
		auto_nmx=atoi(s2);
		return;
	}
	if(strprefix("NPR",s1)) {
		auto_npr=atoi(s2);
		return;
	}
	if(strprefix("NCOL",s1)) {
		auto_ncol=atoi(s2);
		return;
	}
	if(strprefix("DSMIN",s1)) {
		auto_dsmin=atof(s2);
		return;
	}
	if(strprefix("DSMAX",s1)) {
		auto_dsmax=atof(s2);
		return;
	}
	if(strprefix("DS",s1)) {
		auto_ds=atof(s2);
		return;
	}
	if(strprefix("PARMIN",s1)) {
		auto_rl0=atof(s2);
		return;
	}
	if(strprefix("PARMAX",s1)) {
		auto_rl1=atof(s2);
		return;
	}
	if(strprefix("NORMMIN",s1)) {
		auto_a0=atof(s2);
		return;
	}
	if(strprefix("NORMMAX",s1)) {
		auto_a1=atof(s2);
		return;
	}
	if(strprefix("EPSL",s1)) {
		auto_epsl=atof(s2);
		return;
	}
	if(strprefix("EPSU",s1)) {
		auto_epsu=atof(s2);
		return;
	}
	if(strprefix("EPSS",s1)) {
		auto_epss=atof(s2);
		return;
	}
	if(strprefix("RUNNOW",s1)) {
		RunImmediately=atoi(s2);
		return;
	}
	if(strprefix("SEC",s1)) {
		SEc=atoi(s2);
		return;
	}
	if(strprefix("UEC",s1)) {
		UEc=atoi(s2);
		return;
	}
	if(strprefix("SPC",s1)) {
		SPc=atoi(s2);
		return;
	}
	if(strprefix("UPC",s1)) {
		UPc=atoi(s2);
		return;
	}
	if(strprefix("AUTOEVAL",s1)) {
		set_auto_eval_flags(atoi(s2));
		return;
	}
	if(strprefix("AUTOXMAX",s1)) {
		auto_xmax=atof(s2);
		return;
	}
	if(strprefix("AUTOYMAX",s1)) {
		auto_ymax=atof(s2);
		return;
	}
	if(strprefix("AUTOXMIN",s1)) {
		auto_xmin=atof(s2);
		return;
	}
	if(strprefix("AUTOYMIN",s1)) {
		auto_ymin=atof(s2);
		return;
	}
	if(strprefix("AUTOVAR",s1)) {
		find_variable(s2,&i);
		if(i>0) {
			auto_var=i-1;
		}
		return;
	}
	/* postscript options */
	if(strprefix("PS_FONT",s1)) {
		strcpy(PS_FONT,s2);
		return;
	}
	if(strprefix("PS_LineWidth",s1)) {
		PS_LineWidth=atof(s2);
		return;
	}
	if(strprefix("PS_FSIZE",s1)) {
		PS_FontSize=atoi(s2);
		return;
	}
	if(strprefix("PS_COLOR",s1)) {
		PS_ColorFlag=atoi(s2);
		PS_Color=PS_ColorFlag;
		return;
	}
	if(strprefix("TUTORIAL",s1)) {
		if(!(strprefix(s2,"0") || strprefix(s2,"1")))	{
			plintf("TUTORIAL option must be 0 or 1.\n");
			exit(-1);
		}
		DoTutorial=atoi(s2);
		return;
	}
	if(strprefix("S1",s1)) {
		strncpy(SLIDERVAR[0],s2,20);
		SLIDERVAR[0][sizeof(SLIDERVAR[0]) - 1] = '\0';
		return;
	}
	if(strprefix("S2",s1)) {
		strncpy(SLIDERVAR[1],s2,20);
		SLIDERVAR[1][sizeof(SLIDERVAR[0]) - 1] = '\0';
		return;
	}
	if(strprefix("S3",s1)) {
		strncpy(SLIDERVAR[2],s2,20);
		SLIDERVAR[2][sizeof(SLIDERVAR[0]) - 1] = '\0';
		return;
	}
	if(strprefix("SLO1",s1)) {
		SLIDERLO[0]=atof(s2);
		return;
	}
	if(strprefix("SLO2",s1)) {
		SLIDERLO[1]=atof(s2);
		return;
	}
	if(strprefix("SLO3",s1)) {
		SLIDERLO[2]=atof(s2);
		return;
	}
	if(strprefix("SHI1",s1)) {
		SLIDERHI[0]=atof(s2);
		return;
	}
	if(strprefix("SHI2",s1)) {
		SLIDERHI[0]=atof(s2);
		return;
	}
	if(strprefix("SHI3",s1)) {
		SLIDERHI[0]=atof(s2);
		return;
	}
	/* postprocessing options
	 * This is really only relevant for batch jobs as it
	 * writes files then
	 */
	if(strprefix("POSTPROCESS",s1)) {
		post_process=atoi(s2);
		return;
	}
	if(strprefix("HISTLO",s1)) {
		hist_inf.xlo=atof(s2);
		return;
	}
	if(strprefix("HISTHI",s1)) {
		hist_inf.xhi=atof(s2);
		return;
	}
	if(strprefix("HISTBINS",s1)) {
		hist_inf.nbins=atoi(s2);
		return;
	}
	if(strprefix("HISTCOL",s1)) {
		find_variable(s2,&i);
		if(i>(-1)) {
			hist_inf.col=i;
		}
		return;
	}
	if(strprefix("SPECCOL",s1)) {
		find_variable(s2,&i);
		if(i>(-1)) {
			spec_col=i;
		}
		return;
	}
	if(strprefix("SPECCOL2",s1)) {
		find_variable(s2,&i);
		if(i>(-1)) {
			spec_col2=i;
		}
		return;
	}
	if(strprefix("SPECWIDTH",s1)) {
		spec_wid=atoi(s2);
		return;
	}
	if(strprefix("SPECWIN",s1)) {
		spec_win=atoi(s2);
		return;
	}
	if(strprefix("DFGRID",s1)) {
		DF_GRID=atoi(s2);
		return;
	}
	if(strprefix("DFDRAW",s1)) {
		DFBatch=atoi(s2);
		return;
	}
	if(strprefix("NCDRAW",s1)) {
		NCBatch=atoi(s2);
		return;
	}
	/* colorize customizing !! */
	if(strprefix("COLORVIA",s1)) {
		strcpy(ColorVia,s2);
		return;
	}
	if(strprefix("COLORIZE",s1)) {
		ColorizeFlag=atoi(s2);
		return;
	}
	if(strprefix("COLORLO",s1)) {
		ColorViaLo=atof(s2);
		return;
	}
	if(strprefix("COLORHI",s1)) {
		ColorViaHi=atof(s2);
		return;
	}
	plintf("!! Option %s not recognized\n",s1);
}


static void init_X_vals(void) {
	tfBell=1;
	strcpy(big_font_name,"fixed");
	strcpy(small_font_name,"6x13");
	sprintf(UserBlack,"#%s","000000");
	sprintf(UserWhite,"#%s","EDE9E3");
	sprintf(UserMainWinColor,"#%s","808080");
	sprintf(UserDrawWinColor,"#%s","FFFFFF");
	/*
	 * No gradients tends to look cleaner but some
	 * may prefer gradients improved contrast/readability.
	 */
	UserGradients=1;
}


/* --- Static functions --- */
static void check_volterra(void) {
	if (NKernel>0) {
		METHOD=METHOD_VOLTERRA;
	}
}


static void do_intern_set(char *name1, char *value) {
	int i;
	char name[20];
	convert(name1,name);

	i=find_user_name(INITIAL_CONDITION_BOX,name);
	if(i>-1) {
		last_ic[i]=atof(value);
	} else {
		i=find_user_name(PARAMETER_BOX,name);
		if(i>-1) {
			set_val(name,atof(value));
		} else {
			set_option(name,value);
		}
	}
	alloc_meth();
	do_meth();
}


static void split_apart(char *bob, char *name, char *value) {
	int k,i,l;
	l=strlen(bob);
	k=strcspn(bob,"=");
	if(k==l) {
		value[0]=0;
		strcpy(name,bob);
	} else {
		strncpy(name,bob,k);
		name[k]='\0';
		for(i=k+1;i<l;i++) {
			value[i-k-1]=bob[i];
		}
		value[l-k-1]='\0';
	}
}
