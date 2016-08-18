#include "lunch-new.h"

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "adj2.h"
#include "arrayplot.h"
#include "browse.h"
#include "delay_handle.h"
#include "edit_rhs.h"
#include "form_ode.h"
#include "ggets.h"
#include "graf_par.h"
#include "init_conds.h"
#include "integrate.h"
#include "load_eqn.h"
#include "main.h"
#include "many_pops.h"
#include "markov.h"
#include "numerics.h"
#include "parserslow.h"
#include "pp_shoot.h"
#include "storage.h"
#include "struct.h"
#include "solver/volterra2.h"
#include "xpplim.h"

/* --- Forward declarations --- */
static void do_info(FILE *fp);
static void dump_aplot(FILE *fp, int f);
static void dump_eqn(FILE *fp);
static void dump_h_stuff(FILE *fp, int f);
static void dump_range(FILE *fp, int f);
static void dump_shoot_range(FILE *fp, int f);
static void dump_torus(FILE *fp, int f);
static void dump_transpose_info(FILE *fp, int f);
static void io_double(double *z, FILE *fp, int f, char *ss);
static void io_exprs(int f, FILE *fp);
static void io_int(int *i, FILE *fp, int f, char *ss);
static void io_graph(int f, FILE *fp);
static void io_numerics(int f, FILE *fp);
static void io_parameters(int f, FILE *fp);
static void io_string(char *s, int len, FILE *fp, int f);

/* --- Data --- */
static int set_type=0;

/* --- Functions --- */
void do_lunch(int f) {
	int ne,np,ok,temp;
	char bob[DEFAULT_STRING_LENGTH];
	FILE *fp;
	time_t ttt;
	char filename[XPP_MAX_NAME];
	sprintf(filename,"%s.set",this_file);
	if(f==READEM) {
		ping();
		if(!file_selector("Load SET File",filename,"*.set")) {
			return;
		}
		fp=fopen(filename,"r");
		if(fp==NULL) {
			err_msg("Cannot open file");
			return;
		}
		if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
		if(bob[0]=='#') {
			set_type=1;
			io_int(&ne,fp,f," ");
		} else {
			ne=atoi(bob);
			set_type=0;
		}
		/* io_int(&ne,fp,f); */
		io_int(&np,fp,f," ");
		if(ne!=NEQ || np!=NUPAR) {
			err_msg("Incompatible parameters");
			fclose(fp);
			return;
		}
		io_numerics(f,fp);
		if(METHOD==METHOD_VOLTERRA) {
			io_int(&temp,fp,f," ");
			allocate_volterra(temp,1);
			MyStart=1;
		}
		check_delay();
		io_exprs(f,fp);
		io_graph(f,fp);
		if(set_type==1) {
			dump_transpose_info(fp,f);
			dump_h_stuff(fp,f);
			dump_aplot(fp,f);
			dump_torus(fp,f);
			dump_range(fp,f);
		}
		fclose(fp);
		return;
	}
	if(!file_selector("Save SET File",filename,"*.set")) {
		return;
	}
	open_write_file(&fp,filename,&ok);
	if(!ok) {
		return;
	}
	redraw_params();
	ttt=time(0);
	fprintf(fp,"## Set file for %s on %s",this_file,ctime(&ttt));
	io_int(&NEQ,fp,f,"Number of equations and auxiliaries");
	io_int(&NUPAR,fp,f,"Number of parameters");
	io_numerics(f,fp);
	if(METHOD==METHOD_VOLTERRA) {
		io_int(&MaxPoints,fp,f,"Max points for volterra");
	}
	io_exprs(f,fp);
	io_graph(f,fp);
	dump_transpose_info(fp,f);
	dump_h_stuff(fp,f);
	dump_aplot(fp,f);
	dump_torus(fp,f);
	dump_range(fp,f);
	dump_eqn(fp);
	fclose(fp);
}


void file_inf(void) {
	int ok;
	FILE *fp;
	char filename[XPP_MAX_NAME];
	sprintf(filename,"%s.pars",this_file);
	ping();
	if(!file_selector("Save info",filename,"*.pars*")) {
		return;
	}
	open_write_file(&fp,filename,&ok);
	if(!ok) {
		return;
	}
	redraw_params();
	do_info(fp);
	fclose(fp);
}


void io_ic_file(char *fn,int flag) {
	char fnx[XPP_MAX_NAME],c;
	int i,j=0;
	int chk=0;
	FILE *fp;
	char msg[XPP_MAX_NAME];
	for(i=0;i<strlen(fn);i++) {
		c=fn[i];
		if(c!=' ') {
			fnx[j]=c;
			j++;
		}
	}
	fnx[j]=0;
	if(flag==READEM) {
		fp=fopen(fnx,"r");
		if(fp==NULL) {
			err_msg("Cannot open file");
			return;
		}
		for(i=0;i<NODE;i++) {
			chk=fscanf(fp,"%lg",&last_ic[i]);
			if(chk!=1) {
				sprintf(msg,"Expected %d initial conditions but only found %d in %s.",NODE,i,fn);
				err_msg(msg);
				return;
			}
		}
		while (chk != EOF) {
			chk=fscanf(fp,"%lg",&last_ic[i]);
			if(chk!=EOF) {
				sprintf(msg,"Found more than %d initial conditions in %s.",NODE,fn);
				err_msg(msg);
				return;
			}
		}
		fclose(fp);
	}
}


void io_parameter_file(char *fn,int flag) {
	char fnx[XPP_MAX_NAME],c;
	int i,j=0;
	int np;
	FILE *fp;
	time_t ttt;
	for(i=6;i<strlen(fn);i++) {
		c=fn[i];
		if(c!=' ') {
			fnx[j]=c;
			j++;
		}
	}
	fnx[j]=0;
	if(flag==READEM) {
		fp=fopen(fnx,"r");
		if(fp==NULL) {
			err_msg("Cannot open file");
			return;
		}
		io_int(&np,fp,flag," ");
		if(np!=NUPAR) {
			printf("%d",np);
			printf("%d",NUPAR);
			err_msg("Incompatible parameters");
			fclose(fp);
			return;
		}
		io_parameters(flag,fp);
		fclose(fp);
		redo_stuff();
		return;
	}
	fp=fopen(fnx,"w");
	if(fp==NULL) {
		err_msg("Cannot open file");
		return;
	}
	io_int(&NUPAR,fp,flag,"Number params");
	io_parameters(flag,fp);
	ttt=time(0);
	fprintf(fp,"\n\nFile:%s\n%s",this_file, ctime(&ttt));
	fclose(fp);
}


int read_lunch(FILE *fp) {
	int f=READEM,ne,np,temp;
	char bob[DEFAULT_STRING_LENGTH];

	if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
		plintf("Couldnt read file %s", fp);
	}
	if(bob[0]=='#') {
		set_type=1;
		io_int(&ne,fp,f," ");
	} else {
		ne=atoi(bob);
		set_type=0;
	}
	/* io_int(&ne,fp,f); */
	io_int(&np,fp,f," ");
	if(ne!=NEQ || np!=NUPAR) {
		plintf("Set file has incompatible parameters\n");
		return 0;
	}
	io_numerics(f,fp);
	if(METHOD==METHOD_VOLTERRA) {
		io_int(&temp,fp,f," ");
		allocate_volterra(temp,1);
		MyStart=1;
	}
	check_delay();
	io_exprs(f,fp);
	io_graph(f,fp);
	if(set_type==1) {
		dump_transpose_info(fp,f);
		dump_h_stuff(fp,f);
		dump_aplot(fp,f);
		dump_torus(fp,f);
		dump_range(fp,f);
	}
	return 1;
}

/* --- Static functions --- */
static void do_info(FILE *fp) {
	int i;
	/* Must match enum Method. */
	static char *method[]={"Discrete","Euler","Mod. Euler",
						   "Runge-Kutta","Adams","Gear","Volterra","BackEul","QualRK",
						   "Stiff","CVode","DoPri5","DoPri8(3)","Rosenbrock","Symplectic"};
	int div,rem;
	int j;
	double z;
	char bob[200];
	char fstr[15];
	fprintf(fp,"File: %s \n\n Equations... \n",this_file);
	for(i=0;i<NEQ;i++) {
		if(i<NODE && METHOD != METHOD_DISCRETE) {
			strcpy(fstr,"d%s/dT=%s\n");
		}
		if(i<NODE && METHOD == METHOD_DISCRETE) {
			strcpy(fstr,"%s(n+1)=%s\n");
		}
		if(i>=NODE) {
			strcpy(fstr,"%s=%s\n");
		}
		fprintf(fp,fstr,uvar_names[i],ode_names[i]);
	}
	if(FIX_VAR>0) {
		fprintf(fp,"\nwhere ...\n");
		for(i=0;i<FIX_VAR;i++) {
			fprintf(fp,"%s = %s \n",fixinfo[i].name,fixinfo[i].value);
		}
	}
	if(NFUN>0) {
		fprintf(fp, "\nUser-defined functions:\n");
		user_fun_info(fp);
	}
	fprintf(fp,"\n\n Numerical parameters ...\n");
	fprintf(fp,"NJMP=%d  NMESH=%d METHOD=%s EVEC_ITER=%d \n",
			NJMP,NMESH,method[METHOD],EVEC_ITER);
	fprintf(fp,"BVP_EPS=%g,BVP_TOL=%g,BVP_MAXIT=%d \n",
			BVP_EPS,BVP_TOL,BVP_MAXIT);
	fprintf(fp,"DT=%g T0=%g TRANS=%g TEND=%g BOUND=%g DELAY=%g MaxPts=%d\n",
			DELTA_T,T0,TRANS,TEND,BOUND,DELAY,MaxPoints);
	fprintf(fp,"EVEC_ERR=%g, NEWT_ERR=%g HMIN=%g HMAX=%g TOLER=%g \n",
			EVEC_ERR,NEWT_ERR,HMIN,HMAX,TOLER);
	if(POIVAR==0) {
		strcpy(bob,"T");
	} else  {
		strcpy(bob,uvar_names[POIVAR-1]);
	}
	fprintf(fp,"POIMAP=%d POIVAR=%s POIPLN=%g POISGN=%d \n",
			POIMAP,bob,POIPLN,POISGN);

	fprintf(fp,"\n\n Delay strings ...\n");

	for(i=0;i<NODE;i++) {
		fprintf(fp,"%s\n",delay_string[i]);}
	fprintf(fp,"\n\n BCs ...\n");

	for(i=0;i<NODE;i++) {
		fprintf(fp,"0=%s\n",my_bc[i].string);
	}
	fprintf(fp,"\n\n ICs ...\n");

	for(i=0;i<NODE+NMarkov;i++) {
		fprintf(fp,"%s=%.16g\n",uvar_names[i],last_ic[i]);
	}
	fprintf(fp,"\n\n Parameters ...\n");
	div=NUPAR/4;
	rem=NUPAR%4;
	for(j=0;j<div;j++) {
		for(i=0;i<4;i++) {
			get_val(upar_names[i+4*j],&z);
			fprintf(fp,"%s=%.16g   ",upar_names[i+4*j],z);
		}
		fprintf(fp,"\n");
	}
	for(i=0;i<rem;i++) {
		get_val(upar_names[i+4*div],&z);
		fprintf(fp,"%s=%.16g   ",upar_names[i+4*div],z);
	}
	fprintf(fp,"\n");
}


static void dump_aplot(FILE *fp, int f) {
	char bob[DEFAULT_STRING_LENGTH];

	if(f==READEM) {
		if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
	} else {
		fprintf(fp,"# Array plot stuff\n");
		io_string(aplot.name,11,fp,f);
		io_int(&aplot.nacross ,fp,f,"NCols");
		io_int(&aplot.nstart ,fp,f,"Row 1");
		io_int(&aplot.ndown ,fp,f,"NRows");
		io_int(&aplot.nskip ,fp,f,"RowSkip");
		io_double(&aplot.zmin,fp,f,"Zmin");
		io_double(&aplot.zmax,fp,f,"Zmax");
	}
}


static void dump_eqn(FILE *fp) {
	int i;
	char fstr[15];
	fprintf(fp,"RHS etc ...\n");
	for(i=0;i<NEQ;i++) {
		if(i<NODE && METHOD != METHOD_DISCRETE) {
			strcpy(fstr,"d%s/dT=%s\n");
		}
		if(i<NODE && METHOD == METHOD_DISCRETE) {
			strcpy(fstr,"%s(n+1)=%s\n");
		}
		if(i>=NODE) {
			strcpy(fstr,"%s=%s\n");
		}
		fprintf(fp,fstr,uvar_names[i],ode_names[i]);
	}
	if(FIX_VAR>0) {
		fprintf(fp,"\nwhere ...\n");
		for(i=0;i<FIX_VAR;i++) {
			fprintf(fp,"%s = %s \n",fixinfo[i].name,fixinfo[i].value);
		}
	}
	if(NFUN>0) {
		fprintf(fp, "\nUser-defined functions:\n");
		user_fun_info(fp);
	}
}


/* f=1 to read and 0 to write */
static void dump_h_stuff(FILE *fp, int f) {
	char bob[DEFAULT_STRING_LENGTH];
	int i;

	if(f==READEM) {
		if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
	}else {
		fprintf(fp,"# Coupling stuff for H funs\n");
	}
	for(i=0;i<NODE ;i++) {
		io_string(coup_string[i],79,fp,f);
	}
}


static void dump_range(FILE *fp, int f) {
	char bob[DEFAULT_STRING_LENGTH];
	if(f==READEM) {
		if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
	} else {
		fprintf(fp,"# Range information\n");
	}
	io_string(eq_range.item,11,fp,f);
	io_int(&eq_range.col,fp,f,"eq-range stab col");
	io_int(&eq_range.shoot,fp,f,"shoot flag 1=on");
	io_int(&eq_range.steps,fp,f,"eq-range steps");
	io_double(&eq_range.plow,fp,f,"eq_range low");
	io_double(&eq_range.phigh,fp,f,"eq_range high");
	io_string(range.item,11,fp,f);
	io_string(range.item2,11,fp,f);
	io_int(&range.steps,fp,f,"Range steps");
	io_int(&range.cycle,fp,f,"Cycle color 1=on");
	io_int(&range.reset,fp,f,"Reset data 1=on");
	io_int(&range.oldic,fp,f,"Use old I.C.s 1=yes");
	io_double(&range.plow,fp,f,"Par1 low");
	io_double(&range.plow2,fp,f,"Par2 low");
	io_double(&range.phigh,fp,f,"Par1 high");
	io_double(&range.phigh2,fp,f,"Par2 high");
	dump_shoot_range(fp,f);
	if(f==READEM) {
		range.steps2=range.steps;
	}
}


static void dump_shoot_range(FILE *fp, int f) {
	io_string(shoot_range.item,11,fp,f);
	io_int(&shoot_range.side,fp,f,"BVP side");
	io_int(&shoot_range.cycle,fp,f,"color cycle flag 1=on");
	io_int(&shoot_range.steps,fp,f,"BVP range steps");
	io_double(&shoot_range.plow,fp,f,"BVP range low");
	io_double(&shoot_range.phigh,fp,f,"BVP range high");
}


static void dump_torus(FILE *fp, int f) {
	int i;
	char bob[DEFAULT_STRING_LENGTH];
	if(f==READEM) {
		if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
	} else {
		fprintf(fp,"# Torus information \n");
	}
	io_int(&TORUS,fp,f," Torus flag 1=ON");
	io_double(&TOR_PERIOD,fp,f,"Torus period");
	if(TORUS) {
		for(i=0;i<NEQ;i++) {
			io_int(&itor[i],fp,f,uvar_names[i]);
		}
	}
}


static void dump_transpose_info(FILE *fp, int f) {
	char bob[DEFAULT_STRING_LENGTH];
	if(f==READEM) {
		if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
	} else {
		fprintf(fp,"# Transpose variables etc\n");
	}

	io_string(my_trans.firstcol,11,fp,f);
	io_int(&my_trans.ncol,fp,f,"n columns");
	io_int(&my_trans.nrow,fp,f,"n rows");
	io_int(&my_trans.rowskip,fp,f,"row skip");
	io_int(&my_trans.colskip,fp,f,"col skip");
	io_int(&my_trans.row0,fp,f,"row 0");
}


static void io_double(double *z, FILE *fp, int f, char *ss) {
	char bob[DEFAULT_STRING_LENGTH];
	if(f==READEM) {
		if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
		*z=atof(bob);
	} else {
		fprintf(fp,"%.16g  %s\n",*z,ss);
	}
}


static void io_exprs(int f, FILE *fp) {
	int i;
	char temp[DEFAULT_STRING_LENGTH];
	double z;
	if(f==READEM && set_type==1) {
		if(fgets(temp,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		} /* skip a line */
	}
	if(f!=READEM) {
		fprintf(fp,"# Delays\n");
	}
	for(i=0;i<NODE;i++) {
		io_string(delay_string[i],100,fp,f);
	}
	if(f==READEM && set_type==1) {
		if(fgets(temp,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		} /* skip a line */
	}
	if(f!=READEM) {
		fprintf(fp,"# Bndry conds\n");
	}
	for(i=0;i<NODE;i++) {
		io_string(my_bc[i].string,100,fp,f);
	}
	if(f==READEM && set_type==1) {
		if(fgets(temp,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		} /* skip a line */
	}
	if(f!=READEM) {
		fprintf(fp,"# Old ICs\n");
	}
	for(i=0;i<NODE+NMarkov;i++) {
		io_double(&last_ic[i],fp,f,uvar_names[i]);
	}
	if(f==READEM && set_type==1) {
		if(fgets(temp,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		} /* skip a line */
	}
	if(f!=READEM) {
		fprintf(fp,"# Ending  ICs\n");
	}
	for(i=0;i<NODE+NMarkov;i++) {
		io_double(&MyData[i],fp,f,uvar_names[i]);
	}
	if(f==READEM && set_type==1) {
		if(fgets(temp,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		} /* skip a line */
	}
	if(f!=READEM) {
		fprintf(fp,"# Parameters\n");
	}
	for(i=0;i<NUPAR;i++) {
		if(f!=READEM) {
			get_val(upar_names[i],&z);
			io_double(&z,fp,f,upar_names[i]);
		} else {
			io_double(&z,fp,f," ");
			set_val(upar_names[i],z);
		}
	}
	if(f==READEM && Xup) {
		redraw_bcs();
		redraw_ics();
		redraw_delays();
		redraw_params();
	}
}


static void io_graph(int f, FILE *fp) {
	int j,k;
	char temp[DEFAULT_STRING_LENGTH];
	if(f==READEM && set_type==1) {
		if(fgets(temp,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		} /* skip a line */
	}
	if(f!=READEM) {
		fprintf(fp,"# Graphics\n");
	}
	for(j=0;j<3;j++) {
		for(k=0;k<3;k++) {
			io_double(&(MyGraph->rm[k][j]),fp,f,"rm");
		}
	}
	for(j=0;j<MAXPERPLOT;j++) {
		io_int(&(MyGraph->xv[j]),fp,f," ");
		io_int(&(MyGraph->yv[j]),fp,f," ");
		io_int(&(MyGraph->zv[j]),fp,f," ");
		io_int(&(MyGraph->line[j]),fp,f," ");
		io_int(&(MyGraph->color[j]),fp,f," ");
	}

	io_double(&(MyGraph->ZPlane),fp,f," ");
	io_double(&(MyGraph->ZView),fp,f," ");
	io_int(&(MyGraph->PerspFlag),fp,f," ");
	io_int(&(MyGraph->ThreeDFlag),fp,f,"3DFlag");
	io_int(&(MyGraph->TimeFlag),fp,f,"Timeflag");
	io_int(&(MyGraph->ColorFlag),fp,f,"Colorflag");
	io_int(&(MyGraph->grtype),fp,f,"Type");
	io_double(&(MyGraph->color_scale),fp,f,"color scale");
	io_double(&(MyGraph->min_scale),fp,f," minscale");

	io_double(&(MyGraph->xmax),fp,f," xmax");
	io_double(&(MyGraph->xmin),fp,f," xmin");
	io_double(&(MyGraph->ymax),fp,f," ymax");
	io_double(&(MyGraph->ymin),fp,f," ymin");
	io_double(&(MyGraph->zmax),fp,f," zmax");
	io_double(&(MyGraph->zmin),fp,f," zmin");
	io_double(&(MyGraph->xbar),fp,f, " ");
	io_double(&(MyGraph->dx  ),fp,f," ");
	io_double(&(MyGraph->ybar),fp,f," ");
	io_double(&(MyGraph->dy  ),fp,f," ");
	io_double(&(MyGraph->zbar),fp,f," ");
	io_double(&(MyGraph->dz  ),fp,f," ");

	io_double(&(MyGraph->Theta),fp,f," Theta");
	io_double(&(MyGraph->Phi),fp,f, " Phi");
	io_int(&(MyGraph->xshft),fp,f," xshft");
	io_int(&(MyGraph->yshft),fp,f," yshft");
	io_int(&(MyGraph->zshft),fp,f," zshft");
	io_double(&(MyGraph->xlo),fp,f," xlo");
	io_double(&(MyGraph->ylo),fp,f," ylo");
	io_double(&(MyGraph->oldxlo),fp,f," ");
	io_double(&(MyGraph->oldylo),fp,f," ");
	io_double(&(MyGraph->xhi),fp,f," xhi");
	io_double(&(MyGraph->yhi),fp,f," yhi");
	io_double(&(MyGraph->oldxhi),fp,f," ");
	io_double(&(MyGraph->oldyhi),fp,f," ");
	if(f==READEM && Xup) {
		redraw_the_graph();
	}
}


static void io_int(int *i, FILE *fp, int f, char *ss) {
	char bob[DEFAULT_STRING_LENGTH];
	if(f==READEM) {
		if(fgets(bob,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
		*i=atoi(bob);
	} else {
		fprintf(fp,"%d   %s\n",*i,ss);
	}
}


static void io_numerics(int f, FILE *fp) {
	int im = METHOD;
	char *method[]={"Discrete","Euler","Mod. Euler",
					"Runge-Kutta","Adams","Gear","Volterra","BackEul",
					"Qual RK","Stiff","CVode","DorPrin5","DorPri8(3)"};
	char *pmap[]={"Poincare None","Poincare Section","Poincare Max","Period"};
	char temp[DEFAULT_STRING_LENGTH];
	if(f==READEM && set_type==1) {
		if(fgets(temp,DEFAULT_STRING_LENGTH,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		} /* skip a line */
	}
	if(f!=READEM) {
		fprintf(fp,"# Numerical stuff\n");
	}
	io_int(&NJMP,fp,f," nout");
	io_int(&NMESH,fp,f," nullcline mesh");
	io_int(&im,fp,f,method[METHOD]);
	if(f==READEM) {
		do_meth();
		alloc_meth();
	}
	io_double(&TEND,fp,f,"total");
	io_double(&DELTA_T,fp,f,"DeltaT");
	io_double(&T0,fp,f,"T0");
	io_double(&TRANS,fp,f,"Transient");
	io_double(&BOUND,fp,f,"Bound");
	io_double(&HMIN,fp,f,"DtMin");
	io_double(&HMAX,fp,f,"DtMax");
	io_double(&TOLER,fp,f,"Tolerance");
	/* fix stuff concerning the tolerance */
	if(f==READEM) {
		if(set_type==1) {
			io_double(&ATOLER,fp,f,"Abs. Tolerance");
		} else {
			ATOLER=TOLER*10;
		}
	} else {
		io_double(&ATOLER,fp,f,"Abs. Tolerance");
	}
	io_double(&DELAY,fp,f,"Max Delay");
	io_int(&EVEC_ITER,fp,f,"Eigenvector iterates");
	io_double(&EVEC_ERR,fp,f,"Eigenvector tolerance");
	io_double(&NEWT_ERR,fp,f,"Newton tolerance");
	io_double(&POIPLN,fp,f,"Poincare plane");
	io_double(&BVP_TOL,fp,f,"Boundary value tolerance");
	io_double(&BVP_EPS,fp,f,"Boundary value epsilon");
	io_int(&BVP_MAXIT,fp,f,"Boundary value iterates");
	io_int(&POIMAP,fp,f,pmap[POIMAP]);

	io_int(&POIVAR,fp,f,"Poincare variable");
	io_int(&POISGN,fp,f,"Poincare sign");
	io_int(&SOS,fp,f,"Stop on Section");
	io_int(&DelayFlag,fp,f,"Delay flag");
	io_double(&MyTime,fp,f,"Current time");
	io_double(&LastTime,fp,f,"Last Time");
	io_int(&MyStart,fp,f,"MyStart");
	io_int(&INFLAG,fp,f,"INFLAG");
}


static void io_parameters(int f, FILE *fp) {
	int i,index;
	char junk[XPP_MAX_NAME];
	double z;
	for(i=0;i<NUPAR;i++) {
		if(f!=READEM) {
			get_val(upar_names[i],&z);
			io_double(&z,fp,f,upar_names[i]);
		} else {
			io_double(&z,fp,f," ");
			set_val(upar_names[i],z);
			if(!XPPBatch) {
				index=find_user_name(PARAMETER_BOX,upar_names[i]);
				if(index>=0) {
					sprintf(junk,"%.16g",z);
					set_edit_params(&ParamBox,index,junk);
					draw_one_box(ParamBox,index);
				}
			}
		}
	}
	if(!XPPBatch) {
		reset_sliders();
	}
}


static void io_string(char *s, int len, FILE *fp, int f) {
	int i;
	if(f==READEM) {
		if(fgets(s,len,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
		i=0;
		while(i<strlen(s)) {
			if(s[i]=='\n') {
				s[i]=0;
			}
			i++;
		}
	} else {
		fprintf(fp,"%s\n",s);
	}
}
