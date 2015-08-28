#include <stdlib.h> 
 /*       This saves and loads parameters etc when you are out
         to lunch ....

*/
#include <stdio.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "xpplim.h"
#include "struct.h"
#include "shoot.h"


#define READEM 1
#define VOLTERRA 6

double atof();

 extern Window main_win;
 extern GRAPH *MyGraph;
 extern BC_STRUCT my_bc[MAXODE];

typedef struct {
  char *name,*value;
} FIXINFO;

extern FIXINFO fixinfo[MAXODE];
extern int FIX_VAR,NFUN;
 
 extern int NJMP,NMESH,METHOD,NODE,POIMAP,POIVAR,POISGN,SOS,INFLAG,NMarkov;
 extern int NUPAR,NEQ,BVP_MAXIT,EVEC_ITER,DelayFlag,MyStart;
 extern double last_ic[MAXODE],MyData[MAXODE],MyTime,LastTime;
 extern double TEND,DELTA_T,T0,TRANS,BOUND,HMIN,HMAX,TOLER,ATOLER,DELAY;
 extern double POIPLN,EVEC_ERR,NEWT_ERR;
extern double BVP_TOL,BVP_EPS;
extern int MaxPoints;

 extern char upar_names[MAXPAR][11],this_file[100],delay_string[MAXODE][80];
 extern char uvar_names[MAXODE][12]; 
 extern char *ode_names[MAXODE],*fix_names[MAXODE];


file_inf()
{
  int ok;
  FILE *fp;
 /*char filename[256];*/
 char filename[XPP_MAX_NAME];
 sprintf(filename,"%s.pars",this_file);
 ping();
 if(!file_selector("Save info",filename,"*.pars*"))return;
 /* if(new_string("Filename: ",filename)==0)return; */
  open_write_file(&fp,filename,&ok); 
   if(!ok)return;
    redraw_params();
    do_info(fp);
    fclose(fp);
}

 
do_info(fp)
FILE *fp;
{
 int i;
 static char *method[]={"Discrete","Euler","Mod. Euler",
	"Runge-Kutta","Adams","Gear","Volterra","BackEul","QualRK",
         "Stiff","CVode"};
 int div,rem;
 int j;
 double z;
 char bob[200];
 char fstr[15];
 fprintf(fp,"File: %s \n\n Equations... \n",this_file);
 for(i=0;i<NEQ;i++){
   if(i<NODE&&METHOD>0)strcpy(fstr,"d%s/dT=%s\n");
    if(i<NODE&&METHOD==0)strcpy(fstr,"%s(n+1)=%s\n");
     if(i>=NODE)strcpy(fstr,"%s=%s\n");
     fprintf(fp,fstr,uvar_names[i],ode_names[i]);
 }

 if(FIX_VAR>0){
   fprintf(fp,"\nwhere ...\n");
   for(i=0;i<FIX_VAR;i++)
     fprintf(fp,"%s = %s \n",fixinfo[i].name,fixinfo[i].value);
 }
  if(NFUN>0){
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
         if(POIVAR==0)strcpy(bob,"T");
	   else strcpy(bob,uvar_names[POIVAR-1]);
 fprintf(fp,"POIMAP=%d POIVAR=%s POIPLN=%g POISGN=%d \n",
        POIMAP,bob,POIPLN,POISGN);
 
 fprintf(fp,"\n\n Delay strings ...\n");
 
 for(i=0;i<NODE;i++)fprintf(fp,"%s\n",delay_string[i]);
  fprintf(fp,"\n\n BCs ...\n");
 
 for(i=0;i<NODE;i++)fprintf(fp,"0=%s\n",my_bc[i].string);
 fprintf(fp,"\n\n ICs ...\n");
 
 for(i=0;i<NODE+NMarkov;i++)fprintf(fp,"%s=%.16g\n",uvar_names[i],last_ic[i]);
 fprintf(fp,"\n\n Parameters ...\n");
 div=NUPAR/4;
 rem=NUPAR%4;
 for(j=0;j<div;j++){
   for(i=0;i<4;i++)
     {
       get_val(upar_names[i+4*j],&z);
       fprintf(fp,"%s=%.16g   ",upar_names[i+4*j],z);
     }
   fprintf(fp,"\n");
 }
    for(i=0;i<rem;i++){
      get_val(upar_names[i+4*div],&z);
      fprintf(fp,"%s=%.16g   ",upar_names[i+4*div],z);
    }
   
 fprintf(fp,"\n");
}
  


 

do_lunch(f)
int f;
{
 int ne,np,ok,temp;
 char bob[80];
 FILE *fp;
 /*char filename[256];*/
 char filename[XPP_MAX_NAME];
 sprintf(filename,"%s.set",this_file);
 ping();
 if(!file_selector("Read/write set file",filename,"*.set*"))return;
 /* if(new_string("Filename: ",filename)==0)return; */
 if(f==READEM){
   fp=fopen(filename,"r");
   if(fp==NULL){
     err_msg("Cannot open file");
     return;
   }
   io_int(&ne,fp,f);
   io_int(&np,fp,f);
   if(ne!=NEQ||np!=NUPAR){
     err_msg("Incompatible parameters");
     fclose(fp);
     return;
   }
   io_numerics(f,fp);
   if(METHOD==VOLTERRA){
     io_int(&temp,fp,f);
     allocate_volterra(temp,1);
     MyStart=1;
   }
   chk_delay();
   io_exprs(f,fp);
   io_graph(f,fp);
   fclose(fp);
   return;
 }
 open_write_file(&fp,filename,&ok); 
   if(!ok)return;
 redraw_params();
 io_int(&NEQ,fp,f);
 io_int(&NUPAR,fp,f);
 io_numerics(f,fp);
 if(METHOD==VOLTERRA){
     io_int(&MaxPoints,fp,f);
     }
   io_exprs(f,fp);
   io_graph(f,fp);
   fclose(fp);
}
 
 



io_numerics(f,fp)
int f;
FILE *fp;
{
io_int(&NJMP,fp,f);
io_int(&NMESH,fp,f);
io_int(&METHOD,fp,f);
if(f==READEM)do_meth();
io_double(&TEND,fp,f);
io_double(&DELTA_T,fp,f);
io_double(&T0,fp,f);
io_double(&TRANS,fp,f);
io_double(&BOUND,fp,f);
io_double(&HMIN,fp,f);
io_double(&HMAX,fp,f);
io_double(&TOLER,fp,f);
io_double(&DELAY,fp,f);
io_int(&EVEC_ITER,fp,f);
io_double(&EVEC_ERR,fp,f);
io_double(&NEWT_ERR,fp,f);
io_double(&POIPLN,fp,f);
io_double(&BVP_TOL,fp,f);
io_double(&BVP_EPS,fp,f);
io_int(&BVP_MAXIT,fp,f);
io_int(&POIMAP,fp,f);
io_int(&POIVAR,fp,f);
io_int(&POISGN,fp,f);
io_int(&SOS,fp,f);
io_int(&DelayFlag,fp,f);
io_double(&MyTime,fp,f);
io_double(&LastTime,fp,f);
io_int(&MyStart,fp,f);
io_int(&INFLAG,fp,f);
if(f==READEM)
  alloc_meth();
}


io_exprs(f,fp)
int f;
FILE *fp;
{
 int i;
 double z;
 for(i=0;i<NODE;i++)io_string(delay_string[i],100,fp,f);
 for(i=0;i<NODE;i++)io_string(my_bc[i].string,100,fp,f);
 for(i=0;i<NODE+NMarkov;i++)io_double(&last_ic[i],fp,f);
 for(i=0;i<NODE+NMarkov;i++)io_double(&MyData[i],fp,f);
 
 for(i=0;i<NUPAR;i++){
  if(f!=READEM){
    get_val(upar_names[i],&z);
    io_double(&z,fp,f);
  }
  else {
    io_double(&z,fp,f);
    set_val(upar_names[i],z);
  }
}
  
   
 if(f==READEM){
   redraw_bcs();
   redraw_ics();
   redraw_delays();
   redraw_params();
 }
}
   
 



io_graph(f,fp)
int f;
FILE *fp;
{
 int n,j,k;
 for(j=0;j<3;j++)
   for(k=0;k<3;k++)
     io_double(&(MyGraph->rm[k][j]),fp,f);
 for(j=0;j<MAXPERPLOT;j++){
        io_int(&(MyGraph->xv[j]),fp,f);
	io_int(&(MyGraph->yv[j]),fp,f);
	io_int(&(MyGraph->zv[j]),fp,f);
        io_int(&(MyGraph->line[j]),fp,f);
	io_int(&(MyGraph->color[j]),fp,f);
        }

    io_double(&(MyGraph->ZPlane),fp,f);
    io_double(&(MyGraph->ZView),fp,f);
    io_int(&(MyGraph->PerspFlag),fp,f);
    io_int(&(MyGraph->ThreeDFlag),fp,f);
    io_int(&(MyGraph->TimeFlag),fp,f);
    io_int(&(MyGraph->ColorFlag),fp,f);
    io_int(&(MyGraph->grtype),fp,f);
    io_double(&(MyGraph->color_scale),fp,f);
    io_double(&(MyGraph->min_scale),fp,f);

    io_double(&(MyGraph->xmax),fp,f);
    io_double(&(MyGraph->xmin),fp,f);
    io_double(&(MyGraph->ymax),fp,f);
    io_double(&(MyGraph->ymin),fp,f);
    io_double(&(MyGraph->zmax),fp,f);
    io_double(&(MyGraph->zmin),fp,f);
    io_double(&(MyGraph->xbar),fp,f);
    io_double(&(MyGraph->dx  ),fp,f);
    io_double(&(MyGraph->ybar),fp,f);
    io_double(&(MyGraph->dy  ),fp,f);
    io_double(&(MyGraph->zbar),fp,f);
    io_double(&(MyGraph->dz  ),fp,f);

    io_double(&(MyGraph->Theta),fp,f);
    io_double(&(MyGraph->Phi),fp,f);
    io_int(&(MyGraph->xshft),fp,f);
    io_int(&(MyGraph->yshft),fp,f);
    io_int(&(MyGraph->zshft),fp,f);
    io_double(&(MyGraph->xlo),fp,f);
    io_double(&(MyGraph->ylo),fp,f);
    io_double(&(MyGraph->oldxlo),fp,f);
    io_double(&(MyGraph->oldylo),fp,f);
    io_double(&(MyGraph->xhi),fp,f);
    io_double(&(MyGraph->yhi),fp,f);
    io_double(&(MyGraph->oldxhi),fp,f);
    io_double(&(MyGraph->oldyhi),fp,f);
    if(f==READEM)redraw_the_graph();
}

 
io_int(i,fp,f)
int *i,f;
FILE *fp;
{
 char bob[81];
 if(f==READEM){
   fgets(bob,80,fp);
   *i=atoi(bob);
 }
 else
 fprintf(fp,"%d\n",*i);
}

io_double(z,fp,f)
int f;
FILE *fp;
double *z;
{
char bob[81];
 if(f==READEM){
   fgets(bob,80,fp);
   *z=atof(bob);
 }
 else
 fprintf(fp,"%.16g\n",*z); 
}

io_float(z,fp,f)
int f;
FILE *fp;
float *z;
{
 char bob[81];
if(f==READEM){
   fgets(bob,80,fp);
   *z=(float)atof(bob);
 }
 else
 fprintf(fp,"%.16g\n",*z); 
}

io_int_array(k,n,fp,f)
int n,f,*k;
FILE *fp;
{
int i;
 for(i=0;i<n;i++)io_int(fp,&k[i],f);
}

io_double_array(z,n,fp,f)
double *z;
int n,f;
FILE *fp;
{
 int i;
 for(i=0;i<n;i++)io_double(fp,&z[i],f);

}

io_string(s,len,fp,f)
FILE *fp;
char *s;
int f,len;
{
 int i;
 if(f==READEM){
   fgets(s,len,fp);
   i=0;
   while(i<strlen(s)){
     if(s[i]=='\n')s[i]=0;
     i++;
   }
 }
 else 
   fprintf(fp,"%s\n",s);
}
















