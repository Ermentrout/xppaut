#include <stdlib.h> 
#include <string.h>
/* command-line stuff for xpp */
#include <stdio.h>
#define NCMD 32 /* add new commands as needed  */

#define MAKEC 0
#define XORFX 1
#define SILENT 2 
#define CONVERT 3
#define NOICON 4
#define NEWSEED 5
#define ALLWIN 6
#define SETFILE 7
#define MSSTYLE 8
#define PWHITE 9
#define RUNNOW 10
#define BIGF 11
#define SMALLF 12
#define PARFILE 13
#define OUTFILE 14
#define ICFILE 15
#define FCOLOR 16
#define BCOLOR 17
#define BBITMAP 18
#define GRADS 19
#define MINWIDTH 20
#define MINHEIGHT 21
#define MWCOLOR 22
#define DWCOLOR 23
#define BELL 24
#define ITRNSETS 25
#define USET 26
#define RSET 27
#define INCLUDE 28
#define QSETS 29
#define QPARS 30
#define QICS 31

extern char big_font_name[100],small_font_name[100];

extern int RunImmediately;
extern int PaperWhite;
extern int MSStyle;
extern int got_file;
char setfilename[100];
char parfilename[100];
char icfilename[100];
char includefilename[100];

extern char UserBlack[8];
extern char UserWhite[8];
extern char UserMainWinColor[8];
extern char UserDrawWinColor[8];
extern char UserBGBitmap[100];
extern int UserGradients;
extern int UserMinWidth;
extern int UserMinHeight;
extern int UserMinHeight;
extern char UserOUTFILE[256];
extern int tfBell;
extern int use_intern_sets;
int select_intern_sets=0;


#define MAX_INTERN_SET 100
extern int Nintern_set;
int Nintern_2_use;

typedef struct {
  char *name;
  char *does;
  unsigned int use;
} INTERN_SET;


typedef struct {
   char *name;
   struct SET_NAME * next;
} SET_NAME;


SET_NAME *sets2use,*setsNOTuse;

extern INTERN_SET intern_set[MAX_INTERN_SET];


extern char batchout[256];
int loadsetfile=0;
int loadparfile=0;
int loadicfile=0;
int loadincludefile=0;
int querysets=0;
int querypars=0;
int queryics=0;
int dryrun=0;
extern char this_file[100];
extern int XPPBatch;
extern int xorfix;
extern int newseeed;
extern int silent;
extern int allwinvis;
extern int ConvertStyle;
int noicon=1;
int newseed=0;
typedef struct {
  char name[10];
  int len;

} VOCAB;

VOCAB my_cmd[NCMD]=
{
  "-m",3,
  "-xorfix",7,
  "-silent",7,
  "-convert",8,
  "-iconify",7,
  "-newseed",7,
  "-allwin",6,
  "-setfile",7,
  "-ee",3,
  "-white", 6,
  "-runnow",7,
  "-bigfont",8,
  "-smallfont",10,
  "-parfile",8,
  "-outfile",8,
  "-icfile",7,
  "-forecolor",10,
  "-backcolor",10,
  "-backimage",10,
  "-grads",6,
  "-width",6,
  "-height",7,
  "-mwcolor",8,
  "-dwcolor",8,
  "-bell",4,
  "-internset",10,
  "-uset",5,
  "-rset",5,
  "-include",8,
  "-qsets",6,
  "-qpars",6,
  "-qics",5
 };


is_set_name(set,nam)
SET_NAME *set;
char * nam;
{
	if (set==NULL){return(0);}
	SET_NAME *curr;
	
	curr=set;
	
	while(curr)
	{
		if (strcmp(curr->name,nam)==0)
		{
			return(1);
		}
		curr=curr->next;
	}
	
	return(0);
}

add_set(set,nam)
SET_NAME *set;
char * nam;
{
	if (!is_set_name(set,nam))
	{
		SET_NAME *curr;	
		curr = (SET_NAME *)malloc(sizeof(SET_NAME));
        	curr->name = nam;
		curr->next  = set;
		set=curr;
	}
	
	return(set);
}

rm_set(set,nam)
SET_NAME *set;
char *nam;
{
	SET_NAME *curr,*prev;	
	
	if (set==NULL){return;}
	
	curr=set;
	int i=1;
	while(curr)
	{
		if (strcmp(curr->name,nam)==0)
		{
			if (i==1)
			{
				set=curr->next;
			}
			else
			{
				prev->next=curr->next;
			}
			break;
		}
		prev = curr;
		i++;
	}
	
	
	return(set);
}


do_comline(argc,argv)
char **argv;
int argc;
{ 
 int i,k;

 silent = 0;
 got_file=0;
 xorfix=1;
 /*PaperWhite=0;
 */
 setfilename[0]=0;
 parfilename[0]=0;
 icfilename[0]=0;
 includefilename[0]=0;
 for(i=1;i<argc;i++){
   k=parse_it(argv[i]);
   if(k==1){
     strcpy(setfilename,argv[i+1]);
     i++;
     loadsetfile=1;
     
   }
   if(k==2){
     strcpy(small_font_name,argv[i+1]);
     i++;
   }
   if(k==3){
     strcpy(big_font_name,argv[i+1]);
     i++;
   } 
   if(k==4){
     strcat(parfilename,"!load ");
     strcat(parfilename,argv[i+1]);
     i++;
     loadparfile=1;
   }
   if(k==5){
     plintf(argv[i+1]);
     sprintf(batchout,argv[i+1]);
     sprintf(UserOUTFILE,argv[i+1]);
     i++;
   }
   if(k==6){
     strcat(icfilename,argv[i+1]);
     i++;
     loadicfile=1;
   }
   if(k==7){
     if (strlen(argv[i+1]) != 6)
     {
        plintf("Color must be given as hexadecimal string.\n");
	exit(-1);
     }
     set_option("FORECOLOR",argv[i+1]);
     i++;
     
   }
   if(k==8){
     if (strlen(argv[i+1]) != 6)
     {
        plintf("Color must be given as hexadecimal string.\n");
	exit(-1);
     }
     set_option("BACKCOLOR",argv[i+1]);
     i++;
     
   }
   if(k==9){
     /*strcpy(UserBGBitmap,argv[i+1]);
     */
     set_option("BACKIMAGE",argv[i+1]);
     i++;
   }
   if(k==10){
     set_option("GRADS",argv[i+1]);
     i++;
   }
   if(k==11){
     set_option("WIDTH",argv[i+1]);
     i++;
   }if(k==12){
     set_option("HEIGHT",argv[i+1]);
     i++;
   }if(k==13){
     if (strlen(argv[i+1]) != 6)
     {
        plintf("Color must be given as hexadecimal string.\n");
	exit(-1);
     }
     set_option("MWCOLOR",argv[i+1]);
     i++;
   }if(k==14){
     if (strlen(argv[i+1]) != 6)
     {
        plintf("Color must be given as hexadecimal string.\n");
	exit(-1);
     }
     set_option("DWCOLOR",argv[i+1]);
     i++;
   }
   if(k==15){
     set_option("BELL",argv[i+1]);
     i++;
   }
   if(k==16){
     use_intern_sets=atoi(argv[i+1]);
     select_intern_sets=1;
     i++;
   }  
   if(k==17){
     sets2use=add_set(sets2use,argv[i+1]);
     i++;
     select_intern_sets=1;
   }
   if(k==18){
     setsNOTuse=add_set(setsNOTuse,argv[i+1]);
     i++;
     select_intern_sets=1;
   } 
   if(k==19){
     strcpy(includefilename,argv[i+1]);
     i++;
     loadincludefile=1;
   }
 }
}


if_needed_select_sets()
{
	if(!select_intern_sets){return 1;}
	int j;
	for(j=0;j<Nintern_set;j++)
  	{
		intern_set[j].use=use_intern_sets;
		Nintern_2_use+=use_intern_sets;
		
		if (is_set_name(sets2use,intern_set[j].name))
		{
			printf("Internal set %s was included\n",intern_set[j].name);
			if (intern_set[j].use==0){Nintern_2_use++;}
			intern_set[j].use=1;
			
		}
		
		if (is_set_name(setsNOTuse,intern_set[j].name))
		{
			printf("Internal set %s was excluded\n",intern_set[j].name);
			if (intern_set[j].use==1){Nintern_2_use--;}
			intern_set[j].use=0;
		}
	}
	
	printf("A total of %d internal sets will be used\n",Nintern_2_use);
}


if_needed_load_set()
{
  FILE *fp;
  if(!loadsetfile)
    return 1;
  fp=fopen(setfilename,"r");
  if(fp==NULL){
    plintf("Couldn't load %s\n",setfilename);
    return 0;
  }
  read_lunch(fp);
  fclose(fp);
}



if_needed_load_par()
{
  FILE *fp;
  if(!loadparfile)
    return 1;
    plintf("Loading external parameter file: %s\n",parfilename);
    io_parameter_file(parfilename,1);
}


if_needed_load_ic()
{
  FILE *fp;
  if(!loadicfile)
    return 1;
    plintf("Loading external initial condition file: %s\n",icfilename);
    io_ic_file(icfilename,1);
}

parse_it(com)
     char *com;
{
  int j;
  for(j=0;j<NCMD;j++)
  {
  	if(strncmp(com,my_cmd[j].name,my_cmd[j].len)==0)
    	{
    		break;
  	}
  }
  if(j<NCMD){
    switch(j){
    case MAKEC:
      plintf(" C files are no longer part of this version. \n Sorry \n");
      break;
    case SILENT:
      XPPBatch=1;
      break;
    case XORFX:
      xorfix=0;
      break;
    case CONVERT:
      ConvertStyle=1;
      break;
    case NOICON:
      noicon=0;
      break;
    case NEWSEED:
      plintf("Random number seed changed\n");
      newseed=1;
      break;  
    case ALLWIN:
      allwinvis=1;
      break;
    case MSSTYLE:
      MSStyle=1;
      break;
    case PWHITE:
      PaperWhite=1;
      break;
    case RUNNOW:
      RunImmediately=1;
      break;
    case SETFILE:
      return 1;
    case SMALLF:
      return 2;
    case BIGF:
      return 3;
    case PARFILE:
      return 4;
    case OUTFILE:
      return 5; 
    case ICFILE:
      return 6;
    case FCOLOR:
      return 7;
    case BCOLOR:
      return 8;
    case BBITMAP:
      return 9;
    case GRADS:
      return 10;
    case MINWIDTH:
      return 11;
    case MINHEIGHT:
      return 12;
    case MWCOLOR:
      return 13;
    case DWCOLOR:
      return 14;
    case BELL:
      return 15;
    case ITRNSETS:
      return 16;
    case USET:
      return 17;
    case RSET:
      return 18;
    case INCLUDE:
      return 19;
    case QSETS:
      XPPBatch=1;
      querysets=1;
      dryrun=1;
      break;
    case QPARS: 
      XPPBatch=1;
      querypars=1;
      dryrun=1;
      break;
    case QICS:
      XPPBatch=1;
      queryics=1;
      dryrun=1;
      break;
    }
  }
  else {
    if(com[0]=='-'||got_file==1){ 
      plintf("Problem reading option %s\n",com);
      plintf("\nUsage: xppaut filename [options ...]\n\n");
      plintf("Options:\n");
      plintf("  -silent                Batch run without the interface and dump solutions to a file\n");
      plintf("  -xorfix                Work-around for exclusive Or with X on some monitors/graphics setups\n");
      plintf("  -convert               Convert old style ODE files (e.g. phaseplane) to new ODE style\n");
      plintf("  -newseed               Randomizes the random number generator which will often use the same seed\n");
      plintf("  -ee                    Emulates shortcuts of Evil Empire style (MS)\n");
      plintf("  -allwin                Brings XPP up with all the windows visible\n");
      plintf("  -white                 Uses white screen instead of black\n");
      plintf("  -setfile <filename>    Loads the set file before starting up\n");
      plintf("  -runnow                Runs ode file immediately upon startup (implied by -silent)\n");
      plintf("  -bigfont <font>        Use the big font whose filename is given\n");
      plintf("  -smallfont <font>      Use the small font whose filename is given\n");
      plintf("  -parfile <filename>    Load parameters from the named file\n");
      plintf("  -outfile <filename>    Send output to this file (default is output.dat)\n");
      plintf("  -icfile <filename>     Load initial conditions from the named file\n");
      plintf("  -forecolor <######>    Hexadecimal color (e.g. 000000) for foreground\n");
      plintf("  -backcolor <######>    Hexadecimal color (e.g. EDE9E3) for background\n");
      plintf("  -backimage <filename>  Name of bitmap file (.xbm) to load in background\n");
      plintf("  -mwcolor <######>      Hexadecimal color (e.g. 808080) for main window\n");
      plintf("  -dwcolor <######>      Hexadecimal color (e.g. FFFFFF) for drawing window\n");
      plintf("  -grads < 1 | 0 >       Color gradients will | won't be used\n"); 
      plintf("  -width N               Minimum width in pixels of main window\n");
      plintf("  -height N              Minimum height in pixels of main window\n");
      plintf("  -bell < 1 | 0 >        Events will | won't trigger system bell\n");
      plintf("  -internset < 1 | 0 >   Internal sets will | won't be run during batch run\n");
      plintf("  -uset <setname>        Named internal set will be run during batch run\n");
      plintf("  -rset <setname>        Named internal set will not be run during batch run\n");
      plintf("  -include <filename>    Named file will be included (see #include directive)\n");
      plintf("  -qsets                 Query internal sets (output saved to OUTFILE)\n");
      plintf("  -qpars                 Query parameters (output saved to OUTFILE)\n");
      plintf("  -qics                  Query initial conditions (output saved to OUTFILE)\n");
     
      plintf("\n");
      plintf("Environment variables:\n");
      plintf("  XPPHELP                Path to XPPAUT documentation file <xpphelp.html>\n");
      plintf("  XPPBROWSER             Web browser (e.g. /usr/bin/firefox)\n");
      plintf("  XPPSTART               Path to start looking for ODE files\n");
      plintf("\n");
     exit(0);
    }
    else {
      strcpy(this_file,com);
      got_file=1;
    }
  }
  return 0;
}







