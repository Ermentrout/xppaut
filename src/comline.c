/* command-line stuff for xpp */
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "comline.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "form_ode.h"
#include "ggets.h"
#include "integrate.h"
#include "load_eqn.h"
#include "lunch-new.h"
#include "main.h"
#include "xpplim.h"
/* --- Types --- */
typedef struct {
	char *name;
	struct SET_NAME * next;
} SET_NAME;


/* --- Forward declarations --- */
static SET_NAME *add_set(SET_NAME *set, char *nam);
static int is_set_name(SET_NAME *set, char *nam);
static int parse_one(int argc, char **argv);
static void show_help(void);

/* --- Data --- */
static char icfilename[XPP_MAX_NAME];
static char parfilename[XPP_MAX_NAME];
static char readsetfile[XPP_MAX_NAME];
static char setfilename[XPP_MAX_NAME];
static char externaloptionsstring[MAXEXPLEN];

static int externaloptionsflag=0;
static int use_intern_sets = -1;
static SET_NAME *sets2use,*setsNOTuse;

char includefilename[MaxIncludeFiles][XPP_MAX_NAME];
int NincludedFiles=0;
int Nintern_2_use=0;
int querysets=0;
int querypars=0;
int queryics=0;
int dryrun=0;
int noicon=1;
int newseed=0;


/* --- Functions --- */
/**
 * Handles the comando line options in argv[0].
 *
 * @param argc The number of arguments in argv. Must be >= 1.
 * @param argv Argument strings.
 */
void do_comline(int argc, char **argv) {
	silent = 0;
	xorfix=1;

	this_file[0]   = 0;
	setfilename[0] = 0;
	parfilename[0] = 0;
	icfilename[0]  = 0;
	/* Parse the individual inputs.
	 * parse_one returns the number of read inputs, so we have to increment
	 * by its return value
	 */
	for (int i = 1; i < argc; i++) {
		i += parse_one(argc - i, &argv[i]);
	}
}


/**
 * Reads in an external options file.
 *
 * @return 1 if there was either no file or it was correctly read
 *		   0 otherwise
 */
int if_needed_load_ext_options(void) {
	FILE *fp;
	char myopts[MAXEXPLEN];
	char myoptsx[MAXEXPLEN+2];
	if(externaloptionsflag==0) {
		return 1;
	}
	if(externaloptionsflag==1) {
		fp=fopen(readsetfile,"r");
		if(fp==NULL) {
			plintf("%s external set not found\n",readsetfile);
			return 0;
		}
		if(fgets(myopts,MAXEXPLEN,fp)==NULL) {
			plintf("Couldnt read file %s", fp);
		}
		sprintf(myoptsx,"$ %s",myopts);
		plintf("Got this string: {%s}\n",myopts);
		extract_action(myoptsx);
		fclose(fp);
		return 1;
	}
	if(externaloptionsflag==2) {
		sprintf(myoptsx,"$ %s",externaloptionsstring);
		extract_action(myoptsx);
		return 1;
	}
	return 0;
}


/**
 * Reads in an initial conditions file.
 *
 * @return 1 if there was either no file or it was correctly read
 *		   0 otherwise
 */
int if_needed_load_ic(void) {
	if(!icfilename[0]) {
		return 1;
	}
	plintf("Loading external initial condition file: %s\n",icfilename);
	io_ic_file(icfilename, READEM);
	return(1);
}


/**
 * Reads in an external parameter file.
 *
 * @return 1 if there was either no file or it was correctly read
 *		   0 otherwise
 */
int if_needed_load_par(void) {
	if(!parfilename[0]) {
		return 1;
	}
	plintf("Loading external parameter file: %s\n", parfilename);
	io_parameter_file(parfilename, READEM);
	return 1;
}


/**
 * Reads in a settings file.
 *
 * @return 1 if there was either no file or it was correctly read
 *		   0 otherwise
 */
int if_needed_load_set(void) {
	FILE *fp;
	if(!setfilename[0]) {
		return 1;
	}
	fp=fopen(setfilename,"r");
	if(fp == NULL) {
		plintf("Couldn't load %s\n",setfilename);
		return 0;
	}
	read_lunch(fp);
	fclose(fp);
	return 1;
}


/**
 * Reads in selected SET_NAMEs.
 *
 * @return 1 if there was either no file or it was correctly read
 *		   0 otherwise
 */
int if_needed_select_sets(void) {
	if(use_intern_sets == -1 &&
	   sets2use == NULL &&
	   setsNOTuse == NULL) {
		return 1;
	}
	for(int j=0; j < Nintern_set; j++) {
		intern_set[j].use = use_intern_sets;
		Nintern_2_use += use_intern_sets;

		if(is_set_name(sets2use, intern_set[j].name)) {
			plintf("Internal set %s was included\n",intern_set[j].name);
			if(intern_set[j].use==0) {
				Nintern_2_use++;
			}
			intern_set[j].use = 1;
		}

		if(is_set_name(setsNOTuse,intern_set[j].name))	{
			plintf("Internal set %s was excluded\n",intern_set[j].name);
			if(intern_set[j].use==1) {
				Nintern_2_use--;
			}
			intern_set[j].use=0;
		}
	}
	plintf("A total of %d internal sets will be used\n",Nintern_2_use);
	return 1;
}


/* --- Static functions --- */
/**
 * Checks if a given name is part of a SET_NAME.
 *
 * @param set  Pointer to the SET_NAME set
 * @param name Option name that is searched
 * @return 1 if the name is found 0 otherwise.
 */
static int is_set_name(SET_NAME *set, char *name) {
	if(set==NULL) {
		return(0);
	}
	SET_NAME *curr;
	curr=set;
	while(curr) {
		if(strcmp(curr->name,name)==0) {
			return(1);
		}
		curr=(SET_NAME*)curr->next;
	}
	return(0);
}


/**
 * Handles the options in argv[0].
 *
 * @param argc The number of arguments in argv. Must be >= 1.
 * @param argv Argument strings.
 * @return the number of arguments used.
 */
static int parse_one(int argc, char **argv) {
	if (argv[0][0] != '-') {
		strncpy(this_file, argv[0], sizeof(this_file));
		this_file[sizeof(this_file)-1] = '\0';
	}
	/* Options taking no parameter. */
	else if (!strcmp(argv[0],"-mkplot")) {
		MakePlotFlag = 1;
	}
	else if (!strcmp(argv[0],"-silent")) {
		XPPBatch = 1;
	}
	else if (!strcmp(argv[0],"-xorfix")) {
		xorfix = 0;
	}
	else if (!strcmp(argv[0],"-convert")) {
		ConvertStyle = 1;
	}
	else if (!strcmp(argv[0],"-iconify")) {
		noicon = 0;
	}
	else if (!strcmp(argv[0],"-newseed")) {
		newseed = 1;
	}
	else if (!strcmp(argv[0],"-allwin")) {
		allwinvis = 1;
	}
	else if (!strcmp(argv[0],"-ee")) {
		MSStyle = 1;
	}
	else if (!strcmp(argv[0],"-runnow")) {
		RunImmediately = 1;
	}
	else if (!strcmp(argv[0],"-version")) {
		printf("XPPAUT Version %s\n", PACKAGE_VERSION);
		exit(0);
	}
	else if (!strcmp(argv[0],"-noout")) {
		SuppressOut = 1;
	}
	else if (!strcmp(argv[0],"-qsets")) {
		XPPBatch = 1;
		querysets = 1;
		dryrun = 1;
	}
	else if (!strcmp(argv[0],"-qpars")) {
		XPPBatch = 1;
		querypars = 1;
		dryrun = 1;
	}
	else if (!strcmp(argv[0],"-qics")) {
		XPPBatch = 1;
		queryics = 1;
		dryrun = 1;
	} else {
		if (argc < 2) {
			plintf("Problem reading option %s\n", argv[0]);
			show_help();
			exit(1);
		}
		/* Options that take one parameter. */
		else if (!strcmp(argv[0],"-setfile")) {
			strncpy(setfilename, argv[1], sizeof(setfilename));
			setfilename[sizeof(setfilename)-1] = '\0';
		}
		else if (!strcmp(argv[0],"-smallfont")) {
			set_option("SMALLFONT", argv[1]);
		}
		else if (!strcmp(argv[0],"-bigfont")) {
			set_option("BIGFONT", argv[1]);
		}
		else if (!strcmp(argv[0],"-parfile")) {
			snprintf(parfilename, sizeof(parfilename), "!load %s", argv[1]);
		}
		else if (!strcmp(argv[0],"-outfile")) {
			strncpy(batchout, argv[1], sizeof(batchout));
			batchout[sizeof(batchout)-1] = '\0';
			strncpy(UserOUTFILE, argv[1], sizeof(UserOUTFILE));
			UserOUTFILE[sizeof(UserOUTFILE)-1] = '\0';
		}
		else if (!strcmp(argv[0],"-icfile")) {
			strncpy(icfilename, argv[1], sizeof(icfilename));
			icfilename[sizeof(icfilename)-1] = '\0';
		}
		else if (!strcmp(argv[0],"-forecolor")) {
			set_option("FORECOLOR", argv[1]);
		}
		else if (!strcmp(argv[0],"-backcolor")) {
			set_option("BACKCOLOR", argv[1]);
		}
		else if (!strcmp(argv[0],"-backimage")) {
			set_option("BACKIMAGE", argv[1]);
		}
		else if (!strcmp(argv[0],"-grads")) {
			set_option("GRADS", argv[1]);
		}
		else if (!strcmp(argv[0],"-width")) {
			set_option("WIDTH", argv[1]);
		}
		else if (!strcmp(argv[0],"-height")) {
			set_option("HEIGHT", argv[1]);
		}
		else if (!strcmp(argv[0],"-mwcolor")) {
			set_option("MWCOLOR", argv[1]);
		}
		else if (!strcmp(argv[0],"-dwcolor")) {
			set_option("DWCOLOR", argv[1]);
		}
		else if (!strcmp(argv[0],"-bell")) {
			set_option("BELL", argv[1]);
		}
		else if (!strcmp(argv[0],"-internset")) {
			use_intern_sets = atoi(argv[1]);
		}
		else if (!strcmp(argv[0],"-uset")) {
			sets2use = add_set(sets2use, argv[1]);
		}
		else if (!strcmp(argv[0],"-rset")) {
			setsNOTuse = add_set(setsNOTuse, argv[1]);
		}
		else if (!strcmp(argv[0],"-include")) {
			strncpy(includefilename[NincludedFiles], argv[1], sizeof(includefilename[NincludedFiles]));
			includefilename[NincludedFiles][sizeof(includefilename[NincludedFiles])-1] = '\0';
			NincludedFiles++;
		}
		else if (!strcmp(argv[0],"-quiet")) {
			set_option("QUIET", argv[1]);
		}
		else if (!strcmp(argv[0],"-logfile")) {
			set_option("LOGFILE", argv[1]);
		}
		else if (!strcmp(argv[0],"-anifile")) {
			strcpy(anifile, argv[1]);
		}
		else if (!strcmp(argv[0],"-plotfmt")) {
			set_option("PLOTFMT", argv[1]);
		}
		else if (!strcmp(argv[0],"-dfdraw")) {
			set_option("DFDRAW", argv[1]);
		}
		else if (!strcmp(argv[0],"-ncdraw")) {
			set_option("NCDRAW", argv[1]);
		} else {
			plintf("Problem reading option %s\n", argv[0]);
			show_help();
			exit(1);
		}
		return 2;
	}
	return 1;
}


/**
 * Shows xppaut help.
 */
static void show_help(void) {
	plintf("\nUsage: xppaut [filename] [options ...]\n\n");
	plintf("Options:\n");
	plintf("  -silent                Batch run without the interface and "
		   "dump solutions to a file\n");
	plintf("  -xorfix                Work-around for exclusive Or with X on "
		   "some monitors/graphics setups\n");
	plintf("  -convert               Convert old style ODE files (e.g. "
		   "phaseplane) to new ODE style\n");
	plintf("  -newseed               Randomizes the random number generator "
		   "which will often use the same seed\n");
	plintf("  -ee                    Emulates shortcuts of Evil Empire style "
		   "(MS)\n");
	plintf("  -allwin                Brings XPP up with all the windows "
		   "visible\n");
	plintf("  -white                 Uses white screen instead of black\n");
	plintf("  -setfile <filename>    Loads the set file before starting up\n");
	plintf("  -runnow                Runs ode file immediately upon startup "
		   "(implied by -silent)\n");
	plintf("  -bigfont <font>        Use the big font whose filename is "
		   "given\n");
	plintf("  -smallfont <font>      Use the small font whose filename is "
		   "given\n");
	plintf("  -parfile <filename>    Load parameters from the named file\n");
	plintf("  -outfile <filename>    Send output to this file (default is "
		   "output.dat)\n");
	plintf("  -icfile <filename>     Load initial conditions from the named "
		   "file\n");
	plintf("  -forecolor <######>    Hexadecimal color (e.g. 000000) for "
		   "foreground\n");
	plintf("  -backcolor <######>    Hexadecimal color (e.g. EDE9E3) for "
		   "background\n");
	plintf("  -backimage <filename>  Name of bitmap file (.xbm) to load in "
		   "background\n");
	plintf("  -mwcolor <######>      Hexadecimal color (e.g. 808080) for "
		   "main window\n");
	plintf("  -dwcolor <######>      Hexadecimal color (e.g. FFFFFF) for "
		   "drawing window\n");
	plintf("  -grads < 1 | 0 >       Color gradients will | won't be used\n");
	plintf("  -width N               Minimum width in pixels of main window\n");
	plintf("  -height N              Minimum height in pixels of main window\n");
	plintf("  -bell < 1 | 0 >        Events will | won't trigger system bell\n");
	plintf("  -internset < 1 | 0 >   Internal sets will | won't be run "
		   "during batch run\n");
	plintf("  -uset <setname>        Named internal set will be run during "
		   "batch run\n");
	plintf("  -rset <setname>        Named internal set will not be run "
		   "during batch run\n");
	plintf("  -include <filename>    Named file will be included (see "
		   "#include directive)\n");
	plintf("  -qsets                 Query internal sets (output saved to "
		   "OUTFILE)\n");
	plintf("  -qpars                 Query parameters (output saved to "
		   "OUTFILE)\n");
	plintf("  -qics                  Query initial conditions (output saved "
		   "to OUTFILE)\n");
	plintf("  -quiet <1 |0>          Do not print *anything* out to console\n");
	plintf("  -logfile <filename>    Print console output to specified "
		   "logfile \n");
	plintf("  -anifile <filename>    Load an animation code file (.ani) \n");
	plintf("  -plotfmt <svg|ps>      Set Batch plot format\n");
	plintf("  -mkplot                Do a plot in batch mode \n");
	plintf("  -ncdraw 1              Draw nullclines in batch \n");
	plintf("  -dfdraw <1|2>          Draw dfields in batch  \n");
	plintf("  -version               Print XPPAUT version and exit \n");

	plintf("\n");

	plintf("Environment variables:\n");
	plintf("  XPPHELP                Path to XPPAUT documentation file "
		   "<xpphelp.html>\n");
	plintf("  XPPBROWSER             Web browser (e.g. /usr/bin/firefox)\n");
	plintf("  XPPSTART               Path to start looking for ODE files\n");
	plintf("\n");
}

/**
 * Adds a option name to SET_NAME set.
 *
 * @param set the SET_NAME we want to add the option
 * @param name The name we want to add
 * @return The new set including the name
 */
static SET_NAME *add_set(SET_NAME *set, char *name) {
	if(!is_set_name(set,name)) {
		SET_NAME *curr;
		curr = (SET_NAME *)malloc(sizeof(SET_NAME));
		curr->name = (char *)name;
		curr->next  = (struct SET_NAME *)set;
		set=curr;
	}
	return(set);
}
