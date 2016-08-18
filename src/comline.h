#ifndef XPPAUT_COMLINE_H
#define XPPAUT_COMLINE_H

#include "xpplim.h"

/* --- Data --- */
extern char includefilename[MaxIncludeFiles][XPP_MAX_NAME];

extern int dryrun;
extern int newseed;
extern int NincludedFiles;
extern int Nintern_2_use;
extern int noicon;
extern int queryics;
extern int querypars;
extern int querysets;

/* --- Functions --- */
void do_comline(int argc, char **argv);
int if_needed_select_sets(void);
int if_needed_load_set(void);
int if_needed_load_par(void);
int if_needed_load_ic(void);
int if_needed_load_ext_options(void);

#endif /* XPPAUT_COMLINE_H */
