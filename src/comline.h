#ifndef XPPAUT_COMLINE_H
#define XPPAUT_COMLINE_H

#include "ggets.h"
#include "load_eqn.h"

/* --- Types --- */
typedef struct {
  char *name;
  struct SET_NAME * next;
} SET_NAME;

/* --- Data --- */
extern char includefilename[MaxIncludeFiles][XPP_MAX_NAME];

extern int dryrun;
extern int loadincludefile;
extern int NincludedFiles;
extern int Nintern_2_use;
extern int noicon;
extern int queryics;
extern int querypars;
extern int querysets;

/* --- Functions --- */
int is_set_name(SET_NAME *set, char *nam);
SET_NAME *add_set(SET_NAME *set, char *nam);
SET_NAME *rm_set(SET_NAME *set, char *nam);
void do_comline(int argc, char **argv);
int if_needed_select_sets(void);
int if_needed_load_set(void);
int if_needed_load_par(void);
int if_needed_load_ic(void);
int if_needed_load_ext_options(void);
int parse_it(char *com);

#endif /* XPPAUT_COMLINE_H */
