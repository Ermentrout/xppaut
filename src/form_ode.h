#ifndef XPPAUT_FORM_ODE_H
#define XPPAUT_FORM_ODE_H

#include <stdio.h>
#include "newpars.h"
#include "util/strutil.h"
#include "xpplim.h"

/* --- Macros --- */
#define MAXCOMMENTS 500
#define MAXLINES 5000
#define MAXONLY 1000
#define MAXVNAM 33

/* --- Types --- */
typedef struct {
	char *text,*action;
	int aflag;
} ACTION;

typedef struct {
	int *com;
	char *string;
	char *name;
	int side;
} BC_STRUCT;

typedef struct {
	char *name,*value;
} FIXINFO;

/* --- Data --- */
extern char *save_eqn[MAXLINES];
extern char upar_names[MAXPAR][11];
extern char uvar_names[MAXODE][12];

extern int ConvertStyle;
extern int EqType[MAXODE];
extern int FIX_VAR;
extern int N_plist;
extern int n_comments;
extern int NCON;
extern int NCON_START;
extern int NDELAYS;
extern int NEQ_MIN;
extern int NLINES;
extern int NODE;
extern int NSYM_START;
extern int NUPAR;
extern int PrimeStart;

extern int *my_ode[MAXODE];
extern int *plotlist;

extern double default_ic[MAXODE];
extern double default_val[MAXPAR];

extern ACTION comments[MAXCOMMENTS];
extern FILE *convertf;
extern FIXINFO fixinfo[MAXODE];
extern BC_STRUCT my_bc[MAXODE];

/* --- Functions --- */
void create_plot_list(void);
int disc(char *string);
int find_char(char *s1, char *s2, int i0, int *i1);
int get_eqn(FILE *fptr);
int make_eqn(void);
int search_array(char *old, char *new_char, int *i1, int *i2, int *flag);
void strip_saveqn(void);
void subsk(char *big, char *new_char, int k, int flag);

#endif /* XPPAUT_FORM_ODE_H */
