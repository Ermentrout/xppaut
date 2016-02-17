#ifndef XPP_PARSERSLOW_H
#define XPP_PARSERSLOW_H

#include "volterra.h"
#include "xpplim.h"

/* --- Macros --- */
#define MAXARG 20
#define FUN1TYPE 9
#define FUN2TYPE 1
#define VARTYPE 3  /* standard variable */
#define CONTYPE 2  /* standard parameter */
#define UFUNTYPE   24
#define SVARTYPE 4  /* shifted variable */
#define SCONTYPE 32  /* shifted constant  */
#define NETTYPE 6
#define TABTYPE 7
#define USTACKTYPE 8
#define KERTYPE 10
#define VECTYPE 13  /* for vectorized stuff */
#define EVECTYPE 14 /* treat vector like a function */
#define MAXTYPE 20000000  /* this is the maximum number of named stuff */

/* --- Types --- */
typedef struct {
  int narg;
  char args[MAXARG][11];
} UFUN_ARG;

/* --- Data --- */
extern double constants[MAXPAR];
extern int ERROUT;
extern int narg_fun[MAXUFUN];
extern int NFUN;
extern int NKernel;
extern int NSYM;
extern int NTable;
extern int NVAR;
extern int *ufun[MAXUFUN];
extern UFUN_ARG ufun_arg[MAXUFUN];
extern char *ode_names[MAXODE];
extern char *ufun_def[MAXUFUN];
extern char ufun_names[MAXUFUN][12];
extern double variables[MAXODE1];

/* --- Functions --- */
void init_rpn(void);
int get_var_index(char *name);
int add_con(char *name, double value);
int add_kernel(char *name, double mu, char *expr);
int add_var(char *junk, double value);
int add_expr(char *expr, int *command, int *length);
int add_net_name(int index, char *name);
int add_vector_name(int index, char *name);
int add_2d_table(char *name, char *file);
int add_file_table(int index, char *file);
int add_table_name(int index, char *name);
int add_form_table(int index, int nn, double xlo, double xhi, char *formula);
void set_old_arg_names(int narg);
void set_new_arg_names(int narg, char args[10][11]);
int add_ufun_name(char *name, int index, int narg);
void fixup_endfun(int *u, int l, int narg);
int add_ufun_new(int index, int narg, char *rhs, char args[MAXARG][11]);
int add_ufun(char *junk, char *expr, int narg);
int find_lookup(char *name);
int get_param_index(char *name);
int get_val(char *name, double *value);
int set_val(char *name, double value);
void set_ivar(int i, double value);
double get_ivar(int i);
int do_num(char *source, char *num, double *value, int *ind);
void convert(char *source, char *dest);
#endif /* XPP_PARSERSLOW_H */
