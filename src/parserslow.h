#ifndef XPPAUT_PARSERSLOW_H
#define XPPAUT_PARSERSLOW_H

#include "newpars.h"
#include "xpplim.h"

/* --- Macros --- */
#define KN_OLD 1
#define KN 0

/* KERNEL.flag values */
#define CONV 2

#define MXLEN 10
#define VECTYPE 13  /* for vectorized stuff */

#define SETVAR(i,x) if((i)<NVAR)variables[(i)]=(x);
#define GETVAR(i) (i)<NVAR ? variables[(i)]:0.0

/* --- Types --- */
typedef struct {
  double k_n1,k_n,sum,betnn,mu,*al,*cnv;
  int *formula,flag,*kerform;
  char name[20],*expr,*kerexpr;
} KERNEL;

typedef struct {
	char name[MXLEN+1];
	int len;
	int com;
	int arg;
	int pri;
} SYMBOL;

typedef struct {
  int narg;
  char args[MAXARG][11];
} UFUN_ARG;

/* --- Data --- */
extern char *ode_names[MAXODE];
extern char *ufun_def[MAXUFUN];
extern char ufun_names[MAXUFUN][12];

extern double constants[MAXPAR];
extern double *Memory[MAXKER];
extern double variables[MAXODE1];

extern int ERROUT;
extern int narg_fun[MAXUFUN];
extern int MaxPoints;
extern int NFUN;
extern int NKernel;
extern int NSYM;
extern int NTable;
extern int NVAR;
extern int RandSeed;
extern int *ufun[MAXUFUN];

extern KERNEL kernel[MAXKER];
extern UFUN_ARG ufun_arg[MAXUFUN];

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
#endif /* XPPAUT_PARSERSLOW_H */
