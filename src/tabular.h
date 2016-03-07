#ifndef XPPAUT_TABULAR_H
#define XPPAUT_TABULAR_H

#include "xpplim.h"

/* --- Types --- */
/* flag=0 if virgin array,
 * flag=1 if already allocated;
 * flag=2 for function
 *
 * interp=0 for normal interpolation,
 * interp=1 for 'step'
 * interp=2 for cubic spline table
 *
 * xyvals=1 if both x and y vals are needed
 * (xyvals=0 is faster lookup )
 */
typedef struct {
	double xlo,xhi,dx;
	double *y,*x;
	int n,flag,interp,autoeval;
	int xyvals;
	char filename[128],name[12];
} TABULAR;

/* --- Data --- */
extern TABULAR my_table[MAX_TAB];

/* --- Functions --- */
int create_fun_table(int npts, double xlo, double xhi, char *formula, int index);
int get_lookup_len(int i);
void init_table(void);
int load_table(char *filename, int index);
double lookup(double x, int index);
void new_lookup_com(int i);
void redo_all_fun_tables(void);
void set_auto_eval_flags(int f);
void set_table_name(char *name, int index);

#endif /* XPPAUT_TABULAR_H */
