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
void set_auto_eval_flags(int f);
void set_table_name(char *name, int index);
void view_table(int index);
void new_lookup_com(int i);
void new_lookup_ok(void);
double lookupxy(double x, int n, double *xv, double *yv);
double tab_interp(double xlo, double h, double x, double *y, int n, int i);
double lookup(double x, int index);
void init_table(void);
void redo_all_fun_tables(void);
int eval_fun_table(int n, double xlo, double xhi, char *formula, double *y);
int create_fun_table(int npts, double xlo, double xhi, char *formula, int index);
int load_table(char *filename, int index);
int get_lookup_len(int i);

#endif /* XPPAUT_TABULAR_H */

