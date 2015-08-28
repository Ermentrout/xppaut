
#ifndef _tabular_h_
#define _tabular_h_


/* tabular.c */
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


#endif

