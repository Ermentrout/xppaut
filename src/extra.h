#ifndef XPPAUT_EXTRA_H
#define XPPAUT_EXTRA_H

/* --- Data --- */
extern int dll_flag;
extern char dll_fun[256];
extern char dll_lib[256];

/* --- Functions --- */
void load_new_dll(void);
int my_fun(double *in, double *out, int nin, int nout, double *v, double *c);
void auto_load_dll(void);
void do_in_out(void);
void add_export_list(char *in, char *out);
void check_inout(void);
int get_export_count(char *s);
void do_export_list(void);
void parse_inout(char *l, int flag);

#endif /* XPPAUT_EXTRA_H */
