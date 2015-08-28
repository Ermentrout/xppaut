#ifndef _lunch_new_h_
#define _lunch_new_h_


#include <stdio.h>
#include "form_ode.h"



void file_inf(void);
void ps_write_pars(FILE *fp);
void do_info(FILE *fp);
int read_lunch(FILE *fp);
void do_lunch(int f);
void dump_eqn(FILE *fp);
void io_numerics(int f, FILE *fp);
void io_parameter_file(char *fn, int flag);
void io_ic_file(char *fn, int flag);
void io_parameters(int f, FILE *fp);
void io_exprs(int f, FILE *fp);
void io_graph(int f, FILE *fp);
void io_int(int *i, FILE *fp, int f, char *ss);
void io_double(double *z, FILE *fp, int f, char *ss);
void io_float(float *z, FILE *fp, int f, char *ss);
void io_string(char *s, int len, FILE *fp, int f);


#endif
