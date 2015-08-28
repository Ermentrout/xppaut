#ifndef _parserslow_h_
#define _parserslow_h_

#include "volterra.h"
#include "xpplim.h"

#define FUN1TYPE 9
#define FUN2TYPE 1
#define VARTYPE 3
#define CONTYPE 2
#define UFUNTYPE   24
#define SVARTYPE 4
#define SCONTYPE 32
#define NETTYPE 6
#define TABTYPE 7
#define USTACKTYPE 8
#define KERTYPE 10
#define MAXTYPE 20000000 


#define COM(a,b) ((a)*MAXTYPE+(b))


#define MAXARG 20
#define NEGATE 9
#define MINUS 4
#define LPAREN 0
#define RPAREN 1
#define COMMA  2
#define STARTTOK 10
#define ENDTOK 11

#define ENDEXP 999
#define ENDFUN 998
#define STARTDELAY 980
#define DELSYM  42
#define ENDDELAY 996
#define MYIF  995
#define MYELSE 993
#define MYTHEN 994
#define SUMSYM 990
#define ENDSUM 991
#define SHIFTSYM 64
#define ISHIFTSYM 67
#define ENDSHIFT 988
#define SUMINDEX 989
#define LASTTOK MAX_SYMBS-2
#define NUMSYM 987
#define NUMTOK 59
#define CONV 2
#define FIRST_ARG 73
#define ENDDELSHFT 986
#define DELSHFTSYM 65
#define ENDISHIFT 985
#define SETSYM  92
#define ENDSET 981
#define INDX 68
#define INDXVAR 984

/*#define STDSYM 95
*/
#define STDSYM 95

#define INDXCOM 922
#define STARTINDX 70
#define ENDINDX 69



/* #define MXLEN 32 */
#define MXLEN 10

typedef struct
        {
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



#define VECT_ROOT 500 




void init_rpn(void);
void free_ufuns(void);
int duplicate_name(char *junk);
int add_constant(char *junk);
int get_var_index(char *name);
int get_type(int index);
int add_con(char *name, double value);
int add_kernel(char *name, double mu, char *expr);
int add_var(char *junk, double value);
int add_expr(char *expr, int *command, int *length);
int add_net_name(int index, char *name);
int add_vect_name(int index, char *name);
int add_2d_table(char *name, char *file);
int add_file_table(int index, char *file);
int add_table_name(int index, char *name);
int add_form_table(int index, int nn, double xlo, double xhi, char *formula);
void set_old_arg_names(int narg);
void set_new_arg_names(int narg, char args[10][11]);
int add_ufun_name(char *name, int index, int narg);
void fixup_endfun(int *u, int l, int narg);
int add_ufun_new(int index, int narg, char *rhs, char args[20][11]);
int add_ufun(char *junk, char *expr, int narg);
int check_num(int *tok, double value);
int is_ufun(int x);
int is_ucon(int x);
int is_uvar(int x);
int isvar(int y);
int iscnst(int y);
int isker(int y);
int is_kernel(int x);
int is_lookup(int x);
int find_lookup(char *name);
void find_name(char *string, int *index);
int get_param_index(char *name);
int get_val(char *name, double *value);
int set_val(char *name, double value);
void set_ivar(int i, double value);
double get_ivar(int i);
int alg_to_rpn(int *toklist, int *command);
void pr_command(int *command);
void show_where(char *string, int index);
int function_sym(int token);
int unary_sym(int token);
int binary_sym(int token);
int pure_number(int token);
int gives_number(int token);
int check_syntax(int oldtoken, int newtoken);
int make_toks(char *dest, int *my_token);
void tokeninfo(int tok);
int do_num(char *source, char *num, double *value, int *ind);
void convert(char *source, char *dest);
void find_tok(char *source, int *index, int *tok);
double pmod(double x, double y);
void two_args(void);
double bessel_j(double x, double y);
double bessel_y(double x, double y);
double bessi(double nn, double x);
double bessi0(double x);
double bessi1(double x);
char *com_name(int com);
double do_shift(double shift, double variable);
double do_ishift(double shift, double variable);
double do_delay_shift(double delay, double shift, double variable);
double do_delay(double delay, double i);
void one_arg(void);
double normal(double mean, double std);
double max(double x, double y);
double min(double x, double y);

double neg(double z);
double recip(double z);
double heaviside(double z);
double rndom(double z);
double signum(double z);
double dnot(double x);
double dand(double x, double y);
double dor(double x, double y);
double dge(double x, double y);
double dle(double x, double y);
double deq(double x, double y);
double dne(double x, double y);
double dgt(double x, double y);
double dlt(double x, double y);
double evaluate(int *equat);
double eval_rpn(int *equat);

/*  STRING STUFF  */
#ifndef STRUPR
void strupr(char *s);
void strlwr(char *s);
#endif

/*****************************************************/






#endif



























