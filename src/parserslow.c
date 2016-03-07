#include "parserslow.h"

#ifndef HAVE_WCTYPE_H
# include <ctype.h>
#else
# include <wctype.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "comline.h"
#include "delay_handle.h"
#include "ggets.h"
#include "markov.h"
#include "simplenet.h"
#include "strutil.h"
#include "tabular.h"
#include "volterra2.h"
#include "xpplim.h"

/* --- Macros --- */
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

#define NEGATE 9
#define MINUS 4
#define LPAREN 0
#define RPAREN 1
#define COMMA  2
#define STARTTOK 10
#define ENDTOK 11

#define ENDEXP 999
#define ENDFUN 998
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
#define LASTTOK MAX_SYMBS-2
#define NUMSYM 987
#define NUMTOK 59
#define FIRST_ARG 73
#define ENDDELSHFT 986
#define DELSHFTSYM 65
#define ENDISHIFT 985
#define ENDSET 981
#define INDX 68
#define STDSYM 95

#define INDXCOM 922

#define THOUS 10000
#define DOUB_EPS 2.23E-15
#define POP stack[--stack_pointer]
#define PUSH(a) zippy=(a); stack[stack_pointer++]=zippy;

#ifdef NOLGAMMA
double lgamma();
#endif

#ifndef M_PI
# define M_PI	3.14159265358979323846264338327950288
#endif

/* --- Forward declarations --- */
static int duplicate_name(char *junk);
static int add_constant(char *junk);
static int is_ufun(int x);
static int is_ucon(int x);
static int is_uvar(int x);
static int isvar(int y);
static int iscnst(int y);
static int isker(int y);
static int is_lookup(int x);
static void find_name(char *string, int *index);
static int alg_to_rpn(int *toklist, int *command);
static void show_where(char *string, int index);
static int function_sym(int token);
static int unary_sym(int token);
static int binary_sym(int token);
static int pure_number(int token);
static int gives_number(int token);
static int check_syntax(int oldtoken, int newtoken);
static int make_toks(char *dest, int *my_token);
static void tokeninfo(int tok);
static void find_tok(char *source, int *index, int *tok);
static double pmod(double x, double y);

static double bessel_j(double x, double y);
static double bessel_y(double x, double y);
static double bessi(double nn, double x);
static double bessi0(double x);
static double bessi1(double x);
static double do_shift(double shift, double variable);
static double do_ishift(double shift, double variable);
static double do_delay_shift(double delay, double shift, double variable);
static double do_delay(double delay, double i);
static void one_arg(void);
static double max(double x, double y);
static double min(double x, double y);

static double neg(double z);
static double recip(double z);
static double heaviside(double z);
static double rndom(double z);
static double signum(double z);

static double dnot(double x);
static double dand(double x, double y);
static double dor(double x, double y);
static double dge(double x, double y);
static double dle(double x, double y);
static double deq(double x, double y);
static double dne(double x, double y);
static double dgt(double x, double y);
static double dlt(double x, double y);

static void two_args(void);
static double eval_rpn(/* int* */ );


/* --- Data --- */
static double BoxMuller;
static int BoxMullerFlag=0;
static double CurrentIndex=0;
static int SumIndex=1;
static int stack_pointer,uptr;
static double stack[200], ustack[200];
static double zippy;

int RandSeed=12345678;
int ERROUT;
int MaxPoints;
double *Memory[MAXKER];
int NDELAYS=0;
int NKernel;
int NTable;

double constants[MAXPAR];
double variables[MAXODE1];
int *ufun[MAXUFUN];
char *ufun_def[MAXUFUN];
char ufun_names[MAXUFUN][12];
int narg_fun[MAXUFUN];


KERNEL kernel[MAXKER];
UFUN_ARG ufun_arg[MAXUFUN];
static SYMBOL my_symb[MAX_SYMBS] = {
	{"(",1,999,0,1},      /*  0   */
	{")",1,999,0,2},
	{",",1,999,0,3},
	{"+",1,COM(FUN2TYPE,0),0,4},
	{"-",1,COM(FUN2TYPE,1),0,4},
	{"*",1,COM(FUN2TYPE,2),0,6},
	{"/",1,COM(FUN2TYPE,3),0,6},
	{"^",1,COM(FUN2TYPE,5),0,7},
	{"**",2,COM(FUN2TYPE,5),0,7},
	{"~",1,COM(FUN1TYPE,14),0,6},
	{"START",5,-1,0,0},  /* 10  */
	{"END",3,999,0,-1},
	{"ATAN2",5,COM(FUN2TYPE,4),2,10},
	{"MAX",3,COM(FUN2TYPE,6),2,10},
	{"MIN",3,COM(FUN2TYPE,7),2,10},
	{"SIN",3,COM(FUN1TYPE,0),0,10},
	{"COS",3,COM(FUN1TYPE,1),0,10},
	{"TAN",3,COM(FUN1TYPE,2),0,10},
	{"ASIN",4,COM(FUN1TYPE,3),0,10},
	{"ACOS",4,COM(FUN1TYPE,4),0,10},
	{"ATAN",4,COM(FUN1TYPE,5),0,10},  /* 20  */
	{"SINH",4,COM(FUN1TYPE,6),0,10},
	{"TANH",4,COM(FUN1TYPE,7),0,10},
	{"COSH",4,COM(FUN1TYPE,8),0,10},
	{"ABS",3,COM(FUN1TYPE,9),0,10},
	{"EXP",3,COM(FUN1TYPE,10),0,10},
	{"LN",2,COM(FUN1TYPE,11),0,10},
	{"LOG",3,COM(FUN1TYPE,11),0,10},
	{"LOG10",5,COM(FUN1TYPE,12),0,10},
	{"SQRT",4,COM(FUN1TYPE,13),0,10},
	{"HEAV",4,COM(FUN1TYPE,16),0,10},  /*  30 */
	{"SIGN",4,COM(FUN1TYPE,17),0,10},
	{"#$%1",4,COM(USTACKTYPE,0),0,10},
	{"#$%2",4,COM(USTACKTYPE,1),0,10},
	{"#$%3",4,COM(USTACKTYPE,2),0,10},
	{"#$%4",4,COM(USTACKTYPE,3),0,10},
	{"#$%5",4,COM(USTACKTYPE,4),0,10},
	{"#$%6",4,COM(USTACKTYPE,5),0,10},
	{"#$%7",4,COM(USTACKTYPE,6),0,10},
	{"#$%8",4,COM(USTACKTYPE,7),0,10},
	{"FLR",3,COM(FUN1TYPE,18),0,10},  /*  40 */
	{"MOD",3,COM(FUN2TYPE,8),2,10}, /*  41 */
	{"DELAY",5,ENDDELAY,2,10},      /*  42 */   /*  Delay symbol */
	{"RAN",3,COM(FUN1TYPE,19),1,10}, /* 43 */
	{"&",1,COM(FUN2TYPE,9),0,6},  /* logical stuff  */
	{"|",1,COM(FUN2TYPE,10),0,4},
	{">",1,COM(FUN2TYPE,11),0,7},
	{"<",1,COM(FUN2TYPE,12),0,7},
	{"==",2,COM(FUN2TYPE,13),0,7},
	{">=",2,COM(FUN2TYPE,14),0,7},
	{"<=",2,COM(FUN2TYPE,15),0,7}, /*50 */
	{"IF",2,995,1,10},
	{"THEN",4,994,1,10},
	{"ELSE",4,993,1,10},
	{"!=",2,COM(FUN2TYPE,16),0,7},
	{"NOT",3,COM(FUN1TYPE,20),0,6},
	{"NORMAL",6,COM(FUN2TYPE,17),2,10}, /* returns normally dist number */
	{"BESSELJ",7,COM(FUN2TYPE,18),2,10}, /* Bessel J   */
	{"BESSELY",7,COM(FUN2TYPE,19),2,10}, /* Bessel Y */
	{"NXXQQ",5,NUMSYM,0,10},
	{"ERF", 3, COM(FUN1TYPE,21),0,10}, /* 60 */
	{"ERFC",4,COM(FUN1TYPE,22),0,10},
	{"SUM",3,SUMSYM,2,10},
	{"OF",2,ENDSUM,0,10},
	{"SHIFT",5,ENDSHIFT,2,10},
	{"DEL_SHFT",8,ENDDELSHFT,3,10},/* 65 */
	{"HOM_BCS",7,COM(FUN1TYPE,23),0,10},
	{"ISHIFT",6,ENDISHIFT,2,10}, /* 67 */
	{"@",1,INDXCOM,0,10}, /*68 */
	{"]",1,ENDSHIFT,0,10},
	{"[",1,ENDSHIFT,0,10}, /*70 */
	{"POISSON",7,COM(FUN1TYPE,24),0,10}, /* 71 */
	{"SET",3,ENDSET,3,10}, /* 72 */
	{"ARG1",4,COM(USTACKTYPE,0),0,10}, /*  FIXXX ????  */
	{"ARG2",4,COM(USTACKTYPE,1),0,10},
	{"ARG3",4,COM(USTACKTYPE,2),0,10},
	{"ARG4",4,COM(USTACKTYPE,3),0,10},
	{"ARG5",4,COM(USTACKTYPE,4),0,10},
	{"ARG6",4,COM(USTACKTYPE,5),0,10},
	{"ARG7",4,COM(USTACKTYPE,6),0,10},
	{"ARG8",4,COM(USTACKTYPE,7),0,10},
	{"ARG9",4,COM(USTACKTYPE,8),0,10},
	{"ARG10",5,COM(USTACKTYPE,9),0,10},
	{"ARG11",5,COM(USTACKTYPE,10),0,10},
	{"ARG12",5,COM(USTACKTYPE,11),0,10},
	{"ARG13",5,COM(USTACKTYPE,12),0,10},
	{"ARG14",5,COM(USTACKTYPE,13),0,10},
	{"ARG15",5,COM(USTACKTYPE,14),0,10},
	{"ARG16",5,COM(USTACKTYPE,15),0,10},
	{"ARG17",5,COM(USTACKTYPE,16),0,10},
	{"ARG18",5,COM(USTACKTYPE,17),0,10},
	{"ARG19",5,COM(USTACKTYPE,18),0,10},
	{"ARG20",5,COM(USTACKTYPE,19),0,10},
	{"BESSELI",7,COM(FUN2TYPE,20),2,10},/* Bessel I  # 93 */
	{"LGAMMA",6,COM(FUN1TYPE,25),1,10} /* Log Gamma  #94 */
};

int NCON=0,NVAR=0,NFUN=0;
int NSYM=STDSYM;

/*     pointers to functions    */
double (*fun1[50])(/* double */ );
double (*fun2[50])(/* double,double */ );


/*************************
  RPN COMPILER           *
**************************/



/*****************************
*      PARSER.C              *
*
*
*     parses any algebraic expression
*     and converts to an integer array
*     to be interpreted by the rpe_val
*     function.
*
*     the main data structure is a contiguous
*     list of symbols with their priorities
*     and their symbol value
*
*     on the first pass, the expression is converted to
*     a list of integers without any checking except for
*     valid symbols and for numbers
*
*     next this list of integers is converted to an RPN expression
*     for evaluation.
*
*
*  6/95  stuff added to add names to namelist without compilation
*************************************************************/
void init_rpn(void) {
	ERROUT = 1;
	NCON = 0;
	NFUN = 0;
	NVAR = 0;
	NKernel=0;

	MaxPoints=4000;
	NSYM = STDSYM;
	two_args();
	one_arg();
	add_con("PI", M_PI);

	add_con("I'",0.0);
	/*   This is going to be for interacting with the animator */
	SumIndex=NCON-1;
	add_con("mouse_x",0.0);
	add_con("mouse_y",0.0);
	add_con("mouse_vx",0.0);
	add_con("mouse_vy",0.0);

	init_table();
	if (newseed==1) {
		RandSeed=time(0);
	}
	nsrand48(RandSeed);
}


int duplicate_name(char *junk) {
	int i;
	find_name(junk,&i);
	if(i>=0) {
		if(ERROUT) {
			printf("%s is a duplicate name\n",junk);
		}
		return(1);
	}
	return(0);
}


int add_constant(char *junk) {
	int len;
	char string[100];
	if(duplicate_name(junk)==1) {
		return(1);
	}
	if(NCON>=MAXPAR) {
		if(ERROUT) {
			printf("too many constants !!\n");
		}
		return(1);
	}
	convert(junk,string);
	len=strlen(string);
	if(len<1) {
		plintf("Empty parameter - remove spaces\n");
		return 1;
	}
	if(len>MXLEN) {
		len=MXLEN;
	}
	strncpy(my_symb[NSYM].name,string,len);
	my_symb[NSYM].name[len]='\0';
	my_symb[NSYM].len=len;
	my_symb[NSYM].pri=10;
	my_symb[NSYM].arg=0;
	my_symb[NSYM].com=COM(CONTYPE,NCON-1);
	NSYM++;
	return(0);
}


int get_var_index(char *name) {
	int type,com;
	find_name(name,&type);
	if(type<0) {
		return -1;
	}
	com=my_symb[type].com;
	if(is_uvar(com)) {
		return(com%MAXTYPE);
	}
	return(-1);
}


int add_con(char *name, double value) {
	if(NCON>=MAXPAR) {
		if(ERROUT) {
			printf("too many constants !!\n");
		}
		return(1);
	}
	constants[NCON]=value;
	NCON++;
	return(add_constant(name));
}


int add_kernel(char *name, double mu, char *expr) {
	char string[100];
	int len,i,in=-1;
	if(duplicate_name(name)==1) {
		return(1);
	}
	if(NKernel==MAXKER) {
		plintf("Too many kernels..\n");
		return(1);
	}
	if(mu<0 || mu>=1.0) {
		plintf(" mu must lie in [0,1.0) \n");
		return(1);
	}
	convert(name,string);
	len=strlen(string);
	if(len>MXLEN) {
		len=MXLEN;
	}
	strncpy(my_symb[NSYM].name,string,len);
	my_symb[NSYM].name[len]='\0';
	my_symb[NSYM].len=len;
	my_symb[NSYM].pri=10;
	my_symb[NSYM].arg=0;
	my_symb[NSYM].com=COM(KERTYPE,NKernel);
	kernel[NKernel].k_n1=0.0;
	kernel[NKernel].mu=mu;
	kernel[NKernel].k_n=0.0;
	kernel[NKernel].k_n1=0.0;
	kernel[NKernel].flag=0;
	for(i=0;i<strlen(expr);i++) {
		if(expr[i]=='#') {
			in=i;
		}
	}
	if(in==0 || in==(strlen(expr)-1)) {
		plintf("Illegal use of convolution...\n");
		return(1);
	}
	if(in>0) {
		kernel[NKernel].flag=CONV;
		kernel[NKernel].expr=(char *)malloc(strlen(expr)+2-in);
		kernel[NKernel].kerexpr=(char *)malloc(in+1);
		for(i=0;i<in;i++) {
			kernel[NKernel].kerexpr[i]=expr[i];
		}
		kernel[NKernel].kerexpr[in]=0;
		for(i=in+1;i<strlen(expr);i++) {
			kernel[NKernel].expr[i-in-1]=expr[i];
		}
		kernel[NKernel].expr[strlen(expr)-in-1]=0;
		plintf("Convolving %s with %s\n",
			   kernel[NKernel].kerexpr,kernel[NKernel].expr);
	} else {
		kernel[NKernel].expr=(char *)malloc(strlen(expr)+2);
		strcpy(kernel[NKernel].expr,expr);
	}
	NSYM++;
	NKernel++;
	return(0);
}


int add_var(char *junk, double value) {
	char string[100];
	int len;
	if(duplicate_name(junk)==1) {
		return(1);
	}
	if(NVAR>=MAXODE1) {
		if(ERROUT) {
			printf("too many variables !!\n");
		}
		return(1);
	}
	convert(junk,string);
	len=strlen(string);
	if(len>MXLEN) {
		len=MXLEN;
	}
	strncpy(my_symb[NSYM].name,string,len);
	my_symb[NSYM].name[len]='\0';
	my_symb[NSYM].len=len;
	my_symb[NSYM].pri=10;
	my_symb[NSYM].arg=0;
	my_symb[NSYM].com=COM(VARTYPE,NVAR);
	NSYM++;
	variables[NVAR]=value;
	NVAR++;
	return(0);
}


int add_expr(char *expr, int *command, int *length) {
	char dest[1024];
	int my_token[1024];
	int err,i;
	convert(expr,dest);
	err=make_toks(dest,my_token);
	if(err!=0) {
		return(1);
	}
	err = alg_to_rpn(my_token,command);
	if(err!=0) {
		return(1);
	}
	i=0;
	while(command[i]!=ENDEXP) {
		i++;
	}
	*length=i+1;
	return(0);
}


int add_vector_name(int index,char *name) {
	char string[50];
	int len=strlen(name);
	plintf(" Adding vectorizer %s %d \n",name,index);
	if(duplicate_name(name)==1) {
		return(1);
	}
	convert(name,string);
	printf(" 1\n");
	if(len>MXLEN) {
		len=MXLEN;
	}
	strncpy(my_symb[NSYM].name,string,len);
	my_symb[NSYM].name[len]='\0';
	my_symb[NSYM].len=len;
	my_symb[NSYM].pri=10;
	my_symb[NSYM].arg=1;
	my_symb[NSYM].com=COM(VECTYPE,index);
	NSYM++;
	return(0);
}


int add_net_name(int index, char *name) {
	char string[50];
	int len=strlen(name);
	plintf(" Adding net %s %d \n",name,index);
	if(duplicate_name(name)==1) {
		return(1);
	}
	convert(name,string);
	if(len>MXLEN) {
		len=MXLEN;
	}
	strncpy(my_symb[NSYM].name,string,len);
	my_symb[NSYM].name[len]='\0';
	my_symb[NSYM].len=len;
	my_symb[NSYM].pri=10;
	my_symb[NSYM].arg=1;
	my_symb[NSYM].com=COM(NETTYPE,index);
	NSYM++;
	return(0);
}


int add_2d_table(char *name, char *file) {
	plintf(" TWO D NOT HERE YET \n");
	return(1);
}


int add_file_table(int index, char *file) {
	char file2[1000];
	int i2=0,i1=0,n;
	char ch;
	n=strlen(file);
	for(i1=0;i1<n;i1++) {
		ch=file[i1];
		if(((int)ch>31) && ((int)ch<127)) {
			file2[i2]=ch;
			i2++;
		}
	}
	file2[i2]=0;
	if(load_table(file2,index)==0) {
		if(ERROUT) {
			printf("Problem with creating table !!\n");
		}
		return(1);
	}
	return(0);
}


int add_table_name(int index, char *name) {
	char string[50];
	int len=strlen(name);
	if(duplicate_name(name)==1) {
		return(1);
	}
	convert(name,string);
	if(len>MXLEN) {
		len=MXLEN;
	}
	strncpy(my_symb[NSYM].name,string,len);
	my_symb[NSYM].name[len]='\0';
	my_symb[NSYM].len=len;
	my_symb[NSYM].pri=10;
	my_symb[NSYM].arg=1;
	my_symb[NSYM].com=COM(TABTYPE, index);
	set_table_name(name,index);
	NSYM++;
	return(0);
}


/* ADD LOOKUP TABLE   */
int add_form_table(int index, int nn, double xlo, double xhi, char *formula) {
	if(create_fun_table(nn,xlo,xhi,formula,index)==0) {
		if(ERROUT) {
			printf("Problem with creating table !!\n");
		}
		return(1);
	}
	return(0);
}


void set_old_arg_names(int narg) {
	int i;
	for(i=0;i<narg;i++) {
		sprintf(my_symb[FIRST_ARG+i].name,"ARG%d",i+1);
		my_symb[FIRST_ARG+i].len=4;
	}
}


void set_new_arg_names(int narg, char args[10][11]) {
	int i;
	for(i=0;i<narg;i++) {
		strcpy(my_symb[FIRST_ARG+i].name,args[i]);
		my_symb[FIRST_ARG+i].len=strlen(args[i]);
	}
}


/* NEW ADD_FUN for new form_ode code  */
int add_ufun_name(char *name, int index, int narg) {
	char string[50];
	int len=strlen(name);
	if(duplicate_name(name)==1) {
		return(1);
	}
	if(index>=MAXUFUN) {
		if(ERROUT) {
			printf("too many functions !!\n");
		}
		return(1);
	}
	plintf(" Added user fun %s \n",name);
	convert(name,string);
	if(len>MXLEN) {
		len=MXLEN;
	}
	strncpy(my_symb[NSYM].name,string,len);
	my_symb[NSYM].name[len]='\0';
	my_symb[NSYM].len=len;
	my_symb[NSYM].pri=10;
	my_symb[NSYM].arg=narg;
	my_symb[NSYM].com=COM(UFUNTYPE, index);
	NSYM++;
	strcpy(ufun_names[index],name);
	return (0);
}


void fixup_endfun(int *u, int l, int narg) {
	u[l-1]=ENDFUN;
	u[l]=narg;
	u[l+1]=ENDEXP;
}


int add_ufun_new(int index, int narg, char *rhs, char args[MAXARG][11]) {
	int i,l;
	int end;
	if(narg>MAXARG) {
		plintf("Maximal arguments exceeded \n");
		return(1);
	}
	ufun[index]=(int *)malloc(1024);
	if(ufun[index]==NULL) {
		if(ERROUT) {
			printf("not enough memory!!\n");
		}
		return(1);
	}
	ufun_def[index]=(char *)malloc(MAXEXPLEN);
	if(ufun_def[index]==NULL) {
		if(ERROUT) {
			printf("not enough memory!!\n");
		}
		return(1);
	}
	ufun_arg[index].narg=narg;
	for(i=0;i<narg;i++) {
		strcpy(ufun_arg[index].args[i],args[i]);
	}
	set_new_arg_names(narg,args);
	if(add_expr(rhs,ufun[index],&end)==0) {
		ufun[index][end-1]=ENDFUN;
		ufun[index][end]=narg;
		ufun[index][end+1]=ENDEXP;
		strcpy(ufun_def[index],rhs);
		l=strlen(ufun_def[index]);
		ufun_def[index][l]=0;
		narg_fun[index]=narg;
		set_old_arg_names(narg);
		return(0);
	}
	set_old_arg_names(narg);
	if(ERROUT) {
		printf(" ERROR IN FUNCTION DEFINITION\n");
	}
	return(1);
}


/* ADD_UFUN   */
int add_ufun(char *junk, char *expr, int narg) {
	char string[50];
	int i,l;
	int end;
	int len=strlen(junk);

	if(duplicate_name(junk)==1) {
		return(1);
	}
	if(NFUN>=MAXUFUN) {
		if(ERROUT) {
			printf("too many functions !!\n");
		}
		return(1);
	}
	ufun[NFUN]=(int *)malloc(1024);
	if(ufun[NFUN]==NULL) {
		if(ERROUT) {
			printf("not enough memory!!\n");
		}
		return(1);
	}
	ufun_def[NFUN]=(char *)malloc(MAXEXPLEN);
	if(ufun_def[NFUN]==NULL) {
		if(ERROUT) {
			printf("not enough memory!!\n");
		}
		return(1);
	}

	convert(junk,string);
	if(add_expr(expr,ufun[NFUN],&end)==0) {
		if(len>MXLEN) {
			len=MXLEN;
		}
		strncpy(my_symb[NSYM].name,string,len);
		my_symb[NSYM].name[len]='\0';
		my_symb[NSYM].len=len;
		my_symb[NSYM].pri=10;
		my_symb[NSYM].arg=narg;
		my_symb[NSYM].com=COM(UFUNTYPE, NFUN);
		NSYM++;
		ufun[NFUN][end-1]=ENDFUN;
		ufun[NFUN][end]=narg;
		ufun[NFUN][end+1]=ENDEXP;
		strcpy(ufun_def[NFUN],expr);
		l=strlen(ufun_def[NFUN]);
		ufun_def[NFUN][l-1]=0;
		strcpy(ufun_names[NFUN],junk);
		narg_fun[NFUN]=narg;
		for(i=0;i<narg;i++) {
			sprintf(ufun_arg[NFUN].args[i],"ARG%d",i+1);
		}
		NFUN++;
		return(0);
	}
	if(ERROUT) {
		printf(" ERROR IN FUNCTION DEFINITION\n");
	}
	return(1);
}


int is_ufun(int x) {
	if((x/MAXTYPE)==UFUNTYPE) {
		return(1);
	}
	return(0);
}


int is_ucon(int x) {
	if (x/MAXTYPE == CONTYPE) {
		return(1);
	}
	return(0);
}


int is_uvar(int x) {
	if (x/MAXTYPE==VARTYPE) {
		return(1);
	}
	return(0);
}


int isvar(int y) {
	return (y == VARTYPE);
}


int iscnst(int y) {
	return (y == CONTYPE);
}


int isker(int y) {
	return (y == KERTYPE);
}


int is_lookup(int x) {
	if((x/MAXTYPE)==TABTYPE) {
		return(1);
	}
	return(0);
}


int find_lookup(char *name) {
	int index,com;
	find_name(name,&index);
	if(index==-1) {
		return(-1);
	}
	com=my_symb[index].com;
	if(is_lookup(com)) {
		return(com%MAXTYPE);
	}
	return(-1);
}


void find_name(char *string, int *index) {
	char junk[100];
	int i,len;
	convert(string,junk);
	len=strlen(junk);
	for(i=0;i<NSYM;i++) {
		if(len==my_symb[i].len) {
			if(strncmp(my_symb[i].name,junk,len)==0) {
				break;
			}
		}
	}
	if(i<NSYM) {
		*index=i;
	} else {
		*index=-1;
	}
}


int get_param_index(char *name) {
	int type,com;
	find_name(name,&type);
	if(type<0) {
		return(-1);
	}
	com=my_symb[type].com;
	if(is_ucon(com)) {
		return(com % MAXTYPE);
	}
	return(-1);
}


int get_val(char *name, double *value) {
	int type,com;
	*value=0.0;
	find_name(name,&type);
	if(type<0) {
		return(0);
	}
	com=my_symb[type].com;
	if(is_ucon(com)) {
		*value=constants[com % MAXTYPE];
		return(1);
	}
	if(is_uvar(com)) {
		*value=variables[com % MAXTYPE];
		return(1);
	}
	return(0);
}

int set_val(char *name, double value) {
	int type,com;
	find_name(name,&type);
	if(type<0) {
		return(0);
	}
	com=my_symb[type].com;
	if(is_ucon(com)) {
		constants[com % MAXTYPE]=value;
		return(1);
	}
	if(is_uvar(com)) {
		variables[com % MAXTYPE]=value;
		return(1);
	}
	return(0);
}


void set_ivar(int i, double value) {
	SETVAR(i,value);
}


double get_ivar(int i) {
	return(GETVAR(i));
}


int alg_to_rpn(int *toklist, int *command) {
	int tokstak[500],comptr=0,tokptr=0,lstptr=0,temp;
	int ncomma=0;
	int loopstk[100];
	int lptr=0;
	int nif=0,nthen=0,nelse=0;
	int newtok,oldtok;
	int my_com,my_arg,jmp;

	tokstak[0]=STARTTOK;
	tokptr=1;
	oldtok=STARTTOK;
	while(1) {
getnew:
		newtok=toklist[lstptr++];
		if(newtok==DELSYM) {
			temp=my_symb[toklist[lstptr+1]].com;
			/* !! */
			if(is_uvar(temp)) {
				/* ram -- is this right? not sure I understand what was happening here */
				my_symb[LASTTOK].com=COM(SVARTYPE,temp%MAXTYPE); /* create a temporary sybol */
				NDELAYS++;
				toklist[lstptr+1]=LASTTOK;
				my_symb[LASTTOK].pri=10;
			} else {
				printf("Illegal use of DELAY \n");
				return(1);
			}
		}
		/*        check for delshft symbol             */
		if(newtok==DELSHFTSYM) {
			temp=my_symb[toklist[lstptr+1]].com;
			/* !! */
			if(is_uvar(temp)) {
				/* ram -- same issue */
				my_symb[LASTTOK].com=COM(SVARTYPE, temp%MAXTYPE); /* create a temporary sybol */
				NDELAYS++;
				toklist[lstptr+1]=LASTTOK;
				my_symb[LASTTOK].pri=10;
			} else {
				printf("Illegal use of DELAY Shift \n");
				return(1);
			}
		}

		/* check for shift  */
		if(newtok==SHIFTSYM || newtok==ISHIFTSYM) {
			temp=my_symb[toklist[lstptr+1]].com;
			/* !! */
			if(is_uvar(temp) || is_ucon(temp)) {
				/* ram -- same issue */
				if(is_uvar(temp)) {
					my_symb[LASTTOK].com=COM(SVARTYPE, temp%MAXTYPE);
				}
				if(is_ucon(temp)) {
					my_symb[LASTTOK].com=COM(SCONTYPE, temp%MAXTYPE);
				}
				/* create a temporary sybol */
				toklist[lstptr+1]=LASTTOK;
				my_symb[LASTTOK].pri=10;
			} else {
				printf("Illegal use of SHIFT \n");
				return(1);
			}
		}
next:
		if((newtok==ENDTOK) && (oldtok==STARTTOK)) {
			break;
		}
		if(newtok==LPAREN) {
			tokstak[tokptr]=LPAREN;
			tokptr++;
			oldtok=LPAREN;
			goto getnew;
		}
		if(newtok==RPAREN) {
			switch(oldtok) {
			case LPAREN:
				tokptr--;
				oldtok=tokstak[tokptr-1];
				goto getnew;
			case COMMA:
				tokptr--;
				ncomma++;
				oldtok=tokstak[tokptr-1];
				goto next;
			}
		}
		if((newtok==COMMA) && (oldtok==COMMA)) {
			tokstak[tokptr]=COMMA;
			tokptr++;
			goto getnew;
		}
		/* ram -- the THOUS problem */
		/*           if(my_symb[oldtok%THOUS].pri>=my_symb[newtok%THOUS].pri)
		   {
			command[comptr]=my_symb[oldtok%THOUS].com;
		if((my_symb[oldtok%THOUS].arg==2) &&
		(my_symb[oldtok%THOUS].com/MAXTYPE==FUN2TYPE)) */
		if(my_symb[oldtok].pri>=my_symb[newtok].pri) {
			command[comptr]=my_symb[oldtok].com;
			if((my_symb[oldtok].arg==2) &&
			   (my_symb[oldtok].com/MAXTYPE==FUN2TYPE)) {
				ncomma--;
			}
			my_com=command[comptr];
			comptr++;
			if(my_com==NUMSYM) {
				tokptr--;
				command[comptr]=tokstak[tokptr-1];
				comptr++;
				tokptr--;
				command[comptr]=tokstak[tokptr-1];
				comptr++;
			}
			if(my_com==SUMSYM) {
				loopstk[lptr]=comptr;
				comptr++;
				lptr++;
				ncomma-=1;
			}
			if(my_com==ENDSUM) {
				lptr--;
				jmp=comptr-loopstk[lptr]-1;
				command[loopstk[lptr]]=jmp;
			}
			if(my_com==MYIF) {
				loopstk[lptr]=comptr; /* add some space for jump */
				comptr++;
				lptr++;
				nif++;
			}
			if(my_com==MYTHEN) {
				/* First resolve the if jump */
				lptr--;
				jmp=comptr-loopstk[lptr];  /* -1 is old */
				command[loopstk[lptr]]=jmp;
				/* Then set up for the then jump */
				loopstk[lptr]=comptr;
				lptr++;
				comptr++;
				nthen++;
			}
			if(my_com==MYELSE) {
				lptr--;
				jmp=comptr-loopstk[lptr]-1;
				command[loopstk[lptr]]=jmp;
				nelse++;
			}
			if(my_com==ENDDELAY || my_com==ENDSHIFT || my_com==ENDISHIFT) {
				ncomma-=1;
			}
			if(my_com==ENDDELSHFT) {
				ncomma-=2;
			}
			/*    CHECK FOR USER FUNCTION       */
			if(is_ufun(my_com))	{
				my_arg=my_symb[oldtok].arg;
				command[comptr]=my_arg;
				comptr++;
				ncomma=ncomma+1-my_arg;
			}
			/*      USER FUNCTION OKAY          */
			tokptr--;
			oldtok=tokstak[tokptr-1];
			goto next;
		}
		if(newtok==NUMTOK) {
			tokstak[tokptr++]=toklist[lstptr++];
			tokstak[tokptr++]=toklist[lstptr++];
		}
		tokstak[tokptr]=newtok;
		oldtok=newtok;
		tokptr++;
		goto getnew;
	}
	if(ncomma!=0) {
		plintf("Illegal number of arguments\n");
		return(1);
	}
	if((nif!=nelse) || (nif!=nthen)) {
		plintf("If statement missing ELSE or THEN \n");
		return(1);
	}
	command[comptr]=my_symb[ENDTOK].com;
	return(0);
}


void show_where(char *string, int index) {
	char junk[MAXEXPLEN];
	int i;
	/* exit(-1); */
	for(i=0;i<index;i++) {
		junk[i]=' ';
	}
	junk[index]='^';
	junk[index+1]=0;
	plintf("%s\n%s\n",string,junk);
}


/* functions should have ( after them  */
int function_sym(int token) {
	int com=my_symb[token].com;
	int i1=com/MAXTYPE;

	if(i1==FUN1TYPE && !unary_sym(token)) {/* single variable functions */
		return(1);
	}
	if(i1==FUN2TYPE && !binary_sym(token)) {/* two-variable function */
		return(1);
	}
	/* ram this was: if(i1==UFUN || i1==7 || i1==6 || i1==5)return(1); recall: 5 was bad */
	if (i1 == UFUNTYPE || i1 == TABTYPE || i1==VECTYPE || i1 == NETTYPE) {
		return(1);
	}
	if(token==DELSHFTSYM || token==DELSYM || token==SHIFTSYM || token==ISHIFTSYM ||
	   com==MYIF || com==MYTHEN || com==MYELSE || com==SUMSYM || com==ENDSUM) {
		return(1);
	}
	return(0);
}


int unary_sym(int token) {
	/* ram: these are tokens not byte code, so no change here? */
	if(token==9 || token==55) {
		return(1);
	}
	return(0);
}


int binary_sym(int token) {
	/* ram: these are tokens not byte code, so no change here? */
	if(token>2 && token<9) {
		return(1);
	}
	if(token>43 && token<51) {
		return(1);
	}
	if(token==54) {
		return(1);
	}
	return(0);
}


int pure_number(int token) {
	int com=my_symb[token].com;
	int i1=com/MAXTYPE;
	/* !! */
	if(token==NUMTOK || isvar(i1) || iscnst(i1) || isker(i1) || i1==USTACKTYPE || token==INDX) {
		return(1);
	}
	return(0);
}


int gives_number(int token) {
	int com=my_symb[token].com;
	int i1=com/MAXTYPE;
	if(token==INDX || token==NUMTOK) {
		return(1);
	}
	if((i1==FUN1TYPE && !unary_sym(token)) ||
	   (i1==FUN2TYPE && !binary_sym(token))) {
		return(1); /* two-variable function */
	}
	if (i1 == USTACKTYPE || isvar(i1) || iscnst(i1) || i1 == TABTYPE ||
		i1==VECTYPE || i1 == NETTYPE || isker(i1) || i1 == UFUNTYPE) {
		return(1);
	}
	if(com==MYIF || token==DELSHFTSYM || token==DELSYM || token==SHIFTSYM ||
	   token==ISHIFTSYM || com==SUMSYM) {
		return(1);
	}
	return(0);
}


/* 1 is BAD!   */
int check_syntax(int oldtoken, int newtoken) {
	int com2=my_symb[newtoken].com;
	/* if the first symbol or (  or binary symbol then must be unary symbol or
   something that returns a number or another (
*/
	if(unary_sym(oldtoken) || oldtoken==COMMA || oldtoken==STARTTOK
	   || oldtoken==LPAREN || binary_sym(oldtoken)) {
		if(unary_sym(newtoken) || gives_number(newtoken) || newtoken==LPAREN) {
			return(0);
		}
		return(1);
	}
	/* if this is a regular function, then better have ( */
	if(function_sym(oldtoken)) {
		if(newtoken==LPAREN) {
			return(0);
		}
		return(1);
	}

	/* if we have a constant or variable or ) or kernel then better
   have binary symbol or "then" or "else" as next symbol
*/
	if(pure_number(oldtoken)) {
		if(binary_sym(newtoken) || newtoken==RPAREN ||
		   newtoken==COMMA || newtoken==ENDTOK) {
			return(0);
		}
		return(1);
	}
	if(oldtoken==RPAREN) {
		if(binary_sym(newtoken) || newtoken==RPAREN ||
		   newtoken==COMMA || newtoken==ENDTOK) {
			return(0);
		}
		if(com2==MYELSE || com2==MYTHEN || com2==ENDSUM) {
			return(0);
		}
		return(1);
	}
	plintf("Bad token %d \n",oldtoken);
	return(1);
}


/******************************
*    PARSER                   *
******************************/
int make_toks(char *dest, int *my_token) {
	char num[40];
	double value;
	int old_tok=STARTTOK,tok_in=0;
	int index=0,token,nparen=0,lastindex=0;
	/*  WARNING  -- ASSUMES 32 bit int  and 64 bit double  */
	union {
		struct {
			int int1;
			int int2;
		} pieces;
		struct {
			double z;
		} num;
	} encoder;

	while(dest[index]!='\0') {
		lastindex=index;
		find_tok(dest,&index,&token);
		if((token==MINUS) &&
		   ((old_tok==STARTTOK) || (old_tok==COMMA) || (old_tok==LPAREN))) {
			token=NEGATE;
		}
		if(token==LPAREN) {
			++nparen;
		}
		if(token==RPAREN) {
			--nparen;
		}
		if(token==NSYM)	{
			if(do_num(dest,num,&value,&index)) {
				show_where(dest,index);
				return(1);
			}
			/*    new code        3/95      */
			encoder.num.z=value;
			my_token[tok_in++]=NUMTOK;
			my_token[tok_in++]=encoder.pieces.int1;
			my_token[tok_in++]=encoder.pieces.int2;
			if(check_syntax(old_tok,NUMTOK)==1) {
				plintf("Illegal syntax \n");
				show_where(dest,lastindex);
				return(1);
			}
			old_tok=NUMTOK;
		} else {
			my_token[tok_in++]=token;
			if(check_syntax(old_tok,token)==1) {
				plintf("Illegal syntax (Ref:%d %d) \n",old_tok,token);
				show_where(dest,lastindex);
				tokeninfo(old_tok);
				tokeninfo(token);
				return(1);
			}
			old_tok=token;
		}
	}
	my_token[tok_in++]=ENDTOK;
	if(check_syntax(old_tok,ENDTOK)==1) {
		plintf("Premature end of expression \n");
		show_where(dest,lastindex);
		return(1);
	}
	if(nparen!=0) {
		if(ERROUT)printf(" parentheses don't match\n");
		return(1);
	}
	return(0);
}


void tokeninfo(int tok) {
	plintf(" %s %d %d %d %d \n",
		   my_symb[tok].name,my_symb[tok].len,my_symb[tok].com,
		   my_symb[tok].arg,my_symb[tok].pri);
}


int do_num(char *source, char *num, double *value, int *ind) {
	int j=0,i=*ind,error=0;
	int ndec=0,nexp=0,ndig=0;
	char ch,oldch;
	oldch='\0';
	*value=0.0;
	while(1) {
		ch=source[i];
		if(((ch=='+') || (ch=='-')) &&
		   (oldch!='E')) {
			break;
		}
		if((ch=='*')  || (ch=='^') || (ch=='/') || (ch==',') || (ch==')') ||
		   (ch=='\0') || (ch=='|') || (ch=='>') || (ch=='<') || (ch=='&') ||
		   (ch=='=')) {
			break;
		}
		if((ch=='E') || (ch=='.') || (ch=='+') || (ch=='-') || isdigit(ch)) {
			if(isdigit(ch)) {
				ndig++;
			}
			switch(ch) {
			case 'E':
				nexp++;
				if((nexp==2) || (ndig==0)) {
					goto err;
				}
				break;
			case '.':
				ndec++;
				if((ndec==2) || (nexp==1)) {
					goto err;
				}
				break;
			}
			num[j]=ch;
			j++;
			i++;
			oldch=ch;
		} else {
err:
			num[j]=ch;
			j++;
			error=1;
			break;
		}
	}
	num[j]='\0';
	if(error==0) {
		*value=atof(num);
	} else {
		if(ERROUT) {
			printf(" illegal expression: %s\n",num);
		}
	}
	*ind=i;
	return(error);
}


void convert(char *source, char *dest) {
	char ch;
	int i=0,j=0;
	while (1) {
		ch=source[i];
		if(!isspace(ch)) {
			dest[j++]=ch;
		}
		i++;
		if(ch=='\0') {
			break;
		}
	}
	strupr(dest);
}


void find_tok(char *source, int *index, int *tok) {
	int i=*index,maxlen=0,symlen;
	int k,j,my_tok,match;
	my_tok=NSYM;
	for(k=0;k<NSYM;k++) {
		symlen=my_symb[k].len;
		if(symlen<=maxlen) {
			continue;
		}
		match=1;
		for(j=0;j<symlen;j++) {
			if(source[i+j]!=my_symb[k].name[j])	{
				match=0;
				break;
			}
		}
		if(match!=0) {
			my_tok=k;
			maxlen=symlen;
		}
	}
	*index=*index+maxlen;
	*tok=my_tok;
}


double pmod(double x, double y) {
	double z=fmod(x,y);
	if(z<0) {
		z+=y;
	}
	return(z);
}


void two_args(void) {
	fun2[4]=atan2;
	fun2[5]=pow;
	fun2[6]=max;
	fun2[7]=min;
	fun2[8]=pmod; /* This always gives an answer in [0,y) for mod(x,y) */
	fun2[9]=dand;
	fun2[10]=dor;
	fun2[11]=dgt;
	fun2[12]=dlt;
	fun2[13]=deq;
	fun2[14]=dge;
	fun2[15]=dle;
	fun2[16]=dne;
	fun2[17]=normal;
	fun2[18]=bessel_j;
	fun2[19]=bessel_y;
	fun2[20]=bessi;
}


/* These are the Bessel Functions; if you dont have them then return
 * some sort of dummy value or else write a program to compute them
*/

double bessel_j(double x, double y) {
	int n=(int)x;
	return(jn(n,y));
}

double bessel_y(double x, double y) {
	int n=(int)x;
	return(yn(n,y));
}

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

double bessi(double nn, double x) {
	int j,n;
	double  bi,bim,bip,tox,ans;
	n=(int)nn;

	if(n==0) {
		return bessi0(x);
	} else if(n==1) {
		return bessi1(x);
	}
	if (x == 0.0) {
		return 0.0;
	} else {
		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;
		for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if (fabs(bi) > BIGNO) {
				ans *= BIGNI;
				bi *= BIGNI;
				bip *= BIGNI;
			}
			if (j == n) {
				ans=bip;
			}
		}
		ans *= bessi0(x)/bi;
		return x < 0.0  && (n & 1) ? -ans : ans;
	}
}


double bessi0(double x) {
	double  ax,ans;
	double y;
	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
											 +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
											  +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
																				 +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
																													  +y*0.392377e-2))))))));
	}
	return ans;
}

double bessi1(double x) {
	double ax,ans;
	double y;
	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
												   +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	} else {
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
											 -y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
										   +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans *= (exp(ax)/sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

#undef ACC
#undef BIGNO
#undef BIGNI


/*********************************************
		  FANCY DELAY HERE                   *-------------------------<<<
*********************************************/
double do_shift(double shift, double variable) {
	int it, in;
	int i=(int)(variable),ish=(int)shift;

	if(i<0) {
		return(0.0);
	}
	it=i/MAXTYPE;
	in = (i % MAXTYPE) + ish;
	switch(it) {
	case CONTYPE:
		if(in>NCON) {
			return 0.0;
		} else {
			return constants[in];
		}
	case VARTYPE:
		if(in>MAXODE) {
			return 0.0;
		} else {
			return variables[in];
		}
	default:
		plintf("This can't happen: Invalid symbol index for SHIFT: i = %d\n", i);
		return 0.0;
	}
}

double do_ishift(double shift, double variable) {
	return variable+shift;
}


double do_delay_shift(double delay, double shift, double variable) {
	int in;
	int i=(int)(variable);
	int ish=(int)shift;

	if(i<0) {
		return(0.0);
	}
	in=(i % MAXTYPE)+ish;
	if(in>MAXODE) {
		return 0.0;
	}
	if(del_stab_flag>0) {
		if(DelayFlag && delay>0.0) {
			return(get_delay(in-1,delay));
		}
		return(variables[in]);
	}
	return(delay_stab_eval(delay,in));
}


double do_delay(double delay, double i) {
	int variable;
	/* ram - this was a little weird, since i is a double... except I think it's secretely an integer */
	variable = ((int) i) % MAXTYPE;

	if(del_stab_flag>0) {
		if(DelayFlag && delay>0.0) {
			return(get_delay(variable-1,delay));
		}
		return(variables[variable]);
	}
	return(delay_stab_eval(delay,(int)variable));
}


double hom_bcs(int i) {
	return 0.0; /* this is deprecated so no longer used */
}


void one_arg(void) {
	fun1[0]=sin;
	fun1[1]=cos;
	fun1[2]=tan;
	fun1[3]=asin;
	fun1[4]=acos;
	fun1[5]=atan;
	fun1[6]=sinh;
	fun1[7]=tanh;
	fun1[8]=cosh;
	fun1[9]=fabs;
	fun1[10]=exp;
	fun1[11]=log;
	fun1[12]=log10;
	fun1[13]=sqrt;
	fun1[14]=neg;
	fun1[15]=recip;
	fun1[16]=heaviside;
	fun1[17]=signum;
	fun1[18]=floor;
	fun1[19]=rndom;
	fun1[20]=dnot;
	fun1[21]=erf;
	fun1[22]=erfc;
	fun1[23]=hom_bcs;
	fun1[24]=poidev;
	fun1[25]=lgamma;
}


double normal(double mean, double std) {
	double fac,r,v1,v2;
	if(BoxMullerFlag==0) {
		do {
			v1=2.0*ndrand48()-1.0;
			v2=2.0*ndrand48()-1.0;
			r=v1*v1+v2*v2;
		} while(r>=1.0);
		fac=sqrt(-2.0*log(r)/r);
		BoxMuller=v1*fac;
		BoxMullerFlag=1;
		return(v2*fac*std+mean);
	} else {
		BoxMullerFlag=0;
		return(BoxMuller*std+mean);
	}
}


double max(double x, double y) {
	return(((x>y)?x:y));
}


double min(double x, double y) {
	return(((x<y)?x:y));
}


double neg(double z) {
	return(-z);
}

double recip(double z) {
	return(1.00/z);
}


double heaviside(double z) {
	double w=1.0;
	if(z<0) {
		w=0.0;
	}
	return(w);
}


double rndom(double z) {
	return(z*ndrand48());
}


double signum(double z) {
	if(z<0.0) {
		return(-1.0);
	} else if(z>0.0) {
		return(1.0);
	} else {
		return(0.0);
	}
}


/*  logical stuff  */
double dnot(double x) {
	return((double)(x==0.0));
}


double dand(double x, double y) {
	return((double)(x && y));
}


double dor(double x, double y) {
	return((double)(x || y));
}

double dge(double x, double y) {
	return((double)(x>=y));
}

double dle(double x, double y) {
	return((double)(x<=y));
}


double deq(double x, double y) {
	return((double)(x==y));
}


double dne(double x, double y) {
	return((double)(x!=y));
}


double dgt(double x, double y) {
	return((double)(x>y));
}


double dlt(double x, double y) {
	return((double)(x<y));
}
/*              end of logical stuff    */

double evaluate(int *equat) {
	uptr=0;
	stack_pointer=0;
	return(eval_rpn(equat));
}


double eval_rpn(int *equat) {
	int i,it,in,j,*tmpeq;
	int is;

	int low,high,ijmp;
	double temx,temy,temz;
	double sum;
	/*  WARNING  -- ASSUMES 32 bit int  and 64 bit double  */
	union {
		struct {
			int int1;
			int int2;
		} pieces;
		struct {
			double z;
		} num;
	} encoder;

	while((i=*equat++)!=ENDEXP) {
		switch(i) {
		case NUMSYM:
			encoder.pieces.int2=*equat++;
			encoder.pieces.int1=*equat++;
			PUSH(encoder.num.z);
			break;
		case ENDFUN:
			i=*equat++;
			uptr-=i;
			break;
		case MYIF:
			temx=POP;
			ijmp=*equat++;
			if(temx==0.0) {
				equat+=ijmp;
			}
			break;
		case MYTHEN:
			ijmp=*equat++;
			equat+=ijmp;
			break;
		case MYELSE:
			break;
		case ENDDELSHFT:
			temx=POP;
			temy=POP;
			temz=POP;
			PUSH(do_delay_shift(temx,temy,temz));
			break;
		case ENDDELAY:
			temx=POP;
			temy=POP;
			PUSH(do_delay(temx,temy));
			break;
		case ENDSHIFT:
			temx=POP;
			temy=POP;
			PUSH(do_shift(temx,temy));
			break;
		case ENDISHIFT:
			temx=POP;
			temy=POP;
			PUSH(do_ishift(temx,temy));
			break;
		case SUMSYM:
			temx=POP;
			high=(int)temx;
			temx=POP;
			low=(int)temx;
			ijmp=*equat++;
			sum=0.0;
			if(low<=high) {
				for(is=low;is<=high;is++) {
					tmpeq=equat;
					constants[SumIndex]=(double)is;
					sum+=eval_rpn(tmpeq);
				}
			}
			equat+=ijmp;
			PUSH(sum);
			break;
		case ENDSUM:
			return(POP);
		case INDXCOM:
			PUSH(CurrentIndex);
			break;
		default:
		{
			it=i/MAXTYPE;
			in=i%MAXTYPE;
			switch(it) {
			case FUN1TYPE: PUSH(fun1[in](POP));
				break;
			case FUN2TYPE:
			{
				if(in==0) {
					temx=POP;
					temy=POP;
					PUSH(temx+temy);
					goto bye;
				}
				if(in==1) {
					temx=POP;
					temy=POP;
					PUSH(temy-temx);
					goto bye;
				}
				if(in==2) {
					temx=POP;
					temy=POP;
					PUSH(temx*temy);
					goto bye;
				}
				if(in==3) {
					temx=POP;
					if(temx==0.0) {
						temx=DOUB_EPS;
					}
					temy=POP;
					PUSH(temy/temx);
					goto bye;
				}
				temx=POP;
				temy=POP;
				PUSH(fun2[in](temy,temx));
				break;
			}
			case CONTYPE:
				PUSH(constants[in]);
				break;
			case VECTYPE:
				PUSH(vector_value(POP,in));
				break;
			case NETTYPE:
				PUSH(network_value(POP,in));
				break;
			case TABTYPE:
				PUSH(lookup(POP,in));
				break;
			case USTACKTYPE:
				/* ram: so this means ustacks really do need to be of USTACKTYPE */
				PUSH(ustack[uptr-1-in]);
				break;
			case KERTYPE:
				PUSH(ker_val(in));
				break;
			case VARTYPE:
				PUSH(variables[in]);
				break;
				/* indexes for shift and delay operators... */
			case SCONTYPE:
				PUSH((double)(COM(CONTYPE, in)));
				break;
			case SVARTYPE:
				PUSH((double)(COM(VARTYPE, in)));
				break;
			case UFUNTYPE: i=*equat++;
				for(j=0;j<i;j++) {
					ustack[uptr]=POP;
					uptr++;
				}
				PUSH(eval_rpn(ufun[in]));
				break;
			}
bye: j=0;
		}
		}
	}
	return(POP);
}
/* code for log-gamma if you dont have it */


#ifdef NOLGAMMA
double lgamma(double xx) {
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
						  24.01409824083091,-1.231739572450155,
						  0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) {
		ser += cof[j]/++y;
	}
	return -tmp+log(2.5066282746310005*ser/x);
}
#endif
