#ifndef _fmt_h_
#define _fmt_h_

#include "f2c.h"

struct syl
{	int op,p1,p2,p3;
};
#define RET 1
#define REVERT 2
#define GOTO 3
#define X 4
#define SLASH 5
#define STACK 6
#define I 7
#define ED 8
#define NED 9
#define IM 10
#define APOS 11
#define H 12
#define TL 13
#define TR 14
#define T 15
#define COLON 16
#define S 17
#define SP 18
#define SS 19
#define P 20
#define BN 21
#define BZ 22
#define F 23
#define E 24
#define EE 25
#define D 26
#define G 27
#define GE 28
#define L 29
#define A 30
#define AW 31
#define O 32
#define NONL 33
extern struct syl syl[];
extern int pc,parenlvl,revloc;
extern int (*doed)(),(*doned)();
extern int (*dorevert)(),(*donewrec)(),(*doend)();
extern flag cblank,cplus,workdone, nonl;
extern int dummy();
extern char *fmtbuf;
extern int scale;
typedef union
{	real pf;
	doublereal pd;
} ufloat;
typedef union
{	short is;
	char ic;
	int il;
} Uint;
#define GET(x) if((x=(*getn)())<0) return(x)
#define VAL(x) (x!='\n'?x:' ')
#define PUT(x) (*putn)(x)

extern int cursor;


int pars_f(char *s);
char *f_s(char *s, int curloc);
char *f_list(char *s);
char *i_tem(char *s);
int ne_d(char *s, char **p);
int e_d(char *s, char **p);
int op_gen(int a, int b, int c, int d);
char *gt_num(char *s, int *n);
integer do_fio(ftnint *number, char *ptr, ftnlen len);
int en_fio(void);
void fmt_bg(void);
int type_f(int n);
char *ap_end(char *s);



#endif
