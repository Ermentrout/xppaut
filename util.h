#ifndef _util_h_
#define _util_h_

#include "f2c.h"
#include "fio.h"

/* util.c */
void g_char(char *a, ftnlen alen, char *b);
void b_char(char *a, char *b, ftnlen blen);
int inode(char *a);
void mvgbt(int n, int len, char *a, char *b);


#endif

