#ifndef _err_h_
#define _err_h_

#include <stdlib.h>
#include <stdio.h>

#include "f2c.h"
#include "fio.h"

/* err.c */
void fatal(int n, char *s);
void f_init(void);
int canseek(FILE *f);
int nowreading(unit *x);
int nowwriting(unit *x);


#endif
