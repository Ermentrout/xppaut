#ifndef _open_h_
#define _open_h_

#include "f2c.h"
#include "fio.h"
/* open.c */
integer f_open(olist *a);
int fk_open(int seq, int fmt, ftnint n);
int isdev(char *s);


#endif
