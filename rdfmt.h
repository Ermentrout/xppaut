#ifndef _rdfmt_h_
#define _rdfmt_h_

#include <stdlib.h>
#include "f2c.h"
#include "fio.h"
#include "fmt.h"
#include "fp.h"


/* rdfmt.c */
int rd_ed(struct syl *p, char *ptr, ftnlen len);
int rd_ned(struct syl *p);
int rd_I(Uint *n, int w, ftnlen len, register int base);
int rd_L(ftnint *n, int w);
int rd_F(ufloat *p, int w, int d, ftnlen len);
int rd_A(char *p, ftnlen len);
int rd_AW(char *p, int w, ftnlen len);
int rd_H(int n, char *s);
int rd_POS(char *s);


#endif
