#ifndef _wrtfmt_h_
#define _wrtfmt_h_

#include <stdlib.h>
#include "f2c.h"
#include "fio.h"
#include "fmt.h"
/* wrtfmt.c */
int mv_cur(void);
int w_ed(struct syl *p, char *ptr, ftnlen len);
int w_ned(struct syl *p);
int wrt_I(Uint *n, int w, ftnlen len, register int base);
int wrt_IM(Uint *n, int w, int m, ftnlen len);
int wrt_AP(char *s);
int wrt_H(int a, char *s);
int wrt_L(Uint *n, int len, ftnlen sz);
int wrt_A(char *p, ftnlen len);
int wrt_AW(char *p, int w, ftnlen len);
int wrt_G(ufloat *p, int w, int d, int e, ftnlen len);


#endif
