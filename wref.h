#ifndef _wref_h_
#define _wref_h_

#include "fmt.h"

#include "f2c.h"
#include "fio.h"
#include "fp.h"
#include <string.h>

/* wref.c */
int wrt_E(ufloat *p, int w, int d, int e, ftnlen len);
int wrt_F(ufloat *p, int w, int d, ftnlen len);


#endif
