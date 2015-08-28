#ifndef _lread_h_
#define _lread_h_


#include <stdlib.h>
#include "f2c.h"
#include "fio.h"
#include "fmt.h"
#include "lio.h"
#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif


#include "axes2.h"
#include "browse.h"

#include "fp.h"


/* lread.c */
int t_getc(void);
integer e_rsle(void);
int l_read(ftnint *number, char *ptr, ftnlen len, ftnint type);
int l_R(void);
int l_C(void);
int l_L(void);
int l_CHAR(void);
integer s_rsle(cilist *a);
int c_le(cilist *a);
integer do_lio(ftnint *type, ftnint *number, char *ptr, ftnlen len);


#endif
