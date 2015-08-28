
#ifndef _iio_h_
#define _iio_h_

#include "f2c.h"

/* iio.c */
int z_getc(void);
int z_putc(int c);
int z_rnew(void);
integer s_rsfi(icilist *a);
integer s_wsfi(icilist *a);
int c_si(icilist *a);
int z_wnew(void);
integer e_rsfi(void);
integer e_wsfi(void);
int y_ierr(void);


#endif
