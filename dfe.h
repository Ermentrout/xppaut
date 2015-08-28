#ifndef _dfe_h_
#define _dfe_h_


#include "err.h"

#include "f2c.h"
#include "fio.h"
#include "fmt.h"


/* dfe.c */
integer s_rdfe(cilist *a);
integer s_wdfe(cilist *a);
integer e_rdfe(void);
integer e_wdfe(void);
int c_dfe(cilist *a);
int y_rsk(void);
int y_getc(void);
int y_putc(int c);
int y_rev(void);
int y_err(void);
int y_newrec(void);


#endif

