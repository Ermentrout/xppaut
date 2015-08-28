#ifndef _wsfe_h_
#define _wsfe_h_

#include "err.h"
#include "fio.h"
#include "fmt.h"
#include "f2c.h"

/* wsfe.c */
integer s_wsfe(cilist *a);
int x_putc(int c);
void pr_put(int c);
int x_wSL(void);
int xw_end(void);
int xw_rev(void);


#endif

