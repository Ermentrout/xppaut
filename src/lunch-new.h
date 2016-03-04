#ifndef XPPAUT_LUNCH_NEW_H
#define XPPAUT_LUNCH_NEW_H

#include <stdio.h>
#include "form_ode.h"

/* --- Macros --- */
/* flag values */
#define WRITEM 0
#define READEM 1

/* --- Functions --- */
void do_lunch(int f);
void file_inf(void);
void io_ic_file(char *fn, int flag);
void io_parameter_file(char *fn, int flag);
int read_lunch(FILE *fp);

#endif /* XPPAUT_LUNCH_NEW_H */
