#ifndef XPPAUT_EXTRA_H
#define XPPAUT_EXTRA_H

#include "xpplim.h"

/* --- Data --- */
extern int dll_flag;
extern char dll_fun[DEFAULT_STRING_LENGTH];
extern char dll_lib[DEFAULT_STRING_LENGTH];

/* --- Functions --- */
void add_export_list(char *in, char *out);
void auto_load_dll(void);
void do_export_list(void);
void do_in_out(void);
void load_new_dll(void);

#endif /* XPPAUT_EXTRA_H */
