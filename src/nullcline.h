#ifndef XPPAUT_NULLCLINE_H
#define XPPAUT_NULLCLINE_H

#include <stdio.h>


/* --- Data --- */
extern int DF_FLAG;
extern int DF_GRID;
extern int DFBatch;
extern int DOING_DFIELD;
extern int NCBatch;
extern int XNullColor;
extern int YNullColor;

/* --- Functions --- */
void create_new_cline(void);
void direct_field_com(int c);
void do_batch_dfield(void);
void do_batch_nclines(void);
void do_range_clines(void);
void froz_cline_stuff_com(int i);
int get_nullcline_floats(float **v, int *n, int who, int type);
void new_clines_com(int c);
void redraw_dfield(void);
void restore_nullclines(void);
void silent_dfields(void);
void silent_nullclines(void);

#endif /* XPPAUT_NULLCLINE_H */
