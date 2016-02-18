#ifndef XPPAUT_STORAGE_H
#define XPPAUT_STORAGE_H

/* --- Data --- */
extern double *WORK;

extern int IWORK[10000];
extern int MAXSTOR;

/* --- Functions --- */
void init_alloc_info(void);
void alloc_meth(void);
void init_stor(int nrow, int ncol);
void free_storage(int ncol);
int reallocstor(int ncol,int nrow);

#endif /* XPPAUT_STORAGE_H */
