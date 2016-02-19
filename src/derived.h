#ifndef XPPAUT_DERIVED_H
#define XPPAUT_DERIVED_H


/* --- Types --- */
typedef struct {
  int index,*form;
  char *rhs;
  double value;
} DERIVED;

/* --- Functions --- */
void free_derived(void);
int compile_derived(void);
void evaluate_derived(void);
int add_derived(char *name, char *rhs);

#endif /* XPPAUT_DERIVED_H */
