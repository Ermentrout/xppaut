#ifndef XPPAUT_SIMPLENET_H
#define XPPAUT_SIMPLENET_H

/* --- Functions --- */
int add_spec_fun(char *name, char *rhs);
void add_special_name(char *name, char *rhs);
int add_vectorizer(char *name, char *rhs);
void add_vectorizer_name(char *name, char *rhs);
void eval_all_nets(void);
void evaluate_network(int ind);
double network_value(double x, int i);
void update_all_ffts(void);
double vector_value(double x, int i);

#endif /* XPPAUT_SIMPLENET_H */
