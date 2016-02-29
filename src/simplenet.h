#ifndef XPPAUT_SIMPLENET_H
#define XPPAUT_SIMPLENET_H

/* --- Macros --- */
#define MAXW 50

/* --- Types --- */
typedef struct {
	int type,ncon,n;
	char name[20];
	char soname[256],sofun[256];
	int root,root2;
	int f[20];
	int iwgt;
	int *gcom; /* for group commands */
	double *values,*weight,*index,*taud; /* for delays  */
	double *fftr,*ffti,*dr,*di;
	double *wgtlist[MAXW];
} NETWORK;

typedef struct {
	char name[20];
	int root,length,il,ir;
} VECTORIZER;

/* --- Functions --- */
int get_vector_info();
double net_interp(double x, int i);
double network_value(double x, int i);
void init_net(double *v, int n);
int add_spec_fun(char *name, char *rhs);
void add_special_name(char *name, char *rhs);
int add_vectorizer(char *name, char *rhs);
void add_vectorizer_name(char *name, char *rhs);
int is_network(char *s);
void eval_all_nets(void);
void evaluate_network(int ind);
void update_all_ffts(void);
void update_fft(int ind);
void fft_conv(int it, int n, double *values, double *yy, double *fftr, double *ffti, double *dr, double *di);
int gilparse(char *s, int *ind, int *nn);
int g_namelist(char *s, char *root, int *flag, int *i1, int *i2);

#endif /* XPPAUT_SIMPLENET_H */
