#ifndef _histogram_h_
#define _histogram_h_


/* histogram.c */
int two_d_hist(int col1, int col2, int ndat, int n1, int n2, double xlo, double xhi, double ylo, double yhi);
void four_back(void);
void hist_back(void);
void new_four(int nmodes, int col);
int new_2d_hist(void);
void new_hist(int nbins, double zlo, double zhi, int col, int col2, char *condition, int which);
void column_mean(void);
int get_col_info(int *col, char *prompt);
void compute_power(void);
int spectrum(float *data, int nr, int win, int w_type, float *pow);
int cross_spectrum(float *data, float *data2, int nr, int win, int w_type, float *pow, int type);
void compute_sd(void);
void compute_fourier(void);
void compute_correl(void);
void compute_stacor(void);
void mycor(float *x, float *y, int n, double zlo, double zhi, int nbins, float *z, int flag);
void mycor2(float *x, float *y, int n, int nbins, float *z, int flag);
void compute_hist(void);
void sft(float *data, float *ct, float *st, int nmodes, int grid);
void fftxcorr(float *data1, float *data2, int length, int nlag, float *cr, int flag);
void fft(float *data, float *ct, float *st, int nmodes, int length);
void post_process_stuff();
void just_fourier(int flag);
void just_sd(int flag);

#endif

