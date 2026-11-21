/*
 *  lyapunov.c     Lyapunov exponents of the current system by the QR method
 *                 (Numerics -> stocHast -> Liapunov)
 *
 *  The exponents are computed as in Benettin, Galgani, Giorgilli and Strelcyn (1980) and Dieci, Russell
 *  and Van Vleck, SIAM J. Numer. Anal. 34 (1997):  the state y and k tangent vectors Q (n x k) are
 *  integrated together,
 *
 *        y' = f(t,y),      Q' = Df(t,y) Q,
 *
 *  and every tau time units Q is re-orthonormalised by a QR factorisation,  Q <- Q R.  The exponents are
 *
 *        lambda_j = (1/T) * sum over the intervals of  log R_jj .
 *
 *  What makes this different from the older code is that it does not use the integrator chosen in the
 *  Numerics menu (whose tolerance may be loose) and does not work from the stored single precision
 *  orbit.  Everything is done in double precision with its own adaptive Dormand-Prince 5(4) integrator
 *  on the combined system (tolerance 1e-9 by default), so the answer does not depend on Dt or on the
 *  method.  Df q is obtained without a Jacobian, as the central difference
 *
 *        Df(y) q  =  ( f(y + h q) - f(y - h q) ) / (2 h),
 *
 *  with h = 1e-6 (1+|y|)/|q|, which is accurate to about 1e-9.
 *
 *  It starts from the last point of the stored orbit and continues that trajectory.  For maps (Method
 *  Discrete) the same is done with the iteration  y -> F(y)  and the exponents are per unit of time
 *  (per iteration when Dt=1).
 *
 *  Reported:  the exponents, an error estimate from 10 blocks of the run, the sum of the exponents
 *  compared with the average divergence of the flow (when all n are computed: the two must agree),
 *  the Kaplan-Yorke dimension (when it can be found), and, if asked, the running estimates
 *  lambda_j(t) in the data browser to show convergence.
 *
 *  Not available for delay or Volterra equations, Markov variables, global events, Wiener
 *  (stochastic) terms or a Poincare map: the stored data are not a trajectory of an ODE in R^n.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <X11/Xlib.h>
#include "xpplim.h"
#include "struct.h"

#define MAX_LEN_SBOX 25
#define LY_VOLTERRA 6
#define LY_MAXSTEPS 200000000L
#define LY_MAXREC 1000
#define LY_BLOCKS 10

extern float **storage;
extern int storind, NODE, NEQ, METHOD, INFLAG, POIMAP, DelayFlag, NMarkov, NFlags, NWiener, Xup, MAXSTOR;
extern double DELTA_T, MyData[MAXODE];
extern int (*rhs)();

void err_msg(char *string);
int plintf(char *fmt, ...);
void ping(void);
void bottom_msg(int line, char *msg);
int script_is_active(void);
int my_abort(void);
int do_string_box(int n, int row, int col, char *title, char **names, char values[][MAX_LEN_SBOX], int maxchar);
void set_browser_data(float **data, int col0);
void refresh_browser(int length);
int do_range(double *x, int flag);

typedef struct {
  int range, k, save;
  double T, ttrans, tau, tol;
} LYSET;

static LYSET ly = {0, 0, 1, 0.0, 0.0, 0.0, 1.e-9};
static int ly_node = -1;          /* the number of equations the settings were made for */

static float **ly_run = NULL;      /* the running estimates, shown in the browser */
static int ly_nrun = 0, ly_runcols = 0;

/* range over a parameter */
static int LIAP_FLAG = 0, LIAP_K = 1, LIAP_I = 0, LIAP_N = 0;
static float **my_liap = NULL;

/* ---------------------------------------------------------------------- */
/*   the combined system and its integrator                                */
/* ---------------------------------------------------------------------- */

typedef struct {
  int n, k, N;
  double *k1, *k2, *k3, *k4, *k5, *k6, *k7, *y, *zn, *yp, *ym, *fa, *fb;
  double tol, eps;
  long nfev;
} LYW;

static void ly_free(LYW *w)
{
  free(w->k1); free(w->k2); free(w->k3); free(w->k4); free(w->k5); free(w->k6); free(w->k7);
  free(w->y); free(w->zn); free(w->yp); free(w->ym); free(w->fa); free(w->fb);
  memset(w, 0, sizeof(*w));
}

static int ly_alloc(LYW *w, int n, int k)
{
  int N = n * (k + 1);
  double **p[13];
  int i;
  memset(w, 0, sizeof(*w));
  w->n = n; w->k = k; w->N = N; w->eps = 1.e-6;
  p[0] = &w->k1; p[1] = &w->k2; p[2] = &w->k3; p[3] = &w->k4; p[4] = &w->k5; p[5] = &w->k6; p[6] = &w->k7;
  p[7] = &w->y; p[8] = &w->zn;
  p[9] = &w->yp; p[10] = &w->ym; p[11] = &w->fa; p[12] = &w->fb;
  for (i = 0; i < 9; i++) {
    *p[i] = (double *)calloc((size_t)N, sizeof(double));
    if (*p[i] == NULL) { ly_free(w); return 0; }
  }
  for (i = 9; i < 13; i++) {
    *p[i] = (double *)calloc((size_t)n, sizeof(double));
    if (*p[i] == NULL) { ly_free(w); return 0; }
  }
  return 1;
}

static double ly_norm(int n, const double *v)
{
  int i;
  double s = 0.0;
  for (i = 0; i < n; i++) s += v[i] * v[i];
  return sqrt(s);
}

/* Df(y) q by a central difference; f = f(y) is already known in F[0..n-1] */
static void ly_dirderiv(LYW *w, double t, const double *y, const double *q, double *out)
{
  int i, n = w->n;
  double nq = ly_norm(n, q), h;
  if (!(nq > 0.0)) { for (i = 0; i < n; i++) out[i] = 0.0; return; }
  h = w->eps * (1.0 + ly_norm(n, y)) / nq;
  for (i = 0; i < n; i++) { w->yp[i] = y[i] + h * q[i]; w->ym[i] = y[i] - h * q[i]; }
  (*rhs)(t, w->yp, w->fa, n);
  (*rhs)(t, w->ym, w->fb, n);
  for (i = 0; i < n; i++) out[i] = (w->fa[i] - w->fb[i]) / (2.0 * h);
  w->nfev += 2;
}

/* the right-hand side of the combined system Z = (y, q_1, ..., q_k) */
static void ly_F(LYW *w, double t, const double *Z, double *F)
{
  int j, n = w->n;
  (*rhs)(t, (double *)Z, F, n);
  w->nfev++;
  for (j = 0; j < w->k; j++)
    ly_dirderiv(w, t, Z, Z + n * (j + 1), F + n * (j + 1));
}

/* one Dormand-Prince 5(4) step of size h from (t,Z) to Zn; returns the error in units of the tolerance */
static double ly_step(LYW *w, double t, const double *Z, double h, double *Zn)
{
  int i, N = w->N;
  double *y = w->y, err = 0.0, e, sc, a, b;
  double *k1 = w->k1, *k2 = w->k2, *k3 = w->k3, *k4 = w->k4, *k5 = w->k5, *k6 = w->k6, *k7 = w->k7;

  ly_F(w, t, Z, k1);
  for (i = 0; i < N; i++) y[i] = Z[i] + h * (1. / 5 * k1[i]);
  ly_F(w, t + h / 5, y, k2);
  for (i = 0; i < N; i++) y[i] = Z[i] + h * (3. / 40 * k1[i] + 9. / 40 * k2[i]);
  ly_F(w, t + 3 * h / 10, y, k3);
  for (i = 0; i < N; i++) y[i] = Z[i] + h * (44. / 45 * k1[i] - 56. / 15 * k2[i] + 32. / 9 * k3[i]);
  ly_F(w, t + 4 * h / 5, y, k4);
  for (i = 0; i < N; i++)
    y[i] = Z[i] + h * (19372. / 6561 * k1[i] - 25360. / 2187 * k2[i] + 64448. / 6561 * k3[i] - 212. / 729 * k4[i]);
  ly_F(w, t + 8 * h / 9, y, k5);
  for (i = 0; i < N; i++)
    y[i] = Z[i] + h * (9017. / 3168 * k1[i] - 355. / 33 * k2[i] + 46732. / 5247 * k3[i] + 49. / 176 * k4[i]
                       - 5103. / 18656 * k5[i]);
  ly_F(w, t + h, y, k6);
  for (i = 0; i < N; i++)
    Zn[i] = Z[i] + h * (35. / 384 * k1[i] + 500. / 1113 * k3[i] + 125. / 192 * k4[i] - 2187. / 6784 * k5[i]
                        + 11. / 84 * k6[i]);
  ly_F(w, t + h, Zn, k7);
  for (i = 0; i < N; i++) {
    e = h * (71. / 57600 * k1[i] - 71. / 16695 * k3[i] + 71. / 1920 * k4[i] - 17253. / 339200 * k5[i]
             + 22. / 525 * k6[i] - 1. / 40 * k7[i]);
    a = fabs(Z[i]); b = fabs(Zn[i]);
    sc = w->tol + w->tol * (a > b ? a : b);
    e = fabs(e) / sc;
    if (!(e == e)) return 1.e300;
    if (e > err) err = e;
  }
  return err;
}

static double ly_newh(double h, double err)
{
  double fac = (err <= 1.e-12) ? 5.0 : 0.9 * pow(err, -0.2);
  if (fac > 5.0) fac = 5.0;
  if (fac < 0.2) fac = 0.2;
  return h * fac;
}

/* integrate Z from t to t+tau.  *hg carries the step size from one call to the next.  0 = ok,
   -1 = overflow, -2 = too many steps, -3 = the step size collapsed */
static int ly_advance(LYW *w, double t, double *Z, double tau, double *hg)
{
  double tt = 0.0, h = *hg, hh, err;
  int i, last;
  long steps = 0;
  if (!(h > 0.0) || h > tau) h = tau / 4.0;
  while (tt < tau) {
    last = 0; hh = h;
    if (tt + hh >= tau * (1.0 - 1.e-13)) { hh = tau - tt; last = 1; }
    err = ly_step(w, t + tt, Z, hh, w->zn);
    if (err > 1.e299) return -1;
    if (++steps > LY_MAXSTEPS) return -2;
    if (err <= 1.0) {
      for (i = 0; i < w->N; i++) Z[i] = w->zn[i];
      tt = last ? tau : tt + hh;
      {
        double hn = ly_newh(hh, err);
        if (!last) h = hn;
        else if (hh > 0.5 * h && hn < h) h = hn;     /* a step cut short to hit the end says nothing about h */
      }
    } else {
      h = ly_newh(hh, err);
      if (h < 1.e-13 * tau) return -3;
    }
  }
  *hg = h;
  return 0;
}

/* one iteration of a map and of its tangent map: y -> F(y), q -> DF(y) q (rhs returns the new state) */
static void ly_mapstep(LYW *w, double t, double *Z)
{
  int i, j, n = w->n, k = w->k;
  double *fnew = w->k2, *tang = w->k3;     /* scratch of length N >= n*(k+1) */
  (*rhs)(t, Z, fnew, n);
  w->nfev++;
  for (j = 0; j < k; j++) ly_dirderiv(w, t, Z, Z + n * (j + 1), tang + n * j);
  for (i = 0; i < n; i++) Z[i] = fnew[i];
  for (i = 0; i < n * k; i++) Z[n + i] = tang[i];
}

/* ---------------------------------------------------------------------- */
/*   QR                                                                    */
/* ---------------------------------------------------------------------- */

/* orthonormalise the k columns of Q (each of length n) in place by Gram-Schmidt with
   re-orthogonalisation; logr[j] = log of the length of column j after removing the earlier ones.
   0 = ok, -1 = a column vanished or the numbers overflowed */
static int ly_qr(int n, int k, double *Q, double *logr)
{
  int i, j, p, pass;
  double c, r;
  for (j = 0; j < k; j++) {
    double *v = Q + n * j;
    for (pass = 0; pass < 2; pass++)
      for (p = 0; p < j; p++) {
        const double *q = Q + n * p;
        c = 0.0;
        for (i = 0; i < n; i++) c += q[i] * v[i];
        for (i = 0; i < n; i++) v[i] -= c * q[i];
      }
    r = ly_norm(n, v);
    if (!(r > 0.0) || !(r < 1.e300)) return -1;
    for (i = 0; i < n; i++) v[i] /= r;
    logr[j] = log(r);
  }
  return 0;
}

/* a fixed, generic, orthonormal set of starting vectors (a symmetric problem can have an invariant
   coordinate axis, so the axes are not a safe choice) */
static void ly_startvecs(int n, int k, double *Q)
{
  int i, j;
  unsigned long s = 12345UL;
  for (j = 0; j < k; j++)
    for (i = 0; i < n; i++) {
      s = s * 1103515245UL + 12345UL;
      Q[n * j + i] = (double)((s >> 8) & 0xFFFF) / 65536.0 - 0.5;
    }
  {
    double *lr = (double *)calloc((size_t)(k > 0 ? k : 1), sizeof(double));
    if (lr) { ly_qr(n, k, Q, lr); free(lr); }
  }
}

/* the divergence  tr Df(y)  by n central differences */
static double ly_trace(LYW *w, double t, const double *y)
{
  int i, n = w->n;
  double h = w->eps * (1.0 + ly_norm(n, y)), tr = 0.0;
  double *x = w->zn;                       /* scratch of length >= n */
  for (i = 0; i < n; i++) x[i] = y[i];
  for (i = 0; i < n; i++) {
    x[i] = y[i] + h; (*rhs)(t, x, w->fa, n);
    x[i] = y[i] - h; (*rhs)(t, x, w->fb, n);
    x[i] = y[i];
    tr += (w->fa[i] - w->fb[i]) / (2.0 * h);
  }
  w->nfev += 2 * n;
  return tr;
}

/* ---------------------------------------------------------------------- */
/*   the computation                                                       */
/* ---------------------------------------------------------------------- */

typedef struct {
  int k, have_err, have_div, ismap, nint;
  double T, avgdiv, tau;
  double *lam, *err;
} LYRES;

static char ly_msg[300];

static int ly_valid(int batch)
{
  char *m = NULL;
  if (INFLAG == 0 || storind < 2) m = "Liapunov: integrate the equations first; the exponents are computed along the orbit that continues from its last point";
  else if (METHOD == LY_VOLTERRA) m = "Liapunov: not available for Volterra equations";
  else if (DelayFlag) m = "Liapunov: not available for delay equations";
  else if (NMarkov > 0) m = "Liapunov: not available with Markov variables";
  else if (NFlags > 0) m = "Liapunov: not available with global events (flags)";
  else if (NWiener > 0) m = "Liapunov: not available with Wiener (noise) terms";
  else if (POIMAP) m = "Liapunov: the Poincare map is on, so the stored data are not a trajectory";
  else if (DELTA_T < 0.0) m = "Liapunov: integrate forward in time first";
  else if (NODE < 1 || NODE > NEQ) m = "Liapunov: no equations";
  if (m) {
    if (!batch) err_msg(m);
    return 0;
  }
  return 1;
}

static void ly_defaults(void)
{
  double dt = fabs(DELTA_T);
  if (ly_node != NODE) {
    ly.k = (NODE <= 8) ? NODE : 1;
    ly_node = NODE;
    ly.range = 0;
  }
  if (ly.k > NODE) ly.k = NODE;
  if (ly.k < 1) ly.k = 1;
  if (!(ly.tau > 0.0)) ly.tau = (dt > 0.0) ? dt : 1.0;
  if (!(ly.tol > 0.0)) ly.tol = 1.e-9;
}

/* the exponents along the trajectory that continues from the last stored point.
   0 = failed (message in ly_msg), 1 = ok */
static int ly_compute(const LYSET *s, LYRES *R, int keeprun)
{
  int n = NODE, k = s->k, i, j, ismap = (METHOD == 0);
  double t0 = (double)storage[0][storind - 1];
  double dt = fabs(DELTA_T), Tm, tau, hg = 0.0;
  double *Z = NULL, *logr = NULL, *sum = NULL, *bsum = NULL, *bt = NULL;
  double tt, trprev = 0.0, trsum = 0.0;
  LYW W;
  long m, nint, rec_every, ntr = 0;
  int ntrans, nb, bi = 0, bn = 0, rc = 0, nrec = 0, dotrace = 0, ok = 0;

  memset(&W, 0, sizeof(W));
  Tm = (s->T > 0.0) ? s->T : (double)(storage[0][storind - 1] - storage[0][0]);
  if (!(Tm > 0.0)) { sprintf(ly_msg, "Liapunov: the stored orbit has no length; give a time"); return 0; }
  if (!(dt > 0.0)) dt = 1.0;
  tau = (s->tau > 0.0) ? s->tau : dt;
  if (ismap) {                       /* a whole number of iterations per renormalisation */
    double nr = floor(tau / dt + 0.5);
    if (nr < 1.0) nr = 1.0;
    tau = nr * dt;
    nint = (long)floor(Tm / tau + 0.5);
    if (nint < 2) { sprintf(ly_msg, "Liapunov: too few iterations; give a longer time"); return 0; }
    Tm = nint * tau;
  } else {
    if (tau > Tm) tau = Tm;
    nint = (long)ceil(Tm / tau - 1.e-9);
    if (nint < 1) nint = 1;
    tau = Tm / (double)nint;
  }
  if (nint > 100000000L) {
    sprintf(ly_msg, "Liapunov: the renormalisation interval is too small for this time (%ld intervals)", nint);
    return 0;
  }
  ntrans = (s->ttrans > 0.0) ? (int)ceil(s->ttrans / tau) : 0;
  nb = (nint >= 2 * LY_BLOCKS) ? LY_BLOCKS : 0;
  rec_every = nint / LY_MAXREC; if (rec_every < 1) rec_every = 1;
  dotrace = (!ismap && k == n && n <= 60);

  if (!ly_alloc(&W, n, k)) { sprintf(ly_msg, "Liapunov: out of memory"); return 0; }
  W.tol = s->tol;
  Z = (double *)calloc((size_t)W.N, sizeof(double));
  logr = (double *)calloc((size_t)k, sizeof(double));
  sum = (double *)calloc((size_t)k, sizeof(double));
  bsum = (double *)calloc((size_t)((nb > 0 ? nb : 1) * k), sizeof(double));
  bt = (double *)calloc((size_t)(nb > 0 ? nb : 1), sizeof(double));
  R->lam = (double *)calloc((size_t)k, sizeof(double));
  R->err = (double *)calloc((size_t)k, sizeof(double));
  if (!Z || !logr || !sum || !bsum || !bt || !R->lam || !R->err) {
    sprintf(ly_msg, "Liapunov: out of memory");
    goto done;
  }
  for (i = 0; i < n; i++) Z[i] = (double)storage[i + 1][storind - 1];
  ly_startvecs(n, k, Z + n);
  R->k = k; R->T = Tm; R->tau = tau; R->nint = (int)nint; R->ismap = ismap;

  if (keeprun) {
    int ncols = NEQ + 1;
    nrec = (int)(nint / rec_every) + 1;
    if (ly_run) { for (i = 0; i < ly_runcols; i++) free(ly_run[i]); free(ly_run); ly_run = NULL; }
    ly_run = (float **)calloc((size_t)ncols, sizeof(float *));
    ly_runcols = ncols;
    if (ly_run) for (i = 0; i < ncols; i++) {
      ly_run[i] = (float *)calloc((size_t)nrec, sizeof(float));
      if (!ly_run[i]) { for (j = 0; j < i; j++) free(ly_run[j]); free(ly_run); ly_run = NULL; break; }
    }
    ly_nrun = 0;
  }

  /* ----- transient: let the tangent vectors line up ----- */
  tt = t0;
  for (m = 0; m < ntrans; m++) {
    if (ismap) {
      for (j = 0; j < (int)(tau / dt + 0.5); j++) {
        ly_mapstep(&W, tt, Z);
        tt += dt;
        for (i = 0; i < n; i++) if (!(fabs(Z[i]) < 1.e150)) { rc = -1; goto fail; }
      }
    } else {
      rc = ly_advance(&W, tt, Z, tau, &hg);
      if (rc < 0) goto fail;
      tt += tau;
    }
    if (ly_qr(n, k, Z + n, logr) < 0) { rc = -1; goto fail; }
  }

  /* ----- the measurement ----- */
  if (dotrace) trprev = ly_trace(&W, tt, Z);
  for (m = 0; m < nint; m++) {
    if (ismap) {
      for (j = 0; j < (int)(tau / dt + 0.5); j++) {
        ly_mapstep(&W, tt, Z);
        tt += dt;
        for (i = 0; i < n; i++) if (!(fabs(Z[i]) < 1.e150)) { rc = -1; goto fail; }
      }
    } else {
      rc = ly_advance(&W, tt, Z, tau, &hg);
      if (rc < 0) goto fail;
      tt += tau;
    }
    if (ly_qr(n, k, Z + n, logr) < 0) { rc = -1; goto fail; }
    for (j = 0; j < k; j++) {
      sum[j] += logr[j];
      if (nb > 0) bsum[bi * k + j] += logr[j];
    }
    if (nb > 0) {
      bn++;
      if (bn >= nint / nb && bi < nb - 1) { bt[bi] = (double)bn * tau; bi++; bn = 0; }
    }
    if (dotrace) {
      double tr = ly_trace(&W, tt, Z);
      trsum += 0.5 * (trprev + tr) * tau;
      trprev = tr; ntr++;
    }
    if (ly_run && ((m + 1) % rec_every == 0) && ly_nrun < nrec) {
      double el = (double)(m + 1) * tau;
      ly_run[0][ly_nrun] = (float)el;
      for (j = 0; j < k; j++) ly_run[j + 1][ly_nrun] = (float)(sum[j] / el);
      ly_nrun++;
    }
    if (Xup && !script_is_active() && ((m & 63) == 63)) {
      int ch = my_abort();
      if (ch == 27 || ch == '/') { sprintf(ly_msg, "Liapunov: aborted"); goto done; }
    }
  }
  if (nb > 0) bt[bi] = (double)(bn > 0 ? bn : 1) * tau;

  for (j = 0; j < k; j++) R->lam[j] = sum[j] / Tm;

  /* the error estimate from the blocks */
  if (nb > 0) {
    int b;
    R->have_err = 1;
    for (j = 0; j < k; j++) {
      double mean = 0.0, var = 0.0, r;
      for (b = 0; b < nb; b++) { r = bsum[b * k + j] / bt[b]; mean += r; }
      mean /= nb;
      for (b = 0; b < nb; b++) { r = bsum[b * k + j] / bt[b] - mean; var += r * r; }
      var /= (nb - 1);
      R->err[j] = sqrt(var / nb);
    }
  }
  if (dotrace && ntr > 0) { R->have_div = 1; R->avgdiv = trsum / Tm; }
  ok = 1;
  goto done;

fail:
  if (rc == -1) sprintf(ly_msg, "Liapunov: the solution blew up or a tangent vector collapsed (is the orbit on an attractor?)");
  else if (rc == -2) sprintf(ly_msg, "Liapunov: too many steps (a stiff system? try a larger tolerance or a shorter time)");
  else sprintf(ly_msg, "Liapunov: the step size became too small (a stiff system or a singularity)");
done:
  ly_free(&W);
  free(Z); free(logr); free(sum); free(bsum); free(bt);
  if (!ok) { free(R->lam); free(R->err); R->lam = R->err = NULL; }
  return ok;
}

/* the Kaplan-Yorke dimension, or -1 if the exponents computed are not enough */
static double ly_ky(int k, const double *lam, int *known)
{
  int j;
  double s = 0.0;
  *known = 0;
  if (lam[0] < 0.0) { *known = 1; return 0.0; }
  for (j = 0; j < k; j++) {
    if (s + lam[j] < 0.0) {
      *known = 1;
      return (double)j + s / fabs(lam[j]);
    }
    s += lam[j];
  }
  return -1.0;
}

/* ---------------------------------------------------------------------- */
/*   the menu                                                              */
/* ---------------------------------------------------------------------- */

static void ly_report(const LYRES *R)
{
  char line[300], tmp[80];
  int j, known;
  double dky, sum = 0.0;
  plintf("\nLyapunov exponents%s: %d of %d, along %g time units (%d intervals of %g), tolerance %g\n",
         R->ismap ? " (map, per unit time)" : "", R->k, NODE, R->T, R->nint, R->tau, ly.tol);
  for (j = 0; j < R->k; j++) {
    sum += R->lam[j];
    if (R->have_err) plintf("  lambda %d = %12.6g   +- %.2g\n", j + 1, R->lam[j], R->err[j]);
    else plintf("  lambda %d = %12.6g\n", j + 1, R->lam[j]);
  }
  if (R->k > 1) plintf("  sum      = %12.6g\n", sum);
  if (R->have_div) plintf("  average divergence of the flow = %12.6g  (it should equal the sum)\n", R->avgdiv);
  dky = ly_ky(R->k, R->lam, &known);
  if (known) plintf("  Kaplan-Yorke dimension = %.5g\n", dky);
  else if (R->k > 1) plintf("  Kaplan-Yorke dimension: not found (more exponents are needed)\n");
  sprintf(line, "Lyapunov:");
  for (j = 0; j < R->k && j < 4; j++) {
    sprintf(tmp, " %.5g", R->lam[j]);
    strcat(line, tmp);
  }
  if (R->k > 4) strcat(line, " ...");
  if (R->have_err) { sprintf(tmp, "  (+-%.2g)", R->err[0]); strcat(line, tmp); }
  if (R->have_div) { sprintf(tmp, "  sum %.5g, div %.5g", sum, R->avgdiv); strcat(line, tmp); }
  if (known && R->k > 1) { sprintf(tmp, "  D_KY %.4g", dky); strcat(line, tmp); }
  strcpy(ly_msg, line);
}

static int ly_dialog(void)
{
  static char *names[] = {"Range(0/1)", "Nexp", "Time(0=orbit)", "Extra transient", "Renorm interval", "Tolerance", "Save run(0/1)"};
  char values[7][MAX_LEN_SBOX];
  int status, k, range, save;
  double T, tt, tau, tol;
  ly_defaults();
  sprintf(values[0], "%d", ly.range);
  sprintf(values[1], "%d", ly.k);
  sprintf(values[2], "%g", ly.T);
  sprintf(values[3], "%g", ly.ttrans);
  sprintf(values[4], "%g", ly.tau);
  sprintf(values[5], "%g", ly.tol);
  sprintf(values[6], "%d", ly.save);
  status = do_string_box(7, 7, 1, "Lyapunov exponents", names, values, 40);
  if (status == 0) return 0;
  range = atoi(values[0]); k = atoi(values[1]);
  T = atof(values[2]); tt = atof(values[3]); tau = atof(values[4]); tol = atof(values[5]); save = atoi(values[6]);
  if (k < 1) k = 1;
  if (k > NODE) k = NODE;
  if (tt < 0.0) tt = 0.0;
  if (!(tol > 0.0)) tol = 1.e-9;
  if (tol < 1.e-13) tol = 1.e-13;
  if (tol > 1.e-3) tol = 1.e-3;
  ly.range = (range != 0); ly.k = k; ly.T = (T > 0.0) ? T : 0.0; ly.ttrans = tt; ly.tau = tau; ly.tol = tol; ly.save = (save != 0);
  return 1;
}

void do_liapunov(void)
{
  LYRES R;
  int i, j;
  if (!ly_valid(0)) return;
  if (!ly_dialog()) return;
  if (ly.range) {
    double *x = &MyData[0];
    LIAP_FLAG = 1;
    LIAP_K = ly.k;
    do_range(x, 0);
    /* the range is done: put (parameter, lambda_1 ... lambda_k) in the browser's data */
    if (my_liap && LIAP_I > 0) {
      plintf("\nLyapunov exponents over the range (columns 1..%d of the browser hold lambda_1..lambda_%d):\n", LIAP_K, LIAP_K);
      for (i = 0; i < LIAP_I && i < 500; i++) {
        plintf("  %12.6g", (double)my_liap[0][i]);
        for (j = 1; j <= LIAP_K && j <= 6; j++) plintf("  %12.6g", (double)my_liap[j][i]);
        plintf("\n");
      }
      for (i = 0; i < LIAP_I && i < MAXSTOR && i < LIAP_N; i++)
        for (j = 0; j <= LIAP_K && j <= NEQ; j++)
          storage[j][i] = my_liap[j][i];
      storind = LIAP_I;
      refresh_browser(storind);
    }
    if (my_liap) {
      for (j = 0; j <= LIAP_K; j++) free(my_liap[j]);
      free(my_liap);
      my_liap = NULL;
    }
    LIAP_FLAG = 0;
    return;
  }
  memset(&R, 0, sizeof(R));
  ping();
  if (!ly_compute(&ly, &R, ly.save)) { err_msg(ly_msg); return; }
  ly_report(&R);
  if (ly.save && ly_run && ly_nrun > 1) {
    set_browser_data(ly_run, 1);
    refresh_browser(ly_nrun);
  }
  /* In a script an error box is a failure that stops the script, so a result must not be reported that
     way: the full report is already on the terminal, and the summary goes to the message line. */
  if (script_is_active()) bottom_msg(0, ly_msg);
  else err_msg(ly_msg);
  free(R.lam); free(R.err);
}

/* ---------- range over parameters: called from the range integration ---------- */

void alloc_liap(int n)
{
  int j;
  if (LIAP_FLAG == 0) return;
  if (my_liap) { for (j = 0; j <= LIAP_K; j++) free(my_liap[j]); free(my_liap); my_liap = NULL; }
  my_liap = (float **)calloc((size_t)(LIAP_K + 1), sizeof(float *));
  if (!my_liap) { LIAP_FLAG = 0; return; }
  for (j = 0; j <= LIAP_K; j++) {
    my_liap[j] = (float *)calloc((size_t)(n + 1), sizeof(float));
    if (!my_liap[j]) { LIAP_FLAG = 0; return; }
  }
  LIAP_N = n + 1;
  LIAP_I = 0;
}

void do_this_liaprun(int i, double p)
{
  LYRES R;
  LYSET s = ly;
  int j;
  if (LIAP_FLAG == 0 || my_liap == NULL) return;
  if (i >= LIAP_N) return;
  s.k = LIAP_K;
  memset(&R, 0, sizeof(R));
  my_liap[0][i] = (float)p;
  if (!ly_valid(1)) { for (j = 1; j <= LIAP_K; j++) my_liap[j][i] = 0.0f; LIAP_I++; return; }
  if (ly_compute(&s, &R, 0)) {
    for (j = 0; j < LIAP_K; j++) my_liap[j + 1][i] = (float)R.lam[j];
    free(R.lam); free(R.err);
  } else {
    for (j = 1; j <= LIAP_K; j++) my_liap[j][i] = 0.0f;
    plintf("p=%g: %s\n", p, ly_msg);
  }
  LIAP_I++;
}

/* the old entry point: the largest exponent for the stored orbit.  Returns 1 on success. */
int hrw_liapunov(double *liap, int batch, double eps)
{
  LYRES R;
  LYSET s;
  (void)eps;
  if (!ly_valid(batch)) return 0;
  ly_defaults();
  s = ly; s.k = 1; s.T = 0.0; s.ttrans = 0.0;
  memset(&R, 0, sizeof(R));
  if (!ly_compute(&s, &R, 0)) { if (!batch) err_msg(ly_msg); return 0; }
  *liap = R.lam[0];
  free(R.lam); free(R.err);
  return 1;
}
