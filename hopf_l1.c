/*
 *  hopf_l1.c      First Lyapunov coefficient at a Hopf point of  x' = f(x)
 *
 *  Yu. A. Kuznetsov, Elements of Applied Bifurcation Theory, 3rd edn, Springer 2004, Sec. 5.4 (the projection formula):
 *
 *     l1 = Re[ <p,C(q,q,qbar)> - 2<p,B(q,A^-1 B(q,qbar))> + <p,B(qbar,(2i w I - A)^-1 B(q,q))> ] / (2 w)
 *
 *  with  A q = i w q,  A^T p = -i w p,  <p,q> = 1.  A is the Jacobian at the equilibrium x and B, C are the
 *  second and third derivatives of f, here from central differences of the vector field along the directions
 *  in question. l1 < 0: the Hopf bifurcation is supercritical (a stable cycle is born); l1 > 0: subcritical.
 *
 *  Nothing here knows about AUTO's data structures. hopf_l1() takes the equilibrium and uses the right-hand side
 *  rhs() with the parameters as they are now. hopf_l1_par() is for use from the AUTO diagram code: it first puts
 *  AUTO's parameter vector into XPP's parameters (as func() in autpp.c does) and puts them back afterwards.
 *
 *  Only for ordinary differential equations (not maps) with n = NODE equations.
 */

#include <stdio.h>
extern int DelayFlag;
#include <stdlib.h>
#include <string.h>
#include <math.h>

int plintf(char *fmt, ...);
extern int (*rhs)();
extern int NODE, NAutoPar, Auto_index_to_array[8];
extern double constants[];
void eigen(int n, double *a, double *ev, double *work, int *ierr);
void evaluate_derived(void);
void redo_all_fun_tables(void);

/* solve M x = b for an m x m matrix M (row-major, overwritten); b is overwritten by x. 1 = ok */
static int hl_solve(int m, double *M, double *b)
{
  int i, j, k, p;
  double big, t;
  for (k = 0; k < m; k++) {
    p = k; big = fabs(M[k * m + k]);
    for (i = k + 1; i < m; i++) if (fabs(M[i * m + k]) > big) { big = fabs(M[i * m + k]); p = i; }
    if (!(big > 1.e-300)) return 0;
    if (p != k) {
      for (j = 0; j < m; j++) { t = M[k * m + j]; M[k * m + j] = M[p * m + j]; M[p * m + j] = t; }
      t = b[k]; b[k] = b[p]; b[p] = t;
    }
    for (i = k + 1; i < m; i++) {
      t = M[i * m + k] / M[k * m + k];
      if (t == 0.0) continue;
      for (j = k; j < m; j++) M[i * m + j] -= t * M[k * m + j];
      b[i] -= t * b[k];
    }
  }
  for (i = m - 1; i >= 0; i--) {
    t = b[i];
    for (j = i + 1; j < m; j++) t -= M[i * m + j] * b[j];
    b[i] = t / M[i * m + i];
  }
  return 1;
}

/* solve (Br + i Bi)(wr + i wi) = (br + i bi) for n x n real Br, Bi (row-major). 1 = ok */
static int hl_csolve(int n, const double *Br, const double *Bi, const double *br, const double *bi,
                     double *wr, double *wi)
{
  int i, j, m = 2 * n, ok;
  double *M = (double *)malloc(sizeof(double) * (size_t)m * (size_t)m), *b = (double *)malloc(sizeof(double) * (size_t)m);
  if (!M || !b) { free(M); free(b); return 0; }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      M[i * m + j] = Br[i * n + j];           M[i * m + n + j] = -Bi[i * n + j];
      M[(n + i) * m + j] = Bi[i * n + j];     M[(n + i) * m + n + j] = Br[i * n + j];
    }
  for (i = 0; i < n; i++) { b[i] = br[i]; b[n + i] = bi[i]; }
  ok = hl_solve(m, M, b);
  if (ok) for (i = 0; i < n; i++) { wr[i] = b[i]; wi[i] = b[n + i]; }
  free(M); free(b);
  return ok;
}

/* an eigenvector of B for the eigenvalue near i*sig, by inverse iteration. 1 = ok */
static int hl_nullv(int n, const double *B, double sig, double *vr, double *vi)
{
  int i, j, it;
  double d = 1.e-9 * (1.0 + fabs(sig)), s;
  double *Br = (double *)malloc(sizeof(double) * (size_t)n * (size_t)n), *Bi = (double *)malloc(sizeof(double) * (size_t)n * (size_t)n);
  double *wr = (double *)malloc(sizeof(double) * (size_t)n), *wi = (double *)malloc(sizeof(double) * (size_t)n);
  int ok = 1;
  if (!Br || !Bi || !wr || !wi) { free(Br); free(Bi); free(wr); free(wi); return 0; }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) { Br[i * n + j] = B[i * n + j] - (i == j ? d : 0.0); Bi[i * n + j] = (i == j ? -sig : 0.0); }
  for (i = 0; i < n; i++) { vr[i] = 1.0 / sqrt((double)n); vi[i] = 0.0; }
  for (it = 0; it < 4 && ok; it++) {
    if (!hl_csolve(n, Br, Bi, vr, vi, wr, wi)) { ok = 0; break; }
    s = 0.0;
    for (i = 0; i < n; i++) s += wr[i] * wr[i] + wi[i] * wi[i];
    s = sqrt(s);
    if (!(s > 0.0)) { ok = 0; break; }
    for (i = 0; i < n; i++) { vr[i] = wr[i] / s; vi[i] = wi[i] / s; }
  }
  free(Br); free(Bi); free(wr); free(wi);
  return ok;
}

typedef struct { int n; double h; const double *x; double *xs, *f; } HLW;

static void hl_f(HLW *w, const double *dir1, double c1, const double *dir2, double c2, const double *dir3, double c3,
                 double *out)
{
  int k;
  for (k = 0; k < w->n; k++)
    w->xs[k] = w->x[k] + w->h * (c1 * dir1[k] + (dir2 ? c2 * dir2[k] : 0.0) + (dir3 ? c3 * dir3[k] : 0.0));
  (*rhs)(0.0, w->xs, out, w->n);
}

/* B(u,v) and C(u,v,w) for real vectors, by central differences of f */
static void hl_Bre(HLW *w, const double *u, const double *v, double *out)
{
  static const int c[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
  int q, k;
  for (k = 0; k < w->n; k++) out[k] = 0.0;
  for (q = 0; q < 4; q++) {
    hl_f(w, u, c[q][0], v, c[q][1], NULL, 0.0, w->f);
    for (k = 0; k < w->n; k++) out[k] += c[q][2] * w->f[k];
  }
  for (k = 0; k < w->n; k++) out[k] /= 4.0 * w->h * w->h;
}

static void hl_Cre(HLW *w, const double *u, const double *v, const double *z, double *out)
{
  int a, b, c, k;
  for (k = 0; k < w->n; k++) out[k] = 0.0;
  for (a = -1; a <= 1; a += 2) for (b = -1; b <= 1; b += 2) for (c = -1; c <= 1; c += 2) {
    hl_f(w, u, a, v, b, z, c, w->f);
    for (k = 0; k < w->n; k++) out[k] += a * b * c * w->f[k];
  }
  for (k = 0; k < w->n; k++) out[k] /= 8.0 * w->h * w->h * w->h;
}

static int hl_zero(int n, const double *v) { int k; for (k = 0; k < n; k++) if (v[k] != 0.0) return 0; return 1; }

/* complex bilinear form B(U,V), U = Ur + i Ui, ... */
static void hl_Bc(HLW *w, const double *Ur, const double *Ui, const double *Vr, const double *Vi, double *Or, double *Oi)
{
  int k, n = w->n;
  double *t = (double *)malloc(sizeof(double) * (size_t)n);
  for (k = 0; k < n; k++) { Or[k] = 0.0; Oi[k] = 0.0; }
  hl_Bre(w, Ur, Vr, t); for (k = 0; k < n; k++) Or[k] += t[k];
  hl_Bre(w, Ui, Vi, t); for (k = 0; k < n; k++) Or[k] -= t[k];
  hl_Bre(w, Ur, Vi, t); for (k = 0; k < n; k++) Oi[k] += t[k];
  hl_Bre(w, Ui, Vr, t); for (k = 0; k < n; k++) Oi[k] += t[k];
  free(t);
}

/* complex trilinear form C(U,V,W) */
static void hl_Cc(HLW *w, const double *Ur, const double *Ui, const double *Vr, const double *Vi,
                  const double *Wr, const double *Wi, double *Or, double *Oi)
{
  static const double ip[4][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
  int m, k, n = w->n;
  double *t = (double *)malloc(sizeof(double) * (size_t)n);
  for (k = 0; k < n; k++) { Or[k] = 0.0; Oi[k] = 0.0; }
  for (m = 0; m < 8; m++) {
    int a = m & 1, b = (m >> 1) & 1, c = (m >> 2) & 1;
    const double *u = a ? Ui : Ur, *v = b ? Vi : Vr, *z = c ? Wi : Wr;
    if (hl_zero(n, u) || hl_zero(n, v) || hl_zero(n, z)) continue;
    hl_Cre(w, u, v, z, t);
    for (k = 0; k < n; k++) { Or[k] += ip[a + b + c][0] * t[k]; Oi[k] += ip[a + b + c][1] * t[k]; }
  }
  free(t);
}

/* x: the equilibrium (n values). The parameters are those now in force. Returns 1 and sets *l1 (and *omega, the
   Hopf frequency) if x has a pair of eigenvalues +-i*omega near the imaginary axis; 0 if not or if anything fails. */
int hopf_l1x(int n, const double *x, double *l1, double *omega, double *scale)
{
  int i, j, ier = 0, ok = 0, k;
  /* storage: 5 matrices of n*n (A, AT, column-major copy, Br, Bi) and about 30 vectors of n */
  double *A = NULL, *AT = NULL, *cm = NULL, *ev = NULL, *work = NULL, *buf = NULL;
  double *qr, *qi, *p0r, *p0i, *pr, *pi, *qbi, *g11r, *g11i, *g20r, *g20i, *a11, *h20r, *h20i;
  double *c3r, *c3i, *t2r, *t2i, *t3r, *t3i, *zero, *Br, *Bi;
  double om = 0.0, best = 1.e300, big = 0.0, sr, si, s2, re, h;
  HLW w;
  if (n < 2) return 0;

  /* all the work arrays in one block */
  buf = (double *)calloc((size_t)(5 * n * n + 40 * n + 8), sizeof(double));
  if (!buf) return 0;
  A = buf; AT = A + n * n; cm = AT + n * n; Br = cm + n * n; Bi = Br + n * n;
  ev = Bi + n * n; work = ev + 2 * n;
  qr = work + n; qi = qr + n; p0r = qi + n; p0i = p0r + n; pr = p0i + n; pi = pr + n; qbi = pi + n;
  g11r = qbi + n; g11i = g11r + n; g20r = g11i + n; g20i = g20r + n; a11 = g20i + n; h20r = a11 + n; h20i = h20r + n;
  c3r = h20i + n; c3i = c3r + n; t2r = c3i + n; t2i = t2r + n; t3r = t2i + n; t3i = t3r + n; zero = t3i + n;

  for (k = 0; k < n; k++) if (fabs(x[k]) > big) big = fabs(x[k]);
  h = 1.e-3 * (1.0 + big);
  {
    double *xs = (double *)malloc(sizeof(double) * (size_t)n), *f = (double *)malloc(sizeof(double) * (size_t)n);
    double *f1 = (double *)malloc(sizeof(double) * (size_t)n), *f2 = (double *)malloc(sizeof(double) * (size_t)n);
    if (!xs || !f || !f1 || !f2) { free(xs); free(f); free(f1); free(f2); free(buf); return 0; }
    w.n = n; w.h = h; w.x = x; w.xs = xs; w.f = f;

    /* Jacobian A by central differences */
    for (j = 0; j < n; j++) {
      double d = 1.e-6 * (1.0 + fabs(x[j]));
      for (k = 0; k < n; k++) xs[k] = x[k];
      xs[j] = x[j] + d; (*rhs)(0.0, xs, f1, n);
      xs[j] = x[j] - d; (*rhs)(0.0, xs, f2, n);
      for (i = 0; i < n; i++) A[i * n + j] = (f1[i] - f2[i]) / (2.0 * d);
    }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) AT[i * n + j] = A[j * n + i];

    /* eigenvalues of A (XPP's own routine wants a column-major copy) */
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) cm[i + j * n] = A[i * n + j];
    eigen(n, cm, ev, work, &ier);
    if (ier != 0) goto done;
    for (i = 0; i < n; i++) {
      double re_ = ev[2 * i], im_ = ev[2 * i + 1];
      if (im_ > 1.e-8 * (1.0 + big) && fabs(re_) < best) { best = fabs(re_); om = im_; }
    }
    if (!(om > 0.0)) goto done;
    *omega = om;
    if (scale) { *scale = 0.0; for (i = 0; i < n; i++) { double m_ = hypot(ev[2 * i], ev[2 * i + 1]); if (m_ > *scale) *scale = m_; } }

    /* eigenvectors q (A q = i om q) and p (A^T p = -i om p), with <p,q> = conj(p).q = 1 */
    if (!hl_nullv(n, A, om, qr, qi) || !hl_nullv(n, AT, -om, p0r, p0i)) goto done;
    sr = 0.0; si = 0.0;
    for (k = 0; k < n; k++) { sr += p0r[k] * qr[k] + p0i[k] * qi[k]; si += p0r[k] * qi[k] - p0i[k] * qr[k]; }
    s2 = sr * sr + si * si;
    if (!(s2 > 0.0)) goto done;
    for (k = 0; k < n; k++) { pr[k] = (p0r[k] * sr - p0i[k] * si) / s2; pi[k] = (p0r[k] * si + p0i[k] * sr) / s2; }
    for (k = 0; k < n; k++) { qbi[k] = -qi[k]; zero[k] = 0.0; }

    /* g11 = B(q,qbar), a11 = A^-1 g11 */
    hl_Bc(&w, qr, qi, qr, qbi, g11r, g11i);
    for (k = 0; k < n; k++) a11[k] = g11r[k];
    { double *M = (double *)malloc(sizeof(double) * (size_t)n * (size_t)n);
      if (!M) goto done;
      memcpy(M, A, sizeof(double) * (size_t)n * (size_t)n);
      if (!hl_solve(n, M, a11)) { free(M); goto done; }
      free(M); }
    /* g20 = B(q,q), h20 = (2 i om I - A)^-1 g20 */
    hl_Bc(&w, qr, qi, qr, qi, g20r, g20i);
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) { Br[i * n + j] = -A[i * n + j]; Bi[i * n + j] = (i == j ? 2.0 * om : 0.0); }
    if (!hl_csolve(n, Br, Bi, g20r, g20i, h20r, h20i)) goto done;
    /* the three terms */
    hl_Cc(&w, qr, qi, qr, qi, qr, qbi, c3r, c3i);
    hl_Bc(&w, qr, qi, a11, zero, t2r, t2i);
    hl_Bc(&w, qr, qbi, h20r, h20i, t3r, t3i);
    re = 0.0;
    for (k = 0; k < n; k++) {
      double vr = c3r[k] - 2.0 * t2r[k] + t3r[k], vi = c3i[k] - 2.0 * t2i[k] + t3i[k];
      re += pr[k] * vr + pi[k] * vi;
    }
    *l1 = re / (2.0 * om);
    ok = (*l1 == *l1);
  done:
    free(xs); free(f); free(f1); free(f2);
  }
  free(buf);
  return ok;
}

int hopf_l1(int n, const double *x, double *l1, double *omega)
{
  return hopf_l1x(n, x, l1, omega, NULL);
}

/* For the AUTO diagram code: par is AUTO's parameter vector and u the state at the point (first NODE entries).
   The parameters of the model are set from par (as func() in autpp.c does) and restored afterwards, so that the
   parameters you are using are not changed by redrawing a diagram. */
static int hopf_l1x_par(const double *par, const double *u, double *l1, double *omega, double *scale)
{
  int i, ok;
  double saved[8];
  for (i = 0; i < NAutoPar && i < 8; i++) { saved[i] = constants[Auto_index_to_array[i]]; constants[Auto_index_to_array[i]] = par[i]; }
  evaluate_derived();
  redo_all_fun_tables();
  ok = hopf_l1x(NODE, u, l1, omega, scale);
  for (i = 0; i < NAutoPar && i < 8; i++) constants[Auto_index_to_array[i]] = saved[i];
  evaluate_derived();
  redo_all_fun_tables();
  return ok;
}

int hopf_l1_par(const double *par, const double *u, double *l1, double *omega)
{
  double sc;
  return hopf_l1x_par(par, u, l1, omega, &sc);
}

/* ------------------------------------------------------------------------------------------------------------
   For the drawing of two-parameter Hopf curves in the AUTO window.

   The points of a curve arrive one at a time (and again, in order, whenever the diagram is redrawn). For each one
   the first Lyapunov coefficient is found (and remembered, so that redrawing is quick) and then:
     * hopf_curve_point() gives its sign, for the colour of the curve, and tells when l1 has changed sign since the
       previous point of the same curve: a generalised Hopf (Bautin, GH) point, where the position is found by
       interpolating l1 linearly between the two points (frac is the fraction of the way from the previous point).
     * a point counts only if there is a genuine pair of eigenvalues +-i*omega and omega is not tiny compared with
       the other eigenvalues. That keeps neutral saddles, and the neighbourhood of a Bogdanov-Takens or zero-Hopf
       point (where omega goes to 0 and l1 can change sign by going through infinity), from being called GH.
   ------------------------------------------------------------------------------------------------------------ */

typedef struct hlent { double p1, p2, l1, om, sc; int icp1, icp2, ok; struct hlent *next; } HLENT;
#define HL_TAB 2048
static HLENT *hl_tab[HL_TAB];

static unsigned hl_hash(double a, double b)
{
  unsigned long long x, y;
  memcpy(&x, &a, sizeof x); memcpy(&y, &b, sizeof y);
  x ^= y * 0x9E3779B97F4A7C15ULL; x ^= x >> 29; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 32;
  return (unsigned)(x % HL_TAB);
}

/* forget the remembered values (for a new model, or if parameters that AUTO does not hold have been changed) */
void hopf_l1_cache_clear(void)
{
  int i;
  for (i = 0; i < HL_TAB; i++) {
    HLENT *e = hl_tab[i], *nx;
    while (e) { nx = e->next; free(e); e = nx; }
    hl_tab[i] = NULL;
  }
}

static HLENT *hl_find(const double *par, const double *u, int icp1, int icp2)
{
  double p1 = par[icp1], p2 = par[icp2];
  unsigned h = hl_hash(p1, p2);
  HLENT *e;
  for (e = hl_tab[h]; e; e = e->next)
    if (e->p1 == p1 && e->p2 == p2 && e->icp1 == icp1 && e->icp2 == icp2) return e;
  e = (HLENT *)calloc(1, sizeof(HLENT));
  if (!e) return NULL;
  e->p1 = p1; e->p2 = p2; e->icp1 = icp1; e->icp2 = icp2;
  e->ok = hopf_l1x_par(par, u, &e->l1, &e->om, &e->sc);
  e->next = hl_tab[h]; hl_tab[h] = e;
  return e;
}

static int hl_good(const HLENT *e) { return e && e->ok && e->om > 1.e-2 * e->sc; }

static int hl_have = 0;
static double hl_lprev = 0.0, hl_p1prev = 0.0, hl_p2prev = 0.0;
static double hl_printed[40][2];
static int hl_nprinted = 0;

/* first = 1 for the first point of a curve. Returns 1 if l1 changed sign between the previous point and this one. */
extern int DelayFlag;
int hopf_curve_point(int first, const double *par, const double *u, int icp1, int icp2,
                     const char *name1, const char *name2, int *sign, double *frac)
{
  HLENT *e;
  int gh = 0, i;
  *sign = 0; *frac = 0.0;
  /* the first Lyapunov coefficient is computed from the ODE right-hand side
     (Taylor coefficients, no delayed arguments): it is not valid for DDEs */
  if (DelayFlag) { hl_have = 0; return 0; }
  e = hl_find(par, u, icp1, icp2);
  if (first) hl_have = 0;
  if (hl_good(e)) {
    *sign = e->l1 > 0.0 ? 1 : -1;
    if (hl_have && hl_lprev * e->l1 < 0.0) {
      double f = hl_lprev / (hl_lprev - e->l1), g1 = hl_p1prev + f * (e->p1 - hl_p1prev), g2 = hl_p2prev + f * (e->p2 - hl_p2prev);
      gh = 1; *frac = f;
      for (i = 0; i < hl_nprinted; i++) if (hl_printed[i][0] == g1 && hl_printed[i][1] == g2) break;
      if (i == hl_nprinted && hl_nprinted < 40) {
        hl_printed[hl_nprinted][0] = g1; hl_printed[hl_nprinted][1] = g2; hl_nprinted++;
        plintf("GH (generalised Hopf, Bautin) near %s = %.8g, %s = %.8g: the first Lyapunov coefficient changes sign\n",
               name1, g1, name2, g2);
      }
    }
    hl_have = 1; hl_lprev = e->l1; hl_p1prev = e->p1; hl_p2prev = e->p2;
  } else hl_have = 0;
  return gh;
}

/* just the sign of l1 at a point of a Hopf curve (for the PostScript drawing): -1, +1, or 0 if there is none */
int hopf_curve_sign(const double *par, const double *u, int icp1, int icp2)
{
  HLENT *e;
  if (DelayFlag) return 0;
  e = hl_find(par, u, icp1, icp2);
  if (!hl_good(e)) return 0;
  return e->l1 > 0.0 ? 1 : -1;
}
