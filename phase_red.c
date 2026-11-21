/*
 *  phase_red.c      Phase reduction of a stable limit cycle:  the adjoint Z (the infinitesimal phase response curve)
 *                   (stage 1; the isostable quantities, the correction B and the interaction functions H, H2 follow)
 *
 *  For a cycle U(t) of  x' = f(x)  with period T, the adjoint Z(t) solves  Z' = -A(t)^T Z,  A(t) = Df(U(t)),  with
 *  Z(t).U'(t) = 1. Z(0) is the left eigenvector of the monodromy matrix for the multiplier 1.
 *
 *  The method is that of XPP Studio (and is independent of AUTO):
 *    1. the orbit on a fine grid (RK4, Ng steps of size h = T/Ng) with the Jacobian A at every grid point,
 *       by fourth-order central differences;
 *    2. the variational flow  Phi' = A Phi  by RK4 with step 2h (the midpoint is the grid point between),
 *       keeping the one-step matrices Psi_k, so that  Phi(T) = Psi_(P-2) ... Psi_0;
 *    3. Z(0) from the left eigenvector of Phi(T) for the multiplier nearest 1 (inverse iteration), normalised
 *       by Z(0).f(U(0)) = 1;
 *    4. Z carried back with the exact transposes,  Z_k = Psi_k^T Z_(k+1).  Z^T dX is then conserved to rounding
 *       and Z is periodic by construction.
 *  Before step 2 a few Newton steps, on the grid's own map with the period free, make the grid orbit exactly
 *  periodic. (The orbit from (I)nitialconds (P)eriodic is periodic for the true flow, not for RK4 with Ng steps.)
 *
 *  Only for ordinary differential equations, n = number of equations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "phase_red.h"

extern int (*rhs)();
void eigen(int n, double *a, double *ev, double *work, int *ierr);

static void ph_f(int n, const double *y, double *f) { (*rhs)(0.0, (double *)y, f, n); }

/* solve M x = b, M n x n row-major (destroyed), b overwritten. 1 = ok */
static int ph_solve(int n, double *M, double *b)
{
  int i, j, k, p;
  double big, t;
  for (k = 0; k < n; k++) {
    p = k; big = fabs(M[k * n + k]);
    for (i = k + 1; i < n; i++) if (fabs(M[i * n + k]) > big) { big = fabs(M[i * n + k]); p = i; }
    if (!(big > 1.e-300)) return 0;
    if (p != k) {
      for (j = 0; j < n; j++) { t = M[k * n + j]; M[k * n + j] = M[p * n + j]; M[p * n + j] = t; }
      t = b[k]; b[k] = b[p]; b[p] = t;
    }
    for (i = k + 1; i < n; i++) {
      t = M[i * n + k] / M[k * n + k];
      if (t == 0.0) continue;
      for (j = k; j < n; j++) M[i * n + j] -= t * M[k * n + j];
      b[i] -= t * b[k];
    }
  }
  for (i = n - 1; i >= 0; i--) {
    t = b[i];
    for (j = i + 1; j < n; j++) t -= M[i * n + j] * b[j];
    b[i] = t / M[i * n + i];
  }
  return 1;
}

typedef struct { int n, Ng, P; double h, H; double *Y, *F, *A, *end; } PHGRID;

static void ph_gfree(PHGRID *G) { free(G->Y); free(G->F); free(G->A); free(G->end); memset(G, 0, sizeof *G); }

/* one RK4 step of size h from y (in place) */
static void ph_rk4(int n, double *y, double h, double *k1, double *k2, double *k3, double *k4, double *yt)
{
  int i;
  ph_f(n, y, k1); for (i = 0; i < n; i++) yt[i] = y[i] + 0.5 * h * k1[i];
  ph_f(n, yt, k2); for (i = 0; i < n; i++) yt[i] = y[i] + 0.5 * h * k2[i];
  ph_f(n, yt, k3); for (i = 0; i < n; i++) yt[i] = y[i] + h * k3[i];
  ph_f(n, yt, k4); for (i = 0; i < n; i++) y[i] += h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}

/* the RK4 end state after one period on a grid of 200*ceil(N/200) steps, from x (no Jacobians) */
static void ph_endstate(int n, const double *x, double T, int N, double *y)
{
  int K = 200 * ((N + 199) / 200), q, i;
  double h = T / K, *w = (double *)malloc(sizeof(double) * 5 * (size_t)n);
  for (i = 0; i < n; i++) y[i] = x[i];
  for (q = 0; q < K; q++) ph_rk4(n, y, h, w, w + n, w + 2 * n, w + 3 * n, w + 4 * n);
  free(w);
}

/* the orbit grid: Ng = 200*ceil(N/200) steps; Y, F = f(Y) and A = Df(Y) at the Ng+1 points; Y[Ng] is set to x0 */
static int ph_grid(int n, const double *x0, double T, int N, PHGRID *G)
{
  int Ng = 200 * ((N + 199) / 200), nn = n * n, q, g, i, k;
  double h = T / Ng, *y, *w, *yt, *fa, *fb, e, *yg;
  memset(G, 0, sizeof *G);
  G->n = n; G->Ng = Ng; G->P = Ng / 2 + 1; G->h = h; G->H = 2 * h;
  G->Y = (double *)malloc(sizeof(double) * (size_t)(Ng + 1) * (size_t)n);
  G->F = (double *)malloc(sizeof(double) * (size_t)(Ng + 1) * (size_t)n);
  G->A = (double *)malloc(sizeof(double) * (size_t)(Ng + 1) * (size_t)nn);
  G->end = (double *)malloc(sizeof(double) * (size_t)n);
  y = (double *)malloc(sizeof(double) * (size_t)n); w = (double *)malloc(sizeof(double) * 7 * (size_t)n);
  if (!G->Y || !G->F || !G->A || !G->end || !y || !w) { free(y); free(w); ph_gfree(G); return 0; }
  yt = w + 5 * n; fa = w + 5 * n + n; /* w has 7n entries: 5 for rk4, then yt-like scratch below */
  fb = (double *)malloc(sizeof(double) * (size_t)n); yg = (double *)malloc(sizeof(double) * (size_t)n);
  for (i = 0; i < n; i++) y[i] = x0[i];
  for (q = 0; q < Ng; q++) {
    memcpy(G->Y + (size_t)q * n, y, sizeof(double) * (size_t)n);
    ph_rk4(n, y, h, w, w + n, w + 2 * n, w + 3 * n, w + 4 * n);
  }
  memcpy(G->end, y, sizeof(double) * (size_t)n);
  memcpy(G->Y + (size_t)Ng * n, x0, sizeof(double) * (size_t)n);
  for (i = 0; i < (Ng + 1) * n; i++) if (!isfinite(G->Y[i])) { free(y); free(w); free(fb); free(yg); ph_gfree(G); return 0; }
  for (g = 0; g <= Ng; g++) {
    memcpy(yg, G->Y + (size_t)g * n, sizeof(double) * (size_t)n);
    ph_f(n, yg, G->F + (size_t)g * n);
    /* fourth-order central differences: the error of A is amplified around an unstable orbit, so it must be small */
    for (k = 0; k < n; k++) {
      e = 2.e-4 * (1.0 + fabs(yg[k]));
      memcpy(yt, yg, sizeof(double) * (size_t)n);
      { double *Ag = G->A + (size_t)g * nn;
        yt[k] = yg[k] + e; ph_f(n, yt, fa); yt[k] = yg[k] - e; ph_f(n, yt, fb);
        for (i = 0; i < n; i++) Ag[i * n + k] = 8 * (fa[i] - fb[i]);
        yt[k] = yg[k] + 2 * e; ph_f(n, yt, fa); yt[k] = yg[k] - 2 * e; ph_f(n, yt, fb);
        for (i = 0; i < n; i++) Ag[i * n + k] = (Ag[i * n + k] - (fa[i] - fb[i])) / (12 * e);
        yt[k] = yg[k]; }
    }
  }
  free(y); free(w); free(fb); free(yg);
  return 1;
}

/* the step matrices Psi_k of the variational flow (RK4, step H, midpoint = the grid point between), and the product Phi */
static int ph_stepmats(int n, const PHGRID *G, double *PS, double *Phi)
{
  int nn = n * n, P = G->P, k, q, i, c, l, ok = 1;
  double H = G->H, *B = (double *)malloc(sizeof(double) * 8 * (size_t)nn);
  double *S2 = B, *S3 = B + nn, *S4 = B + 2 * nn, *W2 = B + 3 * nn, *W3 = B + 4 * nn, *W4 = B + 5 * nn, *PhiN = B + 6 * nn;
  if (!B) return 0;
  memset(Phi, 0, sizeof(double) * (size_t)nn);
  for (i = 0; i < n; i++) Phi[i * n + i] = 1.0;
#define PH_MUL(g, X, OUT) { const double *Ag_ = G->A + (size_t)(g) * nn; \
    for (i = 0; i < n; i++) for (c = 0; c < n; c++) { double sum_ = 0.0; for (l = 0; l < n; l++) sum_ += Ag_[i * n + l] * (X)[l * n + c]; (OUT)[i * n + c] = sum_; } }
  for (k = 0; k < P - 1; k++) {
    int g = 2 * k;
    double *Ps = PS + (size_t)k * nn;
    const double *A0 = G->A + (size_t)g * nn;
    for (q = 0; q < nn; q++) S2[q] = 0.5 * H * A0[q];
    for (i = 0; i < n; i++) S2[i * n + i] += 1.0;
    PH_MUL(g + 1, S2, W2);
    for (q = 0; q < nn; q++) S3[q] = 0.5 * H * W2[q];
    for (i = 0; i < n; i++) S3[i * n + i] += 1.0;
    PH_MUL(g + 1, S3, W3);
    for (q = 0; q < nn; q++) S4[q] = H * W3[q];
    for (i = 0; i < n; i++) S4[i * n + i] += 1.0;
    PH_MUL(g + 2, S4, W4);
    for (q = 0; q < nn; q++) Ps[q] = H / 6.0 * (A0[q] + 2 * W2[q] + 2 * W3[q] + W4[q]);
    for (i = 0; i < n; i++) Ps[i * n + i] += 1.0;
    for (i = 0; i < n; i++) for (c = 0; c < n; c++) { double sum = 0.0; for (l = 0; l < n; l++) sum += Ps[i * n + l] * Phi[l * n + c]; PhiN[i * n + c] = sum; }
    memcpy(Phi, PhiN, sizeof(double) * (size_t)nn);
  }
#undef PH_MUL
  for (q = 0; q < nn; q++) if (!isfinite(Phi[q])) ok = 0;
  free(B);
  return ok;
}

void phase_free(PHASE *ph)
{
  free(ph->th); free(ph->U); free(ph->F); free(ph->Z); free(ph->I); free(ph->R); free(ph->B);
  memset(ph, 0, sizeof *ph);
}

const char *phase_iso_note(int code)
{
  switch (code) {
  case 0: return "";
  case 1: return "a system of one variable has no nontrivial multiplier";
  case 2: return "the multipliers could not be found";
  case 3: return "the dominant multiplier is one of a complex pair";
  case 4: return "the dominant multiplier is numerically zero (an extremely strongly attracting orbit has no slow amplitude)";
  case 5: return "the dominant multiplier is negative (the isostable has period 2T)";
  case 6: return "the dominant multiplier is too close to 1 (near a bifurcation)";
  case 7: return "the dominant multiplier is not unique";
  case 8: return "the eigenvectors of the dominant multiplier could not be found";
  case 9: return "the isostable functions are not finite (the orbit is too unstable)";
  case 10: return "the isostable response is orthogonal to rho (degenerate)";
  }
  return "";
}

const char *phase_b_note(int code)
{
  switch (code) {
  case 0: return "";
  case 1: return "B is computed for systems of up to 24 variables";
  case 2: return "B could not be found: another multiplier equals 1/nu1";
  }
  return "";
}


/* ------------------------------------------------------------------------------------------------------------------
   Stage 2: the isostable quantities (as in XPP Studio). For the slowest nontrivial multiplier nu = exp(kappa T):
     rho  the periodic Floquet eigenfunction,  rho' = (A - kappa) rho,   so that  x = U(theta) + psi rho(theta),  psi' = kappa psi;
     I    the isostable response,              I'   = kappa I - A^T I,   normalised I.rho = 1;
     B    the first correction to the phase response off the orbit,  grad theta = Z + psi B + ...,
          B' = -(A^T + kappa) B - (D2F[rho])^T Z, the periodic solution.
   rho is carried forward, and I and B backward, with the exact step matrices of the variational flow; the components along
   f (for rho) and Z (for I) are removed at every step, since they are exactly zero.
   ------------------------------------------------------------------------------------------------------------------ */

/* an eigenvector of the n x n row-major matrix Mm (or its transpose) for the eigenvalue near nu, by inverse iteration */
static int ph_inv_iter(int n, const double *Mm, int transpose, double nu, double *w)
{
  static const double shifts[3] = {1.e-9, 1.e-7, 1.e-5};
  int s, i, j, it, nn = n * n;
  double *B = (double *)malloc(sizeof(double) * (size_t)nn), *x = (double *)malloc(sizeof(double) * (size_t)n);
  int found = 0;
  for (s = 0; s < 3 && !found; s++) {
    int good = 1;
    for (i = 0; i < n; i++) w[i] = 1.0 / (i + 1.3);
    for (it = 0; it < 6 && good; it++) {
      double nm = 0.0;
      for (i = 0; i < n; i++) for (j = 0; j < n; j++) B[i * n + j] = (transpose ? Mm[j * n + i] : Mm[i * n + j]) - (i == j ? nu * (1.0 + shifts[s]) : 0.0);
      memcpy(x, w, sizeof(double) * (size_t)n);
      if (!ph_solve(n, B, x)) { good = 0; break; }
      for (i = 0; i < n; i++) nm = fmax(nm, fabs(x[i]));
      if (!(nm > 0.0) || !isfinite(nm)) { good = 0; break; }
      for (i = 0; i < n; i++) w[i] = x[i] / nm;
    }
    if (good) found = 1;
  }
  free(B); free(x);
  return found;
}

/* the Jacobian at y, fourth-order central differences (the same as in the grid) */
static void ph_jac(int n, const double *y, double *out, double *yt, double *fa, double *fb)
{
  int i, k;
  double e;
  for (k = 0; k < n; k++) {
    e = 2.e-4 * (1.0 + fabs(y[k]));
    memcpy(yt, y, sizeof(double) * (size_t)n);
    yt[k] = y[k] + e; ph_f(n, yt, fa); yt[k] = y[k] - e; ph_f(n, yt, fb);
    for (i = 0; i < n; i++) out[i * n + k] = 8 * (fa[i] - fb[i]);
    yt[k] = y[k] + 2 * e; ph_f(n, yt, fa); yt[k] = y[k] - 2 * e; ph_f(n, yt, fb);
    for (i = 0; i < n; i++) out[i * n + k] = (out[i * n + k] - (fa[i] - fb[i])) / (12 * e);
  }
}

/* d/ds A(y + s r) at s = 0 by fourth-order central differences of the Jacobian */
static void ph_djac(int n, const double *y, const double *r, double *out, double *wk)
{
  int nn = n * n, i, q, s;
  double ym = 0.0, eps, *y1 = wk, *yt = wk + n, *fa = wk + 2 * n, *fb = wk + 3 * n, *A1 = wk + 4 * n, *A2 = A1 + nn, *A3 = A2 + nn, *A4 = A3 + nn;
  double *Ao[4]; static const double sg[4] = {1, -1, 2, -2};
  Ao[0] = A1; Ao[1] = A2; Ao[2] = A3; Ao[3] = A4;
  for (i = 0; i < n; i++) ym = fmax(ym, fabs(y[i]));
  eps = 3.e-3 * (1.0 + ym);
  for (s = 0; s < 4; s++) {
    for (i = 0; i < n; i++) y1[i] = y[i] + sg[s] * eps * r[i];
    ph_jac(n, y1, Ao[s], yt, fa, fb);
  }
  for (q = 0; q < nn; q++) out[q] = (8 * (A1[q] - A2[q]) - (A3[q] - A4[q])) / (12 * eps);
}

/* h = (dA)^T z, for dA = d/ds A(y + s r) */
static void ph_hvec(int n, const double *y, const double *z, const double *r, double *h, double *dA, double *wk)
{
  int i, c;
  ph_djac(n, y, r, dA, wk);
  for (c = 0; c < n; c++) { double sm = 0.0; for (i = 0; i < n; i++) sm += z[i] * dA[i * n + c]; h[c] = sm; }
}

static int ph_isostable(int n, const PHGRID *G, const double *PS, const double *Phi, const double *Zs, double T, PHASE *ph)
{
  int P = G->P, nn = n * n, i, j, k, ier = 0, it = 0, i1 = -1, i2 = -1;
  double H = G->H, bd = 1.e300, rate1 = 0.0, rate2 = 0.0, nu, kappa, ek, e1re, e1im, mod1, num, den;
  const double *Fv = G->F, *Y = G->Y, *A = G->A;
  double *cm = NULL, *ev = NULL, *wk = NULL, *vr = NULL, *wl = NULL, *R = NULL, *Iv = NULL, *r = NULL, *Bs = NULL;
  int ok = 0;
  ph->iso = 0; ph->isoSkip = 0; ph->hasB = 0; ph->bNote = 0;
  if (n < 2) { ph->isoSkip = 1; return 0; }
  cm = (double *)malloc(sizeof(double) * (size_t)nn); ev = (double *)malloc(sizeof(double) * 2 * (size_t)n); wk = (double *)malloc(sizeof(double) * (size_t)n);
  vr = (double *)malloc(sizeof(double) * (size_t)n); wl = (double *)malloc(sizeof(double) * (size_t)n); r = (double *)malloc(sizeof(double) * (size_t)n);
  R = (double *)calloc((size_t)P * (size_t)n, sizeof(double)); Iv = (double *)calloc((size_t)P * (size_t)n, sizeof(double));
  if (!cm || !ev || !wk || !vr || !wl || !r || !R || !Iv) goto out;
  /* the multipliers: the eigenvalues of Phi. The one nearest 1 is the trivial one; of the rest, the slowest */
  for (i = 0; i < n; i++) for (j = 0; j < n; j++) cm[i + j * n] = Phi[i * n + j];
  eigen(n, cm, ev, wk, &ier);
  if (ier != 0) { ph->isoSkip = 2; goto out; }
  for (i = 0; i < n; i++) { double d = hypot(ev[2 * i] - 1.0, ev[2 * i + 1]); if (d < bd) { bd = d; it = i; } }
  for (i = 0; i < n; i++) {
    double rt;
    if (i == it) continue;
    rt = fabs(log(hypot(ev[2 * i], ev[2 * i + 1])));
    if (i1 < 0 || rt < rate1) { i2 = i1; rate2 = rate1; i1 = i; rate1 = rt; }
    else if (i2 < 0 || rt < rate2) { i2 = i; rate2 = rt; }
  }
  e1re = ev[2 * i1]; e1im = ev[2 * i1 + 1]; mod1 = hypot(e1re, e1im);
  if (fabs(e1im) > 1.e-7 * mod1) { ph->isoSkip = 3; goto out; }
  if (fabs(e1re) < 1.e-10) { ph->isoSkip = 4; goto out; }
  if (e1re < 0.0) { ph->isoSkip = 5; goto out; }
  if (rate1 < 1.e-6 || fabs(e1re - 1.0) < 1.e-6) { ph->isoSkip = 6; goto out; }
  if (i2 >= 0 && fabs(rate2 - rate1) < 1.e-6 * rate1) { ph->isoSkip = 7; goto out; }
  if (!ph_inv_iter(n, Phi, 0, e1re, vr) || !ph_inv_iter(n, Phi, 1, e1re, wl)) { ph->isoSkip = 8; goto out; }
  /* refine nu by the Rayleigh quotient w^T Phi v / w^T v */
  num = 0.0; den = 0.0;
  for (i = 0; i < n; i++) { double sm = 0.0; for (j = 0; j < n; j++) sm += Phi[i * n + j] * vr[j]; num += wl[i] * sm; den += wl[i] * vr[i]; }
  nu = (den != 0.0 && isfinite(num / den) && num / den > 0.0) ? num / den : e1re;
  kappa = log(nu) / T;
  ek = exp(-kappa * H);
  /* rho forward */
  { double d = 0.0;
    for (i = 0; i < n; i++) d += Zs[i] * vr[i];
    for (i = 0; i < n; i++) R[i] = vr[i] - d * Fv[i]; }
  for (k = 0; k < P - 1; k++) {
    const double *Ps = PS + (size_t)k * nn; double d = 0.0;
    for (i = 0; i < n; i++) { double sm = 0.0; for (j = 0; j < n; j++) sm += Ps[i * n + j] * R[(size_t)k * n + j]; r[i] = ek * sm; d += Zs[(size_t)(k + 1) * n + i] * r[i]; }
    for (i = 0; i < n; i++) R[(size_t)(k + 1) * n + i] = r[i] - d * Fv[(size_t)2 * (k + 1) * n + i];
  }
  /* I backward */
  { double d = 0.0;
    for (i = 0; i < n; i++) d += wl[i] * Fv[(size_t)2 * (P - 1) * n + i];
    for (i = 0; i < n; i++) Iv[(size_t)(P - 1) * n + i] = wl[i] - d * Zs[(size_t)(P - 1) * n + i]; }
  for (k = P - 2; k >= 0; k--) {
    const double *Ps = PS + (size_t)k * nn; double d = 0.0;
    for (j = 0; j < n; j++) { double sm = 0.0; for (i = 0; i < n; i++) sm += Ps[i * n + j] * Iv[(size_t)(k + 1) * n + i]; r[j] = ek * sm; d += r[j] * Fv[(size_t)2 * k * n + j]; }
    for (j = 0; j < n; j++) Iv[(size_t)k * n + j] = r[j] - d * Zs[(size_t)k * n + j];
  }
  for (i = 0; i < P * n; i++) if (!isfinite(R[i]) || !isfinite(Iv[i])) { ph->isoSkip = 9; goto out; }
  /* scale: max |rho| = 1; sign: outward on average; then I.rho = 1 */
  { double s = 0.0, bv = 0.0, d0 = 0.0, sg; int kb = 0, ib = 0;
    for (k = 0; k < P; k++) { double q = 0.0; for (i = 0; i < n; i++) q += R[(size_t)k * n + i] * R[(size_t)k * n + i]; q = sqrt(q); if (q > s) { s = q; kb = k; } }
    for (i = 0; i < n; i++) if (fabs(R[(size_t)kb * n + i]) > bv) { bv = fabs(R[(size_t)kb * n + i]); ib = i; }
    sg = R[(size_t)kb * n + ib] < 0.0 ? -1.0 : 1.0;
    { double *cen = (double *)calloc((size_t)n, sizeof(double)), sum = 0.0, tot = 0.0;
      if (cen) {
        for (k = 0; k < P; k++) for (i = 0; i < n; i++) cen[i] += Y[(size_t)2 * k * n + i] / P;
        for (k = 0; k < P; k++) { double d = 0.0; for (i = 0; i < n; i++) d += R[(size_t)k * n + i] * (Y[(size_t)2 * k * n + i] - cen[i]); sum += d; tot += fabs(d); }
        if (tot > 0.0 && fabs(sum) > 1.e-6 * tot) sg = sum < 0.0 ? -1.0 : 1.0;
        free(cen); } }
    for (i = 0; i < P * n; i++) R[i] *= sg / s;
    for (i = 0; i < n; i++) d0 += Iv[i] * R[i];
    if (!(fabs(d0) > 1.e-14)) { ph->isoSkip = 10; goto out; }
    for (i = 0; i < P * n; i++) Iv[i] /= d0;
  }
  { double errDot = 0.0, errZR = 0.0, errIU = 0.0, errPerR = 0.0, errPerI = 0.0, imax = 0.0;
    for (k = 0; k < P; k++) {
      double dd = 0.0, zr = 0.0, iu = 0.0, in2 = 0.0, fn2 = 0.0;
      for (i = 0; i < n; i++) {
        double Ik = Iv[(size_t)k * n + i], Fk = Fv[(size_t)2 * k * n + i];
        dd += Ik * R[(size_t)k * n + i]; zr += Zs[(size_t)k * n + i] * R[(size_t)k * n + i]; iu += Ik * Fk; in2 += Ik * Ik; fn2 += Fk * Fk;
      }
      errDot = fmax(errDot, fabs(dd - 1.0)); errZR = fmax(errZR, fabs(zr));
      errIU = fmax(errIU, fabs(iu) / (sqrt(in2 * fn2) > 0.0 ? sqrt(in2 * fn2) : 1.0)); imax = fmax(imax, sqrt(in2));
    }
    for (i = 0; i < n; i++) { errPerR = fmax(errPerR, fabs(R[i] - R[(size_t)(P - 1) * n + i])); errPerI = fmax(errPerI, fabs(Iv[i] - Iv[(size_t)(P - 1) * n + i])); }
    errPerI /= (imax > 0.0 ? imax : 1.0);
    ph->nu = nu; ph->kappa = kappa; ph->nuRatio = (i2 >= 0) ? exp(-(rate2 - rate1)) : 0.0;
    /* (the ratio of the moduli of the next multiplier to the first, as Studio reports: mod(rest[1])/mod(e1)) */
    if (i2 >= 0) ph->nuRatio = hypot(ev[2 * i2], ev[2 * i2 + 1]) / mod1;
    ph->errDot = errDot; ph->errZR = errZR; ph->errIU = errIU; ph->errPerR = errPerR; ph->errPerI = errPerI;
  }
  /* B: the periodic solution of  B' = -(A^T + kappa) B - (dA)^T Z,  dA = d/ds DF(U + s rho) */
  if (n <= PHASE_BMAXN) {
    double *hk = (double *)malloc(sizeof(double) * (size_t)P * (size_t)n), *q = (double *)malloc(sizeof(double) * (size_t)(P - 1) * (size_t)n);
    double *scratch = (double *)malloc(sizeof(double) * (size_t)(4 * n + 4 * nn)), *dA = (double *)malloc(sizeof(double) * (size_t)nn);
    double *Phm = (double *)malloc(sizeof(double) * (size_t)nn), *zm = (double *)malloc(sizeof(double) * (size_t)n), *rm = (double *)malloc(sizeof(double) * (size_t)n);
    double *hm = (double *)malloc(sizeof(double) * (size_t)n), *qk = (double *)malloc(sizeof(double) * (size_t)n);
    if (hk && q && scratch && dA && Phm && zm && rm && hm && qk) {
      for (k = 0; k < P; k++) ph_hvec(n, Y + (size_t)2 * k * n, Zs + (size_t)k * n, R + (size_t)k * n, hk + (size_t)k * n, dA, scratch);
      for (k = 0; k < P - 1; k++) {
        const double *A0 = A + (size_t)2 * k * nn, *A2 = A + (size_t)2 * (k + 1) * nn, *Ps = PS + (size_t)k * nn;
        double e2 = exp(kappa * H / 2), e1k = exp(kappa * H);
        for (j = 0; j < n; j++) {     /* Hermite midpoints of Z and rho */
          double z0 = 0.0, z1 = 0.0, r0 = 0.0, r1 = 0.0;
          for (i = 0; i < n; i++) { z0 -= A0[i * n + j] * Zs[(size_t)k * n + i]; z1 -= A2[i * n + j] * Zs[(size_t)(k + 1) * n + i]; }
          for (i = 0; i < n; i++) { r0 += A0[j * n + i] * R[(size_t)k * n + i]; r1 += A2[j * n + i] * R[(size_t)(k + 1) * n + i]; }
          r0 -= kappa * R[(size_t)k * n + j]; r1 -= kappa * R[(size_t)(k + 1) * n + j];
          zm[j] = 0.5 * (Zs[(size_t)k * n + j] + Zs[(size_t)(k + 1) * n + j]) + H / 8 * (z0 - z1);
          rm[j] = 0.5 * (R[(size_t)k * n + j] + R[(size_t)(k + 1) * n + j]) + H / 8 * (r0 - r1);
        }
        ph_hvec(n, Y + (size_t)(2 * k + 1) * n, zm, rm, hm, dA, scratch);
        /* Hermite midpoint of the variational matrix: (I+Ps)/2 + H/8 (A_k - A_{k+1} Ps) */
        for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
          double a1 = 0.0; for (int l = 0; l < n; l++) a1 += A2[i * n + l] * Ps[l * n + j];
          Phm[i * n + j] = 0.5 * ((i == j ? 1.0 : 0.0) + Ps[i * n + j]) + H / 8 * (A0[i * n + j] - a1);
        }
        for (j = 0; j < n; j++) {
          double sm = 0.0, sp = 0.0;
          for (i = 0; i < n; i++) { sm += Phm[i * n + j] * hm[i]; sp += Ps[i * n + j] * hk[(size_t)(k + 1) * n + i]; }
          qk[j] = H / 6 * (hk[(size_t)k * n + j] + 4 * e2 * sm + e1k * sp);
        }
        memcpy(q + (size_t)k * n, qk, sizeof(double) * (size_t)n);
      }
      { double eK = exp(kappa * H);
        double *Bv = (double *)malloc(sizeof(double) * (size_t)P * (size_t)n), *c0 = (double *)malloc(sizeof(double) * (size_t)n), *Km = (double *)malloc(sizeof(double) * (size_t)nn);
        if (Bv && c0 && Km) {
#define PH_BACK(B0) { for (i = 0; i < n; i++) Bv[(size_t)(P - 1) * n + i] = (B0)[i]; \
            for (k = P - 2; k >= 0; k--) { const double *Ps_ = PS + (size_t)k * nn; \
              for (j = 0; j < n; j++) { double sm_ = 0.0; for (i = 0; i < n; i++) sm_ += Ps_[i * n + j] * Bv[(size_t)(k + 1) * n + i]; Bv[(size_t)k * n + j] = eK * sm_ + q[(size_t)k * n + j]; } } }
          for (i = 0; i < n; i++) wk[i] = 0.0;
          PH_BACK(wk);
          for (i = 0; i < n; i++) c0[i] = Bv[i];
          for (i = 0; i < n; i++) for (j = 0; j < n; j++) Km[i * n + j] = (i == j ? 1.0 : 0.0) - nu * Phi[j * n + i];
          if (ph_solve(n, Km, c0)) {      /* c0 now holds B(0) */
            PH_BACK(c0);
            { double ePer = 0.0, bm = 0.0; int fin = 1;
              for (i = 0; i < n; i++) ePer = fmax(ePer, fabs(Bv[i] - c0[i]));
              for (k = 0; k < P * n; k++) { bm = fmax(bm, fabs(Bv[k])); if (!isfinite(Bv[k])) fin = 0; }
              ph->errPerB = ePer / (bm > 0.0 ? bm : 1.0);
              if (fin) { Bs = Bv; Bv = NULL; ph->hasB = 1; } }
          } else ph->bNote = 2;
#undef PH_BACK
        }
        free(Bv); free(c0); free(Km);
      }
    }
    free(hk); free(q); free(scratch); free(dA); free(Phm); free(zm); free(rm); free(hm); free(qk);
  } else ph->bNote = 1;
  ph->iso = 1; ph->I = Iv; ph->R = R; Iv = NULL; R = NULL;
  if (Bs) { ph->B = Bs; Bs = NULL; }
  ok = 1;
out:
  free(cm); free(ev); free(wk); free(vr); free(wl); free(r); free(R); free(Iv); free(Bs);
  return ok;
}

/* one try on a grid of about N steps. 1 = Z found (the checks are in ph->errNorm, ph->errPer) */
static int ph_on_grid(int n, const double *x0in, double Tin, int N, PHASE *ph)
{
  int nn = n * n, i, j, k, P, it, kk, ok = 0;
  double *x0 = (double *)malloc(sizeof(double) * (size_t)n), T = Tin, closure0 = 0.0, closure1 = -1.0, sc0 = 1.0;
  PHGRID G; double *PS = NULL, *Phi = (double *)malloc(sizeof(double) * (size_t)nn);
  PHGRID Gbest; double *PSbest = NULL; double Tbest = T, rmbest = 1.e300; int havebest = 0;
  memcpy(x0, x0in, sizeof(double) * (size_t)n);
  memset(&G, 0, sizeof G); memset(&Gbest, 0, sizeof Gbest);
  if (!x0 || !Phi) goto out;
  if (!ph_grid(n, x0, T, N, &G)) goto out;
  P = G.P;
  PS = (double *)malloc(sizeof(double) * (size_t)(P - 1) * (size_t)nn);
  if (!PS || !ph_stepmats(n, &G, PS, Phi)) goto out;

  /* make the orbit exactly periodic on the grid: Newton on (x, T) for the grid's own map, with the phase condition
     f(x0).(x - x0) = 0. The best iterate is kept. */
  for (i = 0; i < n; i++) { closure0 = fmax(closure0, fabs(G.end[i] - x0[i])); sc0 = fmax(sc0, fabs(x0[i])); }
  if (closure0 > 1.e-11 * sc0 && n <= 40) {
    int M1 = n + 1;
    double *Fref = (double *)malloc(sizeof(double) * (size_t)n), *xp = (double *)malloc(sizeof(double) * (size_t)n);
    double *r = (double *)malloc(sizeof(double) * (size_t)M1), *J = (double *)malloc(sizeof(double) * (size_t)M1 * (size_t)M1);
    double *xt = (double *)malloc(sizeof(double) * (size_t)n), *e2 = (double *)malloc(sizeof(double) * (size_t)n);
    double *Fe = (double *)malloc(sizeof(double) * (size_t)n), *Ph2 = (double *)malloc(sizeof(double) * (size_t)nn);
    double *PhiBest = (double *)malloc(sizeof(double) * (size_t)nn), *dv = (double *)malloc(sizeof(double) * (size_t)M1);
    double Tp = T, rm;
    if (Fref && xp && r && J && xt && e2 && Fe && Ph2 && PhiBest && dv) {
      ph_f(n, x0, Fref);
      memcpy(xp, x0, sizeof(double) * (size_t)n);
      for (it = 0; it < 12; it++) {
        PHGRID Gs; double *PSs;
        memset(&Gs, 0, sizeof Gs);
        if (!ph_grid(n, xp, Tp, N, &Gs)) break;
        PSs = (double *)malloc(sizeof(double) * (size_t)(Gs.P - 1) * (size_t)nn);
        if (!PSs || !ph_stepmats(n, &Gs, PSs, Ph2)) { free(PSs); ph_gfree(&Gs); break; }
        rm = 0.0;
        for (i = 0; i < n; i++) { r[i] = Gs.end[i] - xp[i]; rm = fmax(rm, fabs(r[i])); }
        r[n] = 0.0; for (i = 0; i < n; i++) r[n] += Fref[i] * (xp[i] - x0[i]);
        ph_f(n, Gs.end, Fe);                       /* d(end state)/dT */
        if (rm < rmbest) {
          if (havebest) { ph_gfree(&Gbest); free(PSbest); }
          Gbest = Gs; PSbest = PSs; rmbest = rm; Tbest = Tp; havebest = 1;
          memcpy(PhiBest, Ph2, sizeof(double) * (size_t)nn);
        } else { ph_gfree(&Gs); free(PSs); }
        if (rm < 1.e-13 * sc0) break;
        for (i = 0; i < n; i++) { for (j = 0; j < n; j++) J[i * M1 + j] = Ph2[i * n + j] - (i == j ? 1.0 : 0.0); J[i * M1 + n] = Fe[i]; }
        for (j = 0; j < n; j++) J[n * M1 + j] = Fref[j];
        J[n * M1 + n] = 0.0;
        for (i = 0; i <= n; i++) dv[i] = -r[i];
        if (!ph_solve(M1, J, dv)) break;
        { double lam = 1.0; int moved = 0;           /* a damped step: halve it until the closure residual falls */
          for (kk = 0; kk < 6; kk++) {
            double Tt = Tp + lam * dv[n], rnew = 0.0;
            for (i = 0; i < n; i++) xt[i] = xp[i] + lam * dv[i];
            if (Tt > 0.0) {
              ph_endstate(n, xt, Tt, N, e2);
              for (i = 0; i < n; i++) rnew = fmax(rnew, fabs(e2[i] - xt[i]));
              if (rnew < rm) { memcpy(xp, xt, sizeof(double) * (size_t)n); Tp = Tt; moved = 1; break; }
            }
            lam *= 0.5;
          }
          if (!moved) break; }
      }
    }
    free(Fref); free(xp); free(r); free(J); free(xt); free(e2); free(Fe); free(Ph2); free(dv);
    if (havebest && rmbest < 0.1 * closure0 && PhiBest) {
      ph_gfree(&G); free(PS); G = Gbest; PS = PSbest; T = Tbest; closure1 = rmbest;
      memcpy(Phi, PhiBest, sizeof(double) * (size_t)nn);
      memset(&Gbest, 0, sizeof Gbest); PSbest = NULL; havebest = 0;
    }
    free(PhiBest);
  }
  P = G.P;

  /* the left eigenvector of Phi for the multiplier 1: Phi^T z = z, by inverse iteration with a shift just off 1 */
  {
    static const double shifts[3] = {1.e-9, 1.e-6, 1.e-4};
    double *z = (double *)malloc(sizeof(double) * (size_t)n), *w = (double *)malloc(sizeof(double) * (size_t)n), *B = (double *)malloc(sizeof(double) * (size_t)nn);
    int found = 0, s;
    for (s = 0; s < 3 && !found; s++) {
      int good = 1;
      for (i = 0; i < n; i++) w[i] = 1.0;
      for (it = 0; it < 6 && good; it++) {
        double nm = 0.0;
        for (i = 0; i < n; i++) for (j = 0; j < n; j++) B[i * n + j] = Phi[j * n + i] - (i == j ? 1.0 + shifts[s] : 0.0);
        memcpy(z, w, sizeof(double) * (size_t)n);
        if (!ph_solve(n, B, z)) { good = 0; break; }
        for (i = 0; i < n; i++) nm = fmax(nm, fabs(z[i]));
        if (!(nm > 0.0) || !isfinite(nm)) { good = 0; break; }
        for (i = 0; i < n; i++) w[i] = z[i] / nm;
      }
      if (good) { memcpy(z, w, sizeof(double) * (size_t)n); found = 1; }
    }
    if (found) {
      double sc = 0.0, *Zs = (double *)malloc(sizeof(double) * (size_t)P * (size_t)n), *zc = (double *)malloc(sizeof(double) * (size_t)n), *zn = (double *)malloc(sizeof(double) * (size_t)n);
      for (i = 0; i < n; i++) sc += z[i] * G.F[i];
      if (Zs && zc && zn && fabs(sc) > 1.e-14) {
        for (i = 0; i < n; i++) z[i] /= sc;
        memcpy(Zs + (size_t)(P - 1) * n, z, sizeof(double) * (size_t)n); memcpy(zc, z, sizeof(double) * (size_t)n);
        for (k = P - 2; k >= 0; k--) {
          const double *Ps = PS + (size_t)k * nn;
          for (i = 0; i < n; i++) { double sum = 0.0; for (j = 0; j < n; j++) sum += Ps[j * n + i] * zc[j]; zn[i] = sum; }
          memcpy(zc, zn, sizeof(double) * (size_t)n); memcpy(Zs + (size_t)k * n, zc, sizeof(double) * (size_t)n);
        }
        {
          double errNorm = 0.0, errPer = 0.0, zmax = 0.0;
          int fin = 1;
          for (i = 0; i < P * n; i++) if (!isfinite(Zs[i])) fin = 0;
          if (fin) {
            for (k = 0; k < P; k++) { double d = 0.0; for (i = 0; i < n; i++) d += Zs[(size_t)k * n + i] * G.F[(size_t)2 * k * n + i]; errNorm = fmax(errNorm, fabs(d - 1.0)); }
            for (i = 0; i < n; i++) { errPer = fmax(errPer, fabs(Zs[i] - z[i])); zmax = fmax(zmax, fabs(z[i])); }
            errPer /= (zmax > 0.0 ? zmax : 1.0);
            ph->n = n; ph->P = P; ph->T = T; ph->Ng = G.Ng; ph->closure0 = closure0; ph->closure1 = closure1;
            ph->errNorm = errNorm; ph->errPer = errPer;
            ph->th = (double *)malloc(sizeof(double) * (size_t)P); ph->U = (double *)malloc(sizeof(double) * (size_t)P * (size_t)n);
            ph->F = (double *)malloc(sizeof(double) * (size_t)P * (size_t)n); ph->Z = Zs;
            if (ph->th && ph->U && ph->F) {
              for (k = 0; k < P; k++) {
                ph->th[k] = (double)k / (double)(P - 1);
                memcpy(ph->U + (size_t)k * n, G.Y + (size_t)2 * k * n, sizeof(double) * (size_t)n);
                memcpy(ph->F + (size_t)k * n, G.F + (size_t)2 * k * n, sizeof(double) * (size_t)n);
              }
              ph_isostable(n, &G, PS, Phi, ph->Z, T, ph);
              Zs = NULL; ok = 1;
            }
          }
        }
      }
      free(Zs); free(zc); free(zn);
    }
    free(z); free(w); free(B);
  }
out:
  free(x0); free(Phi); free(PS); ph_gfree(&G);
  if (havebest) { ph_gfree(&Gbest); free(PSbest); }
  return ok;
}

/* If T is m times the period of the cycle through x (m = 2..8), return the largest such m for which the orbit
   closes after T/m (so that the smallest period is found); otherwise 1. */
static int ph_multiple(int n, const double *x, double T, int N)
{
  int m, i;
  double *y = (double *)malloc(sizeof(double) * (size_t)n), sc = 1.0;
  if (!y) return 1;
  for (i = 0; i < n; i++) sc = fmax(sc, fabs(x[i]));
  for (m = 8; m >= 2; m--) {
    double d = 0.0;
    ph_endstate(n, x, T / m, N, y);
    for (i = 0; i < n; i++) d = fmax(d, fabs(y[i] - x[i]));
    if (d < 1.e-6 * sc) { free(y); return m; }
  }
  free(y);
  return 1;
}

/* The adjoint of the cycle through x0 with period T. A grid of about N0 steps (default 800) is tried first and refined
   (doubled) while the checks improve, up to 6400. If T is a whole number of times the true period, the adjoint is for the
   true period (ph->multiple says how many times). The orbit must not be a steady state. Returns 1 and fills *ph. */
static int ph_adjoint_once(int n, const double *x0, double T, int N0, PHASE *ph)
{
  int N, ok = 0;
  double prev = 1.e300;
  PHASE best, cur;
  memset(&best, 0, sizeof best); memset(ph, 0, sizeof *ph);
  for (N = N0 > 0 ? N0 : 800; N <= 6400; N *= 2) {
    double e;
    memset(&cur, 0, sizeof cur);
    if (!ph_on_grid(n, x0, T, N, &cur)) { if (best.th) break; continue; }
    e = cur.errNorm + cur.errPer;
    if (!best.th || e < best.errNorm + best.errPer) { phase_free(&best); best = cur; } else phase_free(&cur);
    if ((best.errNorm < 2.e-5 && best.errPer < 2.e-5) || e > 0.7 * prev) break;
    prev = e;
  }
  if (best.th) { *ph = best; ok = 1; }
  return ok;
}

int phase_adjoint(int n, const double *x0, double T, int N0, PHASE *ph)
{
  int m, k, i;
  double ext = 0.0, sc = 1.0;
  memset(ph, 0, sizeof *ph);
  if (n < 1 || n > PHASE_MAXN || !(T > 0.0)) return 0;
  if (!ph_adjoint_once(n, x0, T, N0, ph)) return 0;
  /* a steady state is not a cycle */
  for (k = 0; k < ph->P; k++) for (i = 0; i < n; i++) { ext = fmax(ext, fabs(ph->U[(size_t)k * n + i] - ph->U[i])); sc = fmax(sc, fabs(ph->U[i])); }
  if (ext < 1.e-6 * sc) { phase_free(ph); return 0; }
  /* the period given may be a multiple of the true one: then use the true one */
  m = ph_multiple(n, ph->U, ph->T, ph->Ng);
  if (m > 1) {
    PHASE again; double *xs = (double *)malloc(sizeof(double) * (size_t)n), Tn = ph->T / m;
    if (xs) {
      memcpy(xs, ph->U, sizeof(double) * (size_t)n);
      if (ph_adjoint_once(n, xs, Tn, N0, &again)) { phase_free(ph); *ph = again; ph->multiple = m; }
      free(xs);
    }
  }
  return 1;
}


/* ------------------------------------------------------------------------------------------------------------------
   Stage 3: H and H2 (as in XPP Studio). h2(t,s) = Z(t).[G_1(U(t),U(s)) rho(t)] + B(t).G(U(t),U(s)),
   h3(t,s) = Z(t).[G_2(U(t),U(s)) rho(s)], G_1 and G_2 the derivatives of G with respect to its two arguments, taken along rho
   by central differences. f_phi(t) = integral_0^inf exp(kappa s) I(t-s).G(U(t-s),U(t-s+phi)) ds is the periodic solution of
   f' = kappa f + I(t).G(U(t),U(t+phi)), found exactly by Fourier series (the coefficient of exp(i w t) is divided by i w - kappa).
   ------------------------------------------------------------------------------------------------------------------ */

/* derivative at 0 of the periodic function sampled at j T/Mp, by Fourier series */
static double ph_spec_d0(const double *a, int Mp, double T)
{
  int m, j;
  double s = 0.0;
  for (m = 1; m < Mp / 2; m++) {
    double im = 0.0;
    for (j = 0; j < Mp; j++) im -= a[j] * sin(2.0 * M_PI * m * j / Mp);
    im /= Mp;
    s += -2.0 * m * im;
  }
  return s * 2.0 * M_PI / T;
}

void phase_h_free(PHASEH *h)
{
  free(h->H); free(h->H2); free(h->A); free(h->Bf); free(h->g1); free(h->g2); free(h->g3); free(h->g23);
  memset(h, 0, sizeof *h);
}

int phase_hfun(const PHASE *ph, const int *idx, int nc, PH_COUPLE G, void *ctx, PHASEH *out)
{
  int n, L, Mp = PHASE_HM, st, k, j, q, i, m, hf, ok = 0;
  double T, kap, hh, umax = 0.0;
  double *Us, *Up, *Um, *Rs, *Q, *H2a, *H3a, *Hs, *F, *o0, *op, *om, *oq, *g0, *d1, *cs, *sn, *qr, *qi, *Ik, *Zk, *Bk;
  memset(out, 0, sizeof *out);
  if (!ph || !ph->iso || !ph->hasB || nc < 1) return 0;
  n = ph->n; L = ph->P - 1; T = ph->T; kap = ph->kappa;
  if (L % Mp != 0) return 0;
  st = L / Mp; hf = L / 2;
  Us = (double *)malloc(sizeof(double) * (size_t)L * n); Up = (double *)malloc(sizeof(double) * (size_t)L * n);
  Um = (double *)malloc(sizeof(double) * (size_t)L * n); Rs = (double *)malloc(sizeof(double) * (size_t)L * n);
  Q = (double *)calloc((size_t)Mp * L, sizeof(double)); H2a = (double *)calloc((size_t)Mp * L, sizeof(double));
  H3a = (double *)calloc((size_t)Mp * L, sizeof(double)); Hs = (double *)calloc((size_t)Mp + 1, sizeof(double));
  F = (double *)calloc((size_t)Mp * L, sizeof(double));
  o0 = (double *)malloc(sizeof(double) * 4 * (size_t)n); op = o0 + n; om = op + n; oq = om + n;
  g0 = (double *)malloc(sizeof(double) * (size_t)Mp * nc); d1 = (double *)malloc(sizeof(double) * (size_t)Mp * nc);
  cs = (double *)malloc(sizeof(double) * (size_t)L); sn = (double *)malloc(sizeof(double) * (size_t)L);
  qr = (double *)malloc(sizeof(double) * (size_t)hf); qi = (double *)malloc(sizeof(double) * (size_t)hf);
  Ik = (double *)malloc(sizeof(double) * (size_t)nc * 3); Zk = Ik + nc; Bk = Zk + nc;
  if (!Us || !Up || !Um || !Rs || !Q || !H2a || !H3a || !Hs || !F || !o0 || !g0 || !d1 || !cs || !sn || !qr || !qi || !Ik) goto done;
  for (k = 0; k < L; k++) for (i = 0; i < n; i++) {
    Us[(size_t)k * n + i] = ph->U[(size_t)k * n + i]; Rs[(size_t)k * n + i] = ph->R[(size_t)k * n + i];
    umax = fmax(umax, fabs(Us[(size_t)k * n + i]));
  }
  hh = 1.e-4 * (1.0 + umax);
  for (k = 0; k < L; k++) for (i = 0; i < n; i++) {
    Up[(size_t)k * n + i] = Us[(size_t)k * n + i] + hh * Rs[(size_t)k * n + i];
    Um[(size_t)k * n + i] = Us[(size_t)k * n + i] - hh * Rs[(size_t)k * n + i];
  }
  for (k = 0; k < L; k++) {
    for (q = 0; q < nc; q++) { Ik[q] = ph->I[(size_t)k * n + idx[q]]; Zk[q] = ph->Z[(size_t)k * n + idx[q]]; Bk[q] = ph->B[(size_t)k * n + idx[q]]; }
    for (j = 0; j < Mp; j++) {
      double qv = 0.0, h2 = 0.0, h3 = 0.0, hs = 0.0;
      m = (k + j * st) % L;
      G(Us + (size_t)k * n, Us + (size_t)m * n, o0, ctx);
      for (q = 0; q < nc; q++) g0[q] = o0[idx[q]];
      G(Up + (size_t)k * n, Us + (size_t)m * n, op, ctx);
      G(Um + (size_t)k * n, Us + (size_t)m * n, om, ctx);
      for (q = 0; q < nc; q++) d1[q] = (op[idx[q]] - om[idx[q]]) / (2.0 * hh);
      G(Us + (size_t)k * n, Up + (size_t)m * n, op, ctx);
      G(Us + (size_t)k * n, Um + (size_t)m * n, om, ctx);
      for (q = 0; q < nc; q++) {
        qv += Ik[q] * g0[q]; hs += Zk[q] * g0[q];
        h2 += Zk[q] * d1[q] + Bk[q] * g0[q];
        h3 += Zk[q] * (op[idx[q]] - om[idx[q]]) / (2.0 * hh);
      }
      Q[(size_t)j * L + k] = qv; H2a[(size_t)j * L + k] = h2; H3a[(size_t)j * L + k] = h3; Hs[j] += hs;
    }
  }
  (void)oq;
  /* f_phi = the periodic solution of f' = kappa f + q_phi, by Fourier series (the Nyquist term is dropped) */
  for (k = 0; k < L; k++) { cs[k] = cos(2.0 * M_PI * k / L); sn[k] = sin(2.0 * M_PI * k / L); }
  for (j = 0; j < Mp; j++) {
    const double *qj = Q + (size_t)j * L; double *f = F + (size_t)j * L, q0 = 0.0;
    for (k = 0; k < L; k++) q0 += qj[k];
    q0 /= L;
    for (m = 1; m < hf; m++) {
      double a = 0.0, b = 0.0;
      for (k = 0; k < L; k++) { int t = (int)(((long)m * k) % L); a += qj[k] * cs[t]; b -= qj[k] * sn[t]; }
      qr[m] = a / L; qi[m] = b / L;
    }
    for (k = 0; k < L; k++) f[k] = q0 / (-kap);
    for (m = 1; m < hf; m++) {
      double w = 2.0 * M_PI * m / T, dr = -kap, di = w, den = dr * dr + di * di;
      double gr = (qr[m] * dr + qi[m] * di) / den, gi = (qi[m] * dr - qr[m] * di) / den;
      for (k = 0; k < L; k++) { int t = (int)(((long)m * k) % L); f[k] += 2.0 * (gr * cs[t] - gi * sn[t]); }
    }
  }
  out->Mp = Mp; out->T = T;
  out->H = (double *)calloc((size_t)Mp + 1, sizeof(double)); out->H2 = (double *)calloc((size_t)Mp + 1, sizeof(double));
  out->A = (double *)calloc((size_t)Mp + 1, sizeof(double)); out->Bf = (double *)calloc((size_t)Mp + 1, sizeof(double));
  out->g1 = (double *)calloc((size_t)Mp + 1, sizeof(double)); out->g2 = (double *)calloc((size_t)Mp + 1, sizeof(double));
  out->g3 = (double *)calloc((size_t)Mp + 1, sizeof(double)); out->g23 = (double *)calloc((size_t)Mp + 1, sizeof(double));
  if (!out->H || !out->H2 || !out->A || !out->Bf || !out->g1 || !out->g2 || !out->g3 || !out->g23) { phase_h_free(out); goto done; }
  for (j = 0; j < Mp; j++) {
    int jm = (Mp - j) % Mp; double a = 0.0, b = 0.0;
    for (k = 0; k < L; k++) { a += F[(size_t)j * L + k] * H2a[(size_t)j * L + k]; b += F[(size_t)jm * L + (k + j * st) % L] * H3a[(size_t)j * L + k]; }
    out->A[j] = a / L; out->Bf[j] = b / L; out->H[j] = Hs[j] / L;
  }
  out->A[Mp] = out->A[0]; out->Bf[Mp] = out->Bf[0]; out->H[Mp] = out->H[0];
  for (j = 0; j <= Mp; j++) {
    int jj = j % Mp, jm = (Mp - jj) % Mp;
    out->g1[j] = out->H[jm] - out->H[jj]; out->g2[j] = out->A[jm] - out->A[jj]; out->g3[j] = out->Bf[jm] - out->Bf[jj];
    out->g23[j] = out->g2[j] + out->g3[j]; out->H2[j] = out->A[j] + out->Bf[j];
  }
  out->gammaA = -ph_spec_d0(out->A, Mp, T); out->gammaB = -ph_spec_d0(out->Bf, Mp, T); out->gamma = ph_spec_d0(out->g23, Mp, T) / 4.0;
  if (fabs(kap) * T < 0.3) out->notes |= 1;
  if (ph->nuRatio > 0.25) out->notes |= 2;
  ok = 1;
done:
  free(Us); free(Up); free(Um); free(Rs); free(Q); free(H2a); free(H3a); free(Hs); free(F); free(o0); free(g0); free(d1);
  free(cs); free(sn); free(qr); free(qi); free(Ik);
  return ok;
}
