/*
 *  periodic_ic.c        (I)nitialconds (P)eriodic
 *
 *  Takes the trajectory in storage, which should have settled onto a stable
 *  limit cycle, and
 *
 *   1. finds the period and the periodic orbit in double precision. The stored
 *      solution is single precision, so it is only used for a first estimate:
 *      a section through a point of the last cycle, and the time to return to it.
 *      The estimate is then improved by following the return map to the section
 *      (integrating from the section to the next crossing of it) until it stops
 *      changing. A stable cycle attracts, so this converges geometrically and
 *      needs no Jacobian, whatever the number of equations.
 *   2. finds the point of the orbit where the variable on the x-axis of the
 *      current view (in a time plot: the variable that is plotted) is largest,
 *   3. makes that point the initial condition, sets the total time (Numerics)
 *      to the period and the transient to 0.
 *
 *  periodic_ic_prepare() does steps 1-3 and returns 1 if it succeeded. The
 *  caller (do_init_data in integrate.c) then runs the equations from there as
 *  for (G)o, so one period is drawn, starting at the peak.
 *
 *  It refuses systems whose solutions are not simple trajectories of an
 *  autonomous ODE: maps, delays, Markov variables, events, Volterra
 *  equations, explicit dependence on t, and a Poincare map.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <X11/Xlib.h>
#include "xpplim.h"
#include "struct.h"

extern float **storage;
extern int storind, NODE, NEQ, METHOD, NJMP, INFLAG, POIMAP, DelayFlag;
extern int NMarkov, NFlags, TORUS, itor[MAXODE];
extern double TEND, TRANS, DELTA_T, T0, last_ic[MAXODE];
extern char uvar_names[MAXODE][12];
extern GRAPH *MyGraph;
extern int (*rhs)();
void err_msg(char *string);
int plintf(char *fmt, ...);
void ping(void);
int script_is_active(void);

#define PIC_VOLTERRA 6
#define PIC_MAXSTEPS 4000000L
#define PIC_NSAMPLE 5000
#define PIC_MAXCYCLES 400
#define PIC_BUDGET 40.0           /* seconds */

/* ---------------------------------------------------------------------- */
/*   an adaptive Dormand-Prince 5(4) integrator for  x' = f(x)             */
/* ---------------------------------------------------------------------- */

typedef struct {
  int n;
  double *k1, *k2, *k3, *k4, *k5, *k6, *k7, *y, *xn;
  double rtol, atol;
  long nfev;
} PICW;

static int pic_alloc(PICW *w, int n)
{
  double **p[9];
  int i;
  p[0] = &w->k1; p[1] = &w->k2; p[2] = &w->k3; p[3] = &w->k4; p[4] = &w->k5;
  p[5] = &w->k6; p[6] = &w->k7; p[7] = &w->y;  p[8] = &w->xn;
  w->n = n; w->nfev = 0; w->rtol = 1.e-11; w->atol = 1.e-13;
  for (i = 0; i < 9; i++) {
    *p[i] = (double *)calloc((size_t)n, sizeof(double));
    if (*p[i] == NULL) return 0;
  }
  return 1;
}

static void pic_free(PICW *w)
{
  free(w->k1); free(w->k2); free(w->k3); free(w->k4); free(w->k5);
  free(w->k6); free(w->k7); free(w->y);  free(w->xn);
}

static void pic_f(PICW *w, const double *x, double *f)
{
  (*rhs)(0.0, (double *)x, f, w->n);
  w->nfev++;
}

/* one step of size h from x to xn; returns the error in units of the tolerance */
static double pic_step(PICW *w, const double *x, double h, double *xn)
{
  int i, n = w->n;
  double *y = w->y, err = 0.0, e, sc, a, b;
  double *k1 = w->k1, *k2 = w->k2, *k3 = w->k3, *k4 = w->k4, *k5 = w->k5, *k6 = w->k6, *k7 = w->k7;

  pic_f(w, x, k1);
  for (i = 0; i < n; i++) y[i] = x[i] + h * (1. / 5 * k1[i]);
  pic_f(w, y, k2);
  for (i = 0; i < n; i++) y[i] = x[i] + h * (3. / 40 * k1[i] + 9. / 40 * k2[i]);
  pic_f(w, y, k3);
  for (i = 0; i < n; i++) y[i] = x[i] + h * (44. / 45 * k1[i] - 56. / 15 * k2[i] + 32. / 9 * k3[i]);
  pic_f(w, y, k4);
  for (i = 0; i < n; i++)
    y[i] = x[i] + h * (19372. / 6561 * k1[i] - 25360. / 2187 * k2[i] + 64448. / 6561 * k3[i] - 212. / 729 * k4[i]);
  pic_f(w, y, k5);
  for (i = 0; i < n; i++)
    y[i] = x[i] + h * (9017. / 3168 * k1[i] - 355. / 33 * k2[i] + 46732. / 5247 * k3[i] + 49. / 176 * k4[i]
                       - 5103. / 18656 * k5[i]);
  pic_f(w, y, k6);
  for (i = 0; i < n; i++)
    xn[i] = x[i] + h * (35. / 384 * k1[i] + 500. / 1113 * k3[i] + 125. / 192 * k4[i] - 2187. / 6784 * k5[i]
                        + 11. / 84 * k6[i]);
  pic_f(w, xn, k7);
  for (i = 0; i < n; i++) {
    e = h * (71. / 57600 * k1[i] - 71. / 16695 * k3[i] + 71. / 1920 * k4[i] - 17253. / 339200 * k5[i]
             + 22. / 525 * k6[i] - 1. / 40 * k7[i]);
    a = fabs(x[i]); b = fabs(xn[i]);
    sc = w->atol + w->rtol * (a > b ? a : b);
    e = fabs(e) / sc;
    if (!(e == e)) return 1.e300;
    if (e > err) err = e;
  }
  return err;
}

static double pic_newh(double h, double err)
{
  double fac = (err <= 1.e-12) ? 5.0 : 0.9 * pow(err, -0.2);
  if (fac > 5.0) fac = 5.0;
  if (fac < 0.2) fac = 0.2;
  return h * fac;
}

/* integrate x over a time T (x is changed in place). 0 = ok */
static int pic_to(PICW *w, double *x, double T, double *hg)
{
  double t = 0.0, h = *hg, hh, err;
  long steps = 0;
  int n = w->n, i, last;
  if (T <= 0.0) return 0;
  if (!(h > 0.0) || h > T) h = T / 10.0;
  while (t < T) {
    last = 0; hh = h;
    if (t + hh >= T * (1.0 - 1.e-14)) { hh = T - t; last = 1; }
    err = pic_step(w, x, hh, w->xn);
    if (err > 1.e299) return -1;
    if (++steps > PIC_MAXSTEPS) return -2;
    if (err <= 1.0) {
      for (i = 0; i < n; i++) x[i] = w->xn[i];
      t = last ? T : t + hh;
      if (!last) h = pic_newh(h, err);
    } else {
      h = pic_newh(hh, err);
      if (h < 1.e-14 * (T > 1.0 ? T : 1.0)) return -3;
    }
  }
  *hg = h;
  return 0;
}

/* integrate from x, which lies on the section  fs.(x-xs)=0, until the flow comes back to that
   plane (crossing it from the negative to the positive side) at a point near xs.
   x becomes the crossing point and *Tc the time taken.  0 = ok */
static double pic_h(int n, const double *x, const double *xs, const double *fs)
{
  int i;
  double s = 0.0;
  for (i = 0; i < n; i++) s += fs[i] * (x[i] - xs[i]);
  return s;
}

static int pic_cross(PICW *w, double *x, double Tmin, double Tmax, const double *xs, const double *fs,
                     const double *scl, double *Tc, double *hg, double *xold)
{
  int n = w->n, i, armed = 0, it;
  double t = 0.0, h = *hg, err, hprev, hnew, dist, d, lo, hi, mid, hm;
  long steps = 0;

  if (!(h > 0.0)) h = Tmax / 100.0;
  hprev = pic_h(n, x, xs, fs);
  while (t < Tmax) {
    for (i = 0; i < n; i++) xold[i] = x[i];
    err = pic_step(w, xold, h, w->xn);
    if (err > 1.e299) return -1;
    if (++steps > PIC_MAXSTEPS) return -2;
    if (err > 1.0) {
      h = pic_newh(h, err);
      if (h < 1.e-14 * Tmax) return -3;
      continue;
    }
    hnew = pic_h(n, w->xn, xs, fs);
    dist = 0.0;
    for (i = 0; i < n; i++) { d = fabs(w->xn[i] - xs[i]) / scl[i]; if (d > dist) dist = d; }
    if (!armed) {
      if (hnew < 0.0) armed = 1;
    } else if (hprev < 0.0 && hnew >= 0.0 && dist < 0.25 && t + h >= Tmin) {
      /* the crossing lies in this step: bisect on the step length */
      lo = 0.0; hi = h;
      for (it = 0; it < 70; it++) {
        mid = 0.5 * (lo + hi);
        pic_step(w, xold, mid, x);
        hm = pic_h(n, x, xs, fs);
        if (hm < 0.0) lo = mid; else hi = mid;
      }
      pic_step(w, xold, hi, x);
      *Tc = t + hi;
      *hg = h;
      return 0;
    }
    for (i = 0; i < n; i++) x[i] = w->xn[i];
    t += h;
    hprev = hnew;
    h = pic_newh(h, err);
  }
  return -4;
}

/* ---------------------------------------------------------------------- */

static void pic_fail(const char *msg, double *a, double *b, double *c, double *d, double *e, double *f, PICW *w)
{
  free(a); free(b); free(c); free(d); free(e); free(f);
  if (w) pic_free(w);
  err_msg((char *)msg);
}

/* Hermite cubic through (t0,x0,f0) and (t1,x1,f1) at the fraction th of the way */
static void pic_hermite(int n, double th, double h, const double *x0, const double *f0,
                        const double *x1, const double *f1, double *out)
{
  int i;
  double h00 = (1 + 2 * th) * (1 - th) * (1 - th), h10 = th * (1 - th) * (1 - th);
  double h01 = th * th * (3 - 2 * th), h11 = th * th * (th - 1);
  for (i = 0; i < n; i++) out[i] = h00 * x0[i] + h10 * h * f0[i] + h01 * x1[i] + h11 * h * f1[i];
}

int periodic_ic_prepare(void)
{
  int n = NODE, N = storind, i, j, k, iv, far, ks, cyc, conv = 0, nbr, ib, kbest, ns;
  double *xe = NULL, *fe = NULL, *amp = NULL, *scl = NULL, *wk = NULL, *wk2 = NULL;
  double *xs, *fs, *x, *xo, *xp, *xj, *bestx;
  double *vs, *fv;
  double maxamp, d, th, ts, T0e, T, Tc, hg, tstart, dT, dx, closure, vbest;
  double ta, tb, tm, g, vk;
  PICW W;
  char msg[160];
  int *brk;

  /* ---- is there something we can work with? ---- */
  if (INFLAG == 0 || N < 2) { err_msg("No prior solution"); return 0; }
  if (METHOD == 0) { err_msg("Periodic: not available for maps"); return 0; }
  if (METHOD == PIC_VOLTERRA) { err_msg("Periodic: not available for Volterra equations"); return 0; }
  if (DelayFlag) { err_msg("Periodic: not available for delay equations"); return 0; }
  if (NMarkov > 0) { err_msg("Periodic: not available with Markov variables"); return 0; }
  if (NFlags > 0) { err_msg("Periodic: not available with global events (flags)"); return 0; }
  if (POIMAP) { err_msg("Periodic: the Poincare map is on, so the stored data are not a trajectory"); return 0; }
  if (DELTA_T < 0.0) { err_msg("Periodic: integrate forward in time first"); return 0; }
  if (n < 1 || n > NEQ) { err_msg("Periodic: no equations"); return 0; }
  if (TORUS) for (i = 0; i < n; i++) if (itor[i]) { err_msg("Periodic: not available with periodic (torus) variables"); return 0; }
  if (N < 30) { err_msg("Periodic: the stored trajectory is too short. Integrate longer first"); return 0; }
  if (!(storage[0][N - 1] > storage[0][0])) { err_msg("Periodic: integrate forward in time first"); return 0; }

  /* the variable that sets the phase: the x-axis variable, or in a time plot the plotted one */
  iv = MyGraph->xv[0] > 0 ? MyGraph->xv[0] - 1 : MyGraph->yv[0] - 1;
  if (iv < 0 || iv >= n) {
    err_msg("Periodic: put a state variable on the x-axis (or plot one against time)");
    return 0;
  }

  memset(&W, 0, sizeof(W));
  if (!pic_alloc(&W, n)) { pic_free(&W); err_msg("Periodic: out of memory"); return 0; }
  xe = (double *)calloc((size_t)n, sizeof(double)); fe = (double *)calloc((size_t)n, sizeof(double));
  amp = (double *)calloc((size_t)n, sizeof(double)); scl = (double *)calloc((size_t)n, sizeof(double));
  wk = (double *)calloc((size_t)(8 * n), sizeof(double));
  wk2 = (double *)calloc((size_t)(2 * n), sizeof(double));
  if (!xe || !fe || !amp || !scl || !wk || !wk2) { pic_fail("Periodic: out of memory", xe, fe, amp, scl, wk, wk2, &W); return 0; }
  xs = wk; fs = wk + n; x = wk + 2 * n; xo = wk + 3 * n; xp = wk + 4 * n; xj = wk + 5 * n; bestx = wk + 6 * n;

  /* ---- is the system autonomous? ---- */
  for (i = 0; i < n; i++) xe[i] = (double)storage[i + 1][N - 1];
  {
    double *fa = wk + 7 * n, *fb = wk2;
    (*rhs)(0.0, xe, fa, n);
    (*rhs)(1.2345, xe, fb, n);
    for (i = 0; i < n; i++)
      if (fabs(fa[i] - fb[i]) > 1.e-9 * (1.0 + fabs(fa[i]))) {
        pic_fail("Periodic: the equations depend on t explicitly, so there is no periodic orbit of this kind",
                 xe, fe, amp, scl, wk, wk2, &W);
        return 0;
      }
  }

  /* ---- first estimate from the stored (single precision) data ---- */
  pic_f(&W, xe, fe);
  maxamp = 0.0;
  for (i = 0; i < n; i++) {
    double lo = 1.e300, hi = -1.e300, v;
    for (k = N / 2; k < N; k++) { v = (double)storage[i + 1][k]; if (v < lo) lo = v; if (v > hi) hi = v; }
    amp[i] = hi - lo;
    if (amp[i] > maxamp) maxamp = amp[i];
  }
  {
    double big = 0.0;
    for (i = 0; i < n; i++) if (fabs(xe[i]) > big) big = fabs(xe[i]);
    if (!(maxamp > 1.e-5 * (1.0 + big))) {
      pic_fail("Periodic: the trajectory is not oscillating (it has settled to a steady state)", xe, fe, amp, scl, wk, wk2, &W);
      return 0;
    }
  }
  for (i = 0; i < n; i++) scl[i] = (amp[i] > 1.e-6 * maxamp) ? amp[i] : maxamp;

  /* going back from the end: the last time the trajectory crossed the section through the final point
     (normal to the flow there) in the forward direction, and was close to that point. A tight tolerance
     first, so that a nearby different loop of a period-2 orbit is not mistaken for the return. */
  ks = -1;
  {
    static const double tol[3] = {0.02, 0.06, 0.15};
    int it;
    for (it = 0; it < 3 && ks < 0; it++) {
      far = 0;
      for (k = N - 2; k > 0; k--) {
        double dd = 0.0, sec0 = 0.0, sec1 = 0.0, d1;
        for (i = 0; i < n; i++) {
          d = ((double)storage[i + 1][k] - xe[i]);
          d1 = ((double)storage[i + 1][k + 1] - xe[i]);
          if (fabs(d) / scl[i] > dd) dd = fabs(d) / scl[i];
          sec0 += d * fe[i]; sec1 += d1 * fe[i];
        }
        if (!far) { if (dd > 0.3) far = 1; continue; }
        if (sec0 < 0.0 && sec1 >= 0.0 && dd < tol[it]) { ks = k; break; }
      }
    }
  }
  if (ks < 0 || N - ks < 6) {
    pic_fail("Periodic: couldn't find a closed orbit at the end of the trajectory. Integrate longer so it settles onto a limit cycle",
             xe, fe, amp, scl, wk, wk2, &W);
    return 0;
  }
  {
    double *x0 = (double *)malloc(sizeof(double) * 4 * (size_t)n), *x1, *f0, *f1;
    double sec0 = 0.0, sec1 = 0.0, hh;
    if (!x0) { pic_fail("Periodic: out of memory", xe, fe, amp, scl, wk, wk2, &W); return 0; }
    x1 = x0 + n; f0 = x1 + n; f1 = f0 + n;
    for (i = 0; i < n; i++) {
      x0[i] = (double)storage[i + 1][ks]; x1[i] = (double)storage[i + 1][ks + 1];
      sec0 += (x0[i] - xe[i]) * fe[i]; sec1 += (x1[i] - xe[i]) * fe[i];
    }
    pic_f(&W, x0, f0); pic_f(&W, x1, f1);
    th = sec0 / (sec0 - sec1);
    hh = (double)storage[0][ks + 1] - (double)storage[0][ks];
    ts = (double)storage[0][ks] + th * hh;
    pic_hermite(n, th, hh, x0, f0, x1, f1, xs);
    free(x0);
  }
  T0e = (double)storage[0][N - 1] - ts;
  if (!(T0e > 0.0)) { pic_fail("Periodic: couldn't find the period", xe, fe, amp, scl, wk, wk2, &W); return 0; }

  /* ---- the return map to the section through xs, followed to its fixed point ---- */
  pic_f(&W, xs, fs);
  for (i = 0; i < n; i++) x[i] = xs[i];
  T = T0e; hg = T0e / 200.0; dT = 1.0; dx = 1.0;
  tstart = (double)clock() / CLOCKS_PER_SEC;
  for (cyc = 1; cyc <= PIC_MAXCYCLES; cyc++) {
    for (i = 0; i < n; i++) xp[i] = x[i];
    i = pic_cross(&W, x, 0.7 * T0e, 1.6 * T0e, xs, fs, scl, &Tc, &hg, xo);
    if (i != 0) {
      pic_fail("Periodic: couldn't follow the orbit around once. The equations may be too stiff for this, or the trajectory has not settled",
               xe, fe, amp, scl, wk, wk2, &W);
      return 0;
    }
    dx = 0.0;
    for (i = 0; i < n; i++) { d = fabs(x[i] - xp[i]) / scl[i]; if (d > dx) dx = d; }
    dT = fabs(Tc - T) / Tc;
    T = Tc;
    if (dx < 2.e-10 && dT < 2.e-10) { conv = 1; break; }
    if ((double)clock() / CLOCKS_PER_SEC - tstart > PIC_BUDGET) break;
  }
  if (!conv && (dx > 1.e-5 || dT > 1.e-5)) {
    pic_fail("Periodic: the orbit has not converged. It attracts too slowly, so integrate longer, or it is not a stable limit cycle",
             xe, fe, amp, scl, wk, wk2, &W);
    return 0;
  }
  /* x is now the fixed point on the section and T the period. But the return map can lock onto a multiple of
     the period (the orbit comes back to the same place after 2, 3, ... cycles), so check whether T/m also closes */
  {
    int m, again = 1;
    while (again) {
      again = 0;
      for (m = 8; m >= 2; m--) {
        double hh = T / m / 50.0, cl = 0.0;
        for (i = 0; i < n; i++) xo[i] = x[i];
        if (pic_to(&W, xo, T / m, &hh) != 0) continue;
        for (i = 0; i < n; i++) { d = fabs(xo[i] - x[i]) / scl[i]; if (d > cl) cl = d; }
        if (cl < 1.e-6) {
          /* the orbit has period T/m: continue the return map on it */
          T = T / m; again = 1;
          for (i = 0; i < n; i++) x[i] = xo[i];
          break;
        }
      }
    }
  }

  /* ---- the maximum of the x-axis variable, found along one period ---- */
  ns = PIC_NSAMPLE;
  vs = (double *)malloc(sizeof(double) * (size_t)(ns + 1));
  fv = (double *)malloc(sizeof(double) * (size_t)(ns + 1));
  brk = (int *)malloc(sizeof(int) * (size_t)(ns + 1));
  if (!vs || !fv || !brk) {
    free(vs); free(fv); free(brk);
    pic_fail("Periodic: out of memory", xe, fe, amp, scl, wk, wk2, &W);
    return 0;
  }
  {
    double dt = T / ns, h1 = dt;
    double *fx = wk2;
    /* sample the orbit, recording d(variable)/dt */
    for (i = 0; i < n; i++) xj[i] = x[i];
    for (k = 0; k <= ns; k++) {
      pic_f(&W, xj, fx);
      fv[k] = fx[iv];
      if (k < ns && pic_to(&W, xj, dt, &h1) != 0) {
        free(vs); free(fv); free(brk);
        pic_fail("Periodic: integration failed along the orbit", xe, fe, amp, scl, wk, wk2, &W);
        return 0;
      }
    }
    /* the maxima: where d(variable)/dt changes from + to - between two samples */
    nbr = 0;
    for (k = 0; k < ns; k++) if (fv[k] > 0.0 && fv[k + 1] <= 0.0) brk[nbr++] = k;
    if (nbr == 0) {
      free(vs); free(fv); free(brk);
      pic_fail("Periodic: the variable has no maximum along the orbit", xe, fe, amp, scl, wk, wk2, &W);
      return 0;
    }
    /* refine each maximum by bisection on d(variable)/dt = 0, and keep the highest */
    vbest = -1.e300; kbest = -1;
    for (ib = 0; ib < nbr; ib++) {
      k = brk[ib];
      for (i = 0; i < n; i++) xj[i] = x[i];
      h1 = dt;
      if (k > 0 && pic_to(&W, xj, k * dt, &h1) != 0) continue;
      ta = 0.0; tb = dt;
      for (i = 0; i < 70; i++) {
        tm = 0.5 * (ta + tb);
        for (j = 0; j < n; j++) xo[j] = xj[j];
        h1 = dt;
        if (pic_to(&W, xo, tm, &h1) != 0) break;
        pic_f(&W, xo, fx);
        g = fx[iv];
        if (g > 0.0) ta = tm; else tb = tm;
      }
      tm = 0.5 * (ta + tb);
      for (j = 0; j < n; j++) xo[j] = xj[j];
      h1 = dt;
      if (pic_to(&W, xo, tm, &h1) != 0) continue;
      vk = xo[iv];
      if (vk > vbest) { vbest = vk; kbest = k; for (j = 0; j < n; j++) bestx[j] = xo[j]; }
    }
    if (kbest < 0) {
      free(vs); free(fv); free(brk);
      pic_fail("Periodic: couldn't locate the maximum", xe, fe, amp, scl, wk, wk2, &W);
      return 0;
    }
  }
  free(vs); free(fv); free(brk);

  /* ---- how well does the orbit close? ---- */
  {
    double hh = T / 200.0;
    for (i = 0; i < n; i++) xo[i] = bestx[i];
    if (pic_to(&W, xo, T, &hh) != 0) closure = -1.0;
    else {
      closure = 0.0;
      for (i = 0; i < n; i++) { d = fabs(xo[i] - bestx[i]) / scl[i]; if (d > closure) closure = d; }
    }
  }

  /* ---- put the result where (G)o will find it ---- */
  for (i = 0; i < n; i++) last_ic[i] = bestx[i];
  TEND = T;
  if (TRANS != 0.0) TRANS = 0.0;

  plintf("Periodic orbit: period T = %.12g  (return map converged after %d cycles, orbit closes to %.2g)\n",
         T, cyc, closure);
  plintf("  start at the maximum of %s, %s = %.10g\n", uvar_names[iv], uvar_names[iv], bestx[iv]);
  plintf("  initial conditions:");
  for (i = 0; i < n && i < 12; i++) plintf(" %s=%.10g", uvar_names[i], bestx[i]);
  if (n > 12) plintf(" ...");
  plintf("\n  total set to %.12g, transient set to 0\n", T);
  {
    double step = fabs(DELTA_T) * (double)(NJMP > 0 ? NJMP : 1), q = T / step, r = fabs(q - floor(q + 0.5));
    if (r > 1.e-6 * (q > 1.0 ? q : 1.0))
      plintf("  note: T / (dt*njmp) = %.6g is not a whole number, so the stored orbit will stop %.4g short of one period\n",
             q, (q - floor(q)) * step);
  }
  sprintf(msg, "Period T=%.8g. Starting at max of %s; total set to T", T, uvar_names[iv]);
  if (script_is_active())
    plintf("%s\n", msg);
  /* a script gets a note, not an error box */
  else err_msg(msg);

  free(xe); free(fe); free(amp); free(scl); free(wk); free(wk2);
  pic_free(&W);
  return 1;
}
