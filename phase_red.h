#ifndef _phase_red_h_
#define _phase_red_h_
/* Phase reduction of a stable limit cycle (see phase_red.c). Stage 1: the adjoint Z. */
#define PHASE_MAXN 60
typedef struct {
  int n, P, Ng;                 /* equations, points on one period (P-1 intervals), grid steps */
  int multiple;                 /* m > 1 if the period given was m times the true period (the adjoint is for the true one) */
  double T;                     /* the period (changed slightly if the grid orbit had to be closed) */
  double closure0, closure1;    /* grid-orbit closure error before and after the Newton polish (-1: not done) */
  double errNorm, errPer;       /* max |Z.f - 1| over the orbit, and the periodicity error of Z */
  double *th;                   /* phase 0..1, P values */
  double *U, *F, *Z;            /* orbit, vector field and adjoint, P rows of n */
  /* stage 2: the isostable quantities, for the slowest nontrivial Floquet multiplier nu = exp(kappa T) when it is real, positive and unique */
  int iso, isoSkip;             /* iso = 1 if computed; otherwise isoSkip is a reason (phase_iso_note) */
  double nu, kappa, nuRatio;    /* nuRatio = |next multiplier| / |nu| (in the sense of the Floquet exponents) */
  double errDot, errZR, errIU, errPerR, errPerI;   /* max |I.rho-1|, max |Z.rho|, I.f, periodicity of rho and of I */
  double *I, *R;                /* the isostable response I and the Floquet eigenfunction rho, P rows of n */
  int hasB, bNote;              /* B computed, or why not */
  double errPerB;
  double *B;                    /* the first correction to the phase response, P rows of n (up to PHASE_BMAXN variables) */
} PHASE;
#define PHASE_BMAXN 24
const char *phase_iso_note(int code);   /* the reason for isoSkip */
const char *phase_b_note(int code);     /* the reason for bNote */
int phase_adjoint(int n, const double *x0, double T, int N0, PHASE *ph);
void phase_free(PHASE *ph);

/* Stage 3: the interaction functions of two identical cells coupled through G(X, X'), X the cell and X' its partner.
   G is supplied as a function that returns the coupling in every equation (n values) for the state of the cell and of the partner.
   With phases theta_i = t + phi_i, the phase difference phi = phi_2 - phi_1 obeys  phi' = eps [g1 + eps (g2 + g3)](phi):
      H(phi)  = < Z(t).G(U(t), U(t+phi)) >,                 g1(phi) = H(-phi) - H(phi)
      H2(phi) = A(phi) + Bf(phi),  A = < f_phi h2 >,  Bf = < f_(-phi)(t+phi) h3 >,   g2 = A(-phi) - A(phi), g3 = Bf(-phi) - Bf(phi)
   (Wilson and Ermentrout, Phys. Rev. Lett. 123, 164101 (2019); see phase_red.c). Arrays have PHASE_HM+1 values at phi = j T/PHASE_HM. */
#define PHASE_HM 100
typedef void (*PH_COUPLE)(const double *own, const double *partner, double *out, void *ctx);
typedef struct {
  int Mp; double T;
  double *H, *H2, *A, *Bf, *g1, *g2, *g3, *g23;
  double gamma, gammaA, gammaB;   /* gamma = (g2'(0) + g3'(0))/4: the coefficient of the fourth derivative for diffusive coupling */
  int notes;                      /* bit 1: kappa*T small (expansion breaks down), bit 2: next multiplier not much smaller */
} PHASEH;
/* idx: the nc equations that receive coupling (the others must give 0). Needs I, rho and B in ph, and (P-1) a multiple of PHASE_HM. */
int phase_hfun(const PHASE *ph, const int *idx, int nc, PH_COUPLE G, void *ctx, PHASEH *out);
void phase_h_free(PHASEH *h);
#endif
