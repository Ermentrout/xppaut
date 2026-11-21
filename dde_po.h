#ifndef _dde_po_h_
#define _dde_po_h_

/*  Periodic orbits of delay differential equations (constant delays).

    Orbit x(s), s in [0,1] (time rescaled by the period T), is represented
    by a continuous piecewise polynomial: ntst intervals, degree ncol on
    each, nodes equally spaced inside an interval, collocation at the Gauss
    points.  The delayed value x(t-tau) is x(s - tau/T) read off the same
    periodic polynomial, so tau > T works.  The Newton and the pseudo-arclength
    continuation use a dense Jacobian.

    The model is supplied through three callbacks, and the delayed values
    are fetched inside rhs() by calling ddepo_delay(tau,var), where var is
    the 1-based index of the state variable. */

typedef void   (*DDEPO_RHS)(void *ctx,double t,const double *x,double *f);
typedef void   (*DDEPO_SETPAR)(void *ctx,int ip,double v);
typedef double (*DDEPO_GETPAR)(void *ctx,int ip);

/* called for each accepted point; return nonzero to stop.
   status bits: 1 = parameter turning point (fold) between this and the
   previous point, 2 = last point, 4 = stopped because the amplitude collapsed
   (a Hopf point is near), 8 = the mesh was refined to reach this point,
   16 = a special point located between the previous and the next point (period doubling,
   torus or fold, see P->event); it is reported before the next regular point */
typedef int (*DDEPO_REPORT)(void *rctx,int istep,double par,double T,
			    const double *X,int status);

typedef struct {
  /* set by the caller */
  int n,ntst,ncol;
  DDEPO_RHS rhs; DDEPO_SETPAR setpar; DDEPO_GETPAR getpar; void *ctx;
  int ip;                 /* index of the continuation parameter          */
  double tol;             /* Newton tolerance (default 1e-8)              */
  int maxit;              /* Newton iterations (default 8)                */
  int verbose;
  double dtol;            /* mesh refinement: double NTST when the relative defect between the
                             collocation points exceeds dtol (default 3e-4; <0 = never refine) */
  int maxunk;             /* limit on ntst*ncol*n when refining (default 3000) */
  int nrefine;            /* number of refinements done (output) */
  int stop;               /* why ddepo_continue ended (output): 0 end of the parameter range,
                             1 nmax steps, 2 amplitude collapsed (Hopf nearby), 3 step size
                             fell below dsmin (Newton did not converge), 4 stopped by the callback */
  int locate;             /* 1: locate period doublings and tori along the branch (needs the Floquet multipliers) */
  double theta;           /* output of ddepo_tr_curve: the rotation angle at the current point */
  int event;              /* output, valid while the callback gets status bit 16: 1 period doubling, 2 torus (or neutral saddle), 3 fold */
  int mu_ok,nmu;          /* output: the nontrivial multipliers of the last tested point, by decreasing modulus */
  double mu_r[16],mu_i[16],mu_triv;
  int (*poll)(void *);    /* optional: called before every step (and retry); nonzero stops the run (stop code 5) */
  void *pollctx;
  int nunres;             /* points accepted although the mesh could not be refined further (output) */
  double maxdefect;       /* largest defect among those points (output) */
  double amin;            /* stop when the orbit amplitude falls below amin times the
                             largest amplitude seen (Hopf nearby); <0 disables; default 1e-3 */
  /* solution */
  double *X;              /* (ntst*ncol) nodes, n values each             */
  double T;               /* period                                       */
  /* private */
  int nmesh;
  double *Xd;             /* reference derivative, phase condition        */
  int hopf;               /* 1 after ddepo_hopf_guess: the next continue starts at a Hopf point */
  double *phi;            /* mode shape of the first-order Hopf orbit, at the mesh nodes */
  double *gz,*gw;         /* Gauss points and weights on [0,1]            */
  double *L,*Lp;          /* Lagrange values / derivatives at Gauss pts   */
} DDEPO;

int  ddepo_alloc(DDEPO *P,int n,int ntst,int ncol);
void ddepo_free(DDEPO *P);
/* fill X from a function orbit(ctx,s,x), s in [0,1), and set the period */
void ddepo_set_guess(DDEPO *P,double T,void (*orbit)(void *,double,double *),void *octx);
/* Start at a Hopf point of a steady state u0 (n values).  Delta(i omega) = i omega I - Cr + i S
   is the characteristic matrix there, with Cr and S in column-major layout (the matrices of
   auto_delay_hopf_matrices).  ddepo_hopf_vector returns the critical eigenvector va + i vb
   (largest component real, equal to 1) and the relative residual |Delta v|/|Delta|; a value
   much larger than 1e-3 means this is not a Hopf point.  ddepo_hopf_guess then sets the
   first-order orbit  u0 + eps Re(v e^{2 pi i s}),  T = 2 pi/omega, and the next call to
   ddepo_continue corrects it with the amplitude fixed (period and parameter free) and follows
   the branch with growing amplitude; the parameter direction is whatever the branch does. */
double ddepo_hopf_vector(int n,const double *Cr,const double *S,double omega,double *va,double *vb);
void ddepo_hopf_guess(DDEPO *P,const double *u0,double omega,const double *va,const double *vb,double eps);
/* Newton at the current parameter value; returns iterations used, -1 on failure */
int  ddepo_correct(DDEPO *P);
/* pseudo-arclength continuation in parameter P->ip; dir=+1/-1 initial direction.
   Returns the number of points, or <0: -1 memory/size, -2 Newton failed at the start,
   -3 tangent failed, -4 step size fell below dsmin at the first step, -5 start outside [pmin,pmax].
   The mesh is refined automatically (see dtol). */
int  ddepo_continue(DDEPO *P,double pmin,double pmax,int dir,double ds,
		    double dsmin,double dsmax,int nmax,DDEPO_REPORT cb,void *rctx);
/* relative defect of the equation between the collocation points (mesh error indicator) */
double ddepo_defect(DDEPO *P);
/* Floquet multipliers.  The monodromy operator is discretized by collocation time stepping of
   the variational equation along the orbit's own mesh over one period, with the history window
   of one longest delay; its eigenvalues approximate the leading multipliers.
   ddepo_multipliers returns, sorted by decreasing modulus, up to nmax of the raw eigenvalues
   (for a periodic orbit one of them is the trivial multiplier 1); returns their number, or -1.
   ddepo_floquet removes the trivial multiplier (the one closest to 1) and returns the nontrivial
   ones; *triv_err is |trivial multiplier - 1|, which indicates the discretization error. */
int ddepo_multipliers(DDEPO *P,double *mur,double *mui,int nmax);
int ddepo_floquet(DDEPO *P,double *mur,double *mui,int nmax,double *triv_err);
/* Period doubling branch switching.  ddepo_floquet_mode: the critical Floquet eigenfunction (multiplier
   closest to mu_target, normally -1) of the orbit in P, as nmesh*n node values over one period
   (anti-periodic when mu=-1), scaled to maximum 1; also returns the multiplier and a consistency
   measure.  ddepo_pd_guess: P2 (same n and ncol, twice the intervals) gets the first-order orbit of
   double period; the next ddepo_continue on P2 corrects it at fixed amplitude and follows the
   doubled branch (like a Hopf start).  Both return 0 on success. */
int ddepo_floquet_mode(DDEPO *P,double mu_target,double *v,double *mu,double *mismatch);
int ddepo_pd_guess(DDEPO *P2,const DDEPO *P1,const double *v,double eps);
/* Two-parameter continuation of a fold of cycles.  P holds the fold orbit (X,T) of a branch followed in
   parameter P->ip, with both parameters set to their values at the fold; ip2 is the second parameter.
   The fold is followed in (p1,p2) until p1 leaves [pmin,pmax] or p2 leaves [qmin,qmax] (the curve is ended
   on the boundary), nmax steps, or Newton fails; dir=+1 follows p1 upward at the start.  The callback
   gets every point (status 2 on the last); P->stop says why it ended (0 range, 1 nmax, 3 Newton failure,
   4 callback, 5 user).  Returns the number of points, or <0 on failure at the start (-1 size or memory,
   -2 Newton, -3 tangent, -4 null vector).  Unknowns: the orbit, T, p1, p2 and the null vector of the
   collocation Jacobian: about 2*(variables*NTST*NCOL) of them, with dense algebra. */
typedef int (*ddepo_cb2)(void *rctx,int istep,double p1,double p2,double T,const double *X,int status);
int ddepo_fold_curve(DDEPO *P,int ip2,double pmin,double pmax,double qmin,double qmax,int dir,
		     double ds,double dsmin,double dsmax,int nmax,ddepo_cb2 cb,void *rctx);
/* The same for a period doubling: P holds the orbit at the doubling (multiplier -1).  The extra unknown is the
   antiperiodic null vector of the collocation operator (about N unknowns, not N+1).  Same arguments and return
   values as ddepo_fold_curve (-4: no multiplier near -1 or no null vector). */
int ddepo_pd_curve(DDEPO *P,int ip2,double pmin,double pmax,double qmin,double qmax,int dir,
		   double ds,double dsmin,double dsmax,int nmax,ddepo_cb2 cb,void *rctx);
/* The same for a torus (Neimark-Sacker) point: P holds the orbit at the torus point, with a complex pair of
   multipliers exp(+-i theta) on the unit circle.  The extra unknowns are the angle theta and the complex null
   vector of the twisted operator (Y(s+1) = exp(i theta) Y(s)): about 2N+1 of them.  The angle of the current point
   is in P->theta in the callback (it is not reduced modulo 2 pi: it varies continuously along the curve).  A curve
   can end at a strong resonance (theta near 0 or pi, where the null space doubles); then Newton fails
   (P->stop = 3) with P->theta near a multiple of pi.  Same arguments and return values as ddepo_fold_curve. */
int ddepo_tr_curve(DDEPO *P,int ip2,double pmin,double pmax,double qmin,double qmax,int dir,
		   double ds,double dsmin,double dsmax,int nmax,ddepo_cb2 cb,void *rctx);
/* eigenvalues of a general real n x n matrix (row-major, destroyed); returns 0 on success */
int ddepo_eigenvalues(int n,double *a,double *wr,double *wi);
/* residual norm (max-norm) of the collocation equations at the current solution */
double ddepo_residual(DDEPO *P);
/* value of the orbit component c (0-based) at s in [0,1) */
double ddepo_eval(const DDEPO *P,int c,double s);

/* the hook: value of state var (1-based) at time t-tau on the orbit being solved */
double ddepo_delay(double tau,int var);

#endif
