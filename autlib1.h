#ifndef _autlib1_h_
#define _autlib1_h_


#include "wsfe.h"
#include "f2c.h"
#include "autevd.h"
#include "auto_x11.h"
#include "sfe.h"
/* autlib1.c */
int dfinit_(void);
int autoae_(doublereal *w, integer *iw, integer *itp, integer *ncpp, int (*funi)(void), int (*stpnt)(void));
int autobv_(doublereal *w, integer *iw, integer *itp, integer *ncpp, int (*funi)(void), int (*bcni)(void), int (*icni)(void), int (*stpnt)(void), int (*fnbpbv)(void));
int wae_(integer *itp, integer *lw, integer *liw);
int wbv_(integer *itp, integer *lw, integer *liw);
int wsae_(integer *itp, integer *ncpp, integer *lf, integer *ldfdu, integer *ldfdp, integer *laa, integer *lstud, integer *lstu, integer *lstrl, integer *lstrld, integer *lrhs, integer *ldu, integer *ludot, integer *lu, integer *luold, integer *lsmat, integer *lrnllv, integer *lu1, integer *lev, integer *lwkev, integer *lw, integer *lir, integer *lic, integer *liw, integer *m1aa, integer *m1stbf, integer *ndim2, integer *igenwts);
int cnstnt_(void);
int init1_(integer *itp, integer *igenwts);
int wsbv_(integer *itp, integer *ncpp, integer *lf, integer *ldfdu, integer *ldfdp, integer *laa, integer *lbb, integer *lcc, integer *ldd, integer *lups, integer *luldps, integer *lupldp, integer *ludtps, integer *lwbrbd, integer *lrhsa, integer *lrhsd, integer *ltint, integer *luintt, integer *ldups, integer *leqf, integer *luneq, integer *ltm, integer *ldtm, integer *ltm2, integer *lu, integer *lubc0, integer *lubc1, integer *ldbc, integer *luicd, integer *lficd, integer *ldicd, integer *lw, integer *litm, integer *lial, integer *lir, integer *lic, integer *liwbr, integer *liw, integer *m1aa, integer *m2aa, integer *m1bb, integer *m2bb, integer *m1cc, integer *m1dd, integer *m1u, integer *m1bc, integer *m1ic, integer *lp0, integer *lp1, integer *lpoin, integer *lev, integer *lwkev, integer *lsmat, integer *lrnllv, integer *igenwts);
/*int cnrlae_(int (*funi)(void), int (*stpnt)(void), integer *ibr, integer *m1aa, doublereal *aa, integer *m1stbf, doublereal *stud, doublereal *stu, doublereal *strl, doublereal *strld, doublereal *u, doublereal *rhs, doublereal *du, doublereal *udot, doublereal *uold, integer *ndim2, doublereal *smat, doublereal *rnllv, doublereal *f, integer *m1df, doublereal *dfdu, doublereal *dfdp, doublereal *u1, doublecomplex *ev, doublereal *wkev, integer *ir, integer *ic);
*/
int stpnus_(integer *ibr, doublereal *u, integer *ndm2, doublereal *smat, doublereal *dfdu, doublereal *dfuxx, doublereal *dfdp, doublereal *v, doublereal *f, integer *ir, integer *ic);
int stpnae_(integer *ibr, doublereal *u, integer *ndm2, doublereal *smat, doublereal *dfdu, doublereal *dfuxx, doublereal *dfdp, doublereal *v, doublereal *f, integer *ir, integer *ic);
int stprae_(int (*funi)(void), integer *istop, doublereal *rds, integer *nit, integer *ibr, integer *ntot, integer *m1aa, doublereal *aa, doublereal *rhs, doublereal *du, doublereal *udot, doublereal *uold, doublereal *u, doublereal *f, integer *m1df, doublereal *dfdu, doublereal *dfdp, integer *ir, integer *ic);
int contae_(doublereal *rds, doublereal *udot, doublereal *uold, doublereal *u);
int solvae_(int (*funi)(void), integer *istop, doublereal *rds, integer *nit, integer *ibr, integer *ntot, integer *m1aa, doublereal *aa, doublereal *rhs, doublereal *du, doublereal *uold, doublereal *u, doublereal *f, integer *m1df, doublereal *udot, doublereal *dfdu, doublereal *dfdp, integer *ir, integer *ic);
/*int lcspae_(doublereal (*fncs)(void), int (*funi)(void), integer *istop, integer *itp, doublereal *qcs, integer *ibr, integer *ntot, integer *m1aa, doublereal *aa, doublereal *rhs, doublereal *du, doublereal *udot, doublereal *uold, doublereal *u, doublereal *f, integer *m1df, doublereal *dfdu, doublereal *dfdp, doublecomplex *ev, doublereal *wkev, integer *ir, integer *ic);
*/
doublereal fnbpae_(logical *chng, int (*funi)(void), integer *m1aa, doublereal *aa, doublereal *u, doublereal *uold, doublereal *udot, doublereal *rhs, integer *m1df, doublereal *dfdu, doublereal *dfdp, doublecomplex *ev, doublereal *wkev, integer *ir, integer *ic, integer *ibr, integer *ntot);
doublereal fnlpae_(logical *chng, int (*funi)(void), integer *m1aa, doublereal *aa, doublereal *u, doublereal *uold, doublereal *udot, doublereal *rhs, integer *m1df, doublereal *dfdu, doublereal *dfdp, doublecomplex *ev, doublereal *wkev, integer *ir, integer *ic, integer *ibr, integer *ntot);
doublereal fnhbae_(logical *chng, int (*funi)(void), integer *m1aa, doublereal *aa, doublereal *u, doublereal *uold, doublereal *udot, doublereal *rhs, integer *m1df, doublereal *dfdu, doublereal *dfdp, doublecomplex *ev, doublereal *wkev, integer *ir, integer *ic, integer *ibr, integer *ntot);
doublereal fnuzae_(logical *chng, int (*funi)(void), integer *m1aa, doublereal *aa, doublereal *u, doublereal *uold, doublereal *udot, doublereal *rhs, integer *m1df, doublereal *dfdu, doublereal *dfdp, doublecomplex *ev, doublereal *wkev, integer *ir, integer *ic, integer *ibr, integer *ntot);
int stbif_(integer *nbif, integer *m1aa, doublereal *aa, integer *m1stbf, doublereal *stud, doublereal *stu, doublereal *strl, doublereal *strld, doublereal *du, doublereal *udot, doublereal *u, integer *m1df, doublereal *dfdu, doublereal *dfdp, integer *ir, integer *ic);
int swpnt_(integer *nbif, integer *ipos, doublereal *rds, integer *m1stbf, doublereal *stud, doublereal *stu, doublereal *strl, doublereal *strld, doublereal *udot, doublereal *u);
int swprc_(int (*funi)(void), integer *istop, integer *ibr, integer *ntot, integer *nit, integer *m1aa, doublereal *aa, doublereal *rhs, doublereal *du, doublereal *udot, doublereal *uold, doublereal *u, doublereal *u1, doublereal *f, integer *m1df, doublereal *dfdu, doublereal *dfdp, doublereal *rds, integer *ir, integer *ic);
int sthd_(void);
int headng_(integer *iunit, integer *n1, integer *n2);
void cnvrt_(char *ret_val, ftnlen ret_val_len, integer *i);
int stplae_(integer *istop, integer *itp, integer *nit, integer *ntot, integer *lab, integer *ibr, doublereal *u);
int wrline_(integer *ibr, integer *ntot, integer *itp, integer *lab, doublereal *vaxis, doublereal *u);
int wrtsp8_(integer *itp, integer *ntot, integer *lab, integer *ibr, doublereal *u);
int msh_(doublereal *tm);
int genwts_(integer *ncol);
int cpnts_(integer *ncol, doublereal *zm);
int cntdif_(integer *n, doublereal *d);
int wint_(integer *n, doublereal *wi);
int adptds_(doublereal *rds, integer *nit);
int adapt_(integer *nold, integer *ncold, integer *nnew, integer *ncnew, doublereal *tm, doublereal *dtm, integer *m1u, doublereal *ups, doublereal *vps, doublereal *tint, doublereal *uintt, doublereal *eqf, doublereal *uneq, doublereal *wrksp, doublereal *tm2, integer *itm, integer *ial);
int interp_(integer *ndim, integer *n, integer *nc, doublereal *tm, integer *m1u, doublereal *ups, integer *n1, integer *nc1, doublereal *tm1, doublereal *ups1, doublereal *tm2, integer *itm1);
int newmsh_(integer *m1u, doublereal *ups, integer *nold, integer *ncold, doublereal *tmold, doublereal *dtmold, integer *nnew, doublereal *tmnew, doublereal *eqf, doublereal *uneq, doublereal *wrksp, integer *ial, integer *iper);
int ordr_(integer *n, doublereal *tm, integer *n1, doublereal *tm1, integer *itm1);
int intwts_(integer *n, doublereal *z, doublereal *x, doublereal *wts);
int eqdf_(integer *ntst, integer *ndim, integer *ncol, doublereal *dtm, integer *m1u, doublereal *ups, doublereal *eqf, doublereal *wrksp, integer *iper);
int eig2_(integer *ndim, integer *m1a, doublereal *a, doublecomplex *ev, doublereal *wkev, integer *ier);
int eig_(integer *ndim, integer *m1a, doublereal *a, doublecomplex *ev, doublereal *wkev, integer *ier);
int nlvc_(integer *n, integer *m, integer *k, doublereal *a, doublereal *u, integer *ir, integer *ic);
int nrmlz_(integer *ndim, doublereal *v);
doublereal pi_(doublereal *r);
int ge_(integer *n, integer *m1a, doublereal *a, integer *nrhs, integer *m1u, doublereal *u, integer *m1f, doublereal *f, integer *ir, integer *ic);
int newlab_(integer *isw, integer *ibr, integer *lab);
int findl3_(integer *irs, integer *itp, integer *nfpar, logical *found);
int readl3_(integer *ips, integer *ibr, doublereal *u, doublereal *par);
int skip3_(integer *nskip, logical *eof3);
doublereal rinpr_(integer *ndim1, integer *m1u, doublereal *ups, doublereal *vps, doublereal *dtm);
doublereal rnrmsq_(integer *ndim1, integer *m1u, doublereal *ups, doublereal *dtm);
doublereal rintg_(integer *m1u, integer *ic, doublereal *ups, doublereal *dtm);
doublereal rnrm2_(integer *m1u, integer *ic, doublereal *ups, doublereal *dtm);
doublereal rmxups_(integer *m1u, integer *i, doublereal *ups);
doublereal rmnups_(integer *m1u, integer *i, doublereal *ups);
int scalebb_(integer *m1u, doublereal *dvps, doublereal *rld, doublereal *dtm);




#endif
