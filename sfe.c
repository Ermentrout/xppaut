/* sequential formatted external common routines*/
#include "err.h"

#include "sfe.h"

#include "open.h"
#include "fmt.h"


extern char *fmtbuf;

integer e_rsfe()
{	int n;
	n=en_fio();
	if (cf == stdout)
		fflush(stdout);
	else if (cf == stderr)
		fflush(stderr);
	fmtbuf=NULL;
	return(n);
}
integer c_sfe(a) cilist *a; /* check */
{	unit *p;
	if(a->ciunit >= MXUNIT || a->ciunit<0)
		err(a->cierr,101,"startio");
	p = &units[a->ciunit];
	if(p->ufd==NULL && fk_open(SEQ,FMT,a->ciunit)) err(a->cierr,114,"sfe")
	if(!p->ufmt) err(a->cierr,102,"sfe")
	return(0);
}
integer e_wsfe()
{	return(e_rsfe());
}
