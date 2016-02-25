#ifndef XPPAUT_XPPLIM_H
#define XPPAUT_XPPLIM_H

#include <stdio.h>

/* --- Macros --- */
#define MAXODE 5000
#define MAXODE1 4999
#define MAXDELAY 50
#define MAXPRIMEVAR (MAXODE-10)/2
#define MAXPAR 400
#define MAXFLAG 2000
#define MAX_SYMBS 10000
#define MAXUFUN 50
#define MAX_TAB 50
#define MAXKER 50
#define MAXNET 50
#define MAXMARK 200
#define MAXVEC 100
#define MAX_ANI_LINES 2000
#define MAX_INTERN_SET 500
#define MAX_LEN_SBOX 25
#define MaxIncludeFiles 10


/* --- Macros --- */
/*
The acutual max filename length is determined by the
FILENAME_MAX (see <stdio.h>), and usually 4096 -- but
this is huge and usually overkill.  On the otherhand
the old Xpp default string buffer size of 100 is a bit
restricitive for lengths of filenames. You could also
set this define in the Makefile or at compile time to
override the below definition.
*/

#if FILENAME_MAX < 300
# define XPP_MAX_NAME 300
#else
# define XPP_MAX_NAME FILENAME_MAX
#endif

#endif /* XPPAUT_XPPLIM_H */
