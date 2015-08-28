/* STARTUP PROCEDURE FOR UNIX FORTRAN PROGRAMS */

#include "stdio.h"
#include "signal.h"

#ifndef SIGIOT
#define SIGIOT SIGABRT
#endif

#ifdef NO__STDC
#define ONEXIT onexit
extern void f_exit();
#else
#ifdef __STDC__
#include "stdlib.h"
extern void f_exit(void);
#ifndef NO_ONEXIT
#define ONEXIT atexit
extern int atexit(void (*)(void));
#endif
#else
#ifndef NO_ONEXIT
#define ONEXIT onexit
extern void f_exit();
#endif
#endif
#endif

/*static void sigdie(s, kill)
register char *s;
int kill;
{*/
/* print error message, then clear buffers */
/*fflush(stderr);
fprintf(stderr, "%s\n", s);
f_exit();
fflush(stderr);

if(kill)
	{*/
	/* now get a core */
	/*signal(SIGIOT, 0);
	abort();
	}
else
	exit(1);
}*/

/*static void sigfdie(n)
int n;
{
sigdie("Floating Exception", 1);
}



static void sigidie(n)
int n;
{
sigdie("IOT Trap", 1);
}


static void sigqdie(n)
int n;
{
sigdie("Quit signal", 1);
}



static void sigindie(n)
int n;
{
sigdie("Interrupt", 0);
}



static void sigtdie(n)
int n;
{
sigdie("Killed", 0);
}
*/

int xargc;
char **xargv;

