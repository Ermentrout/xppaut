/* Derived parameter stuff !!  */
#include "derived.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "calc.h"
#include "form_ode.h"
#include "ggets.h"
#include "parserslow.h"

/* --- Macros --- */
#define MAXDERIVED 200

/* --- Types --- */
typedef struct {
	int index,*form;
	char *rhs;
	double value;
} DERIVED;

/* --- Data --- */
static DERIVED derived[MAXDERIVED];
static int nderived=0;

/* --- Functions --- */
/* this adds a derived quantity  */
int add_derived(char *name, char *rhs) {
	int n=strlen(rhs)+2;
	int i0;
	if(nderived>=MAXDERIVED) {
		plintf(" Too many derived constants! \n");
		return(1);
	}
	i0=nderived;
	derived[i0].rhs=(char *)malloc(n);
	/* save the right hand side */
	strcpy(derived[i0].rhs,rhs);
	/* this is the constant to which it addresses */
	derived[i0].index=NCON;
	/* add the name to the recognized symbols */
	plintf(" derived constant[%d] is %s = %s\n",NCON,name,rhs);
	nderived++;
	return(add_con(name,0.0));
}


/* This compiles all of the formulae. t is called only once during the session */
int compile_derived(void) {
	int i,k;
	int f[DEFAULT_STRING_LENGTH],n;
	for(i=0;i<nderived;i++) {
		if(add_expr(derived[i].rhs,f,&n)==1) {
			plintf(" Bad right-hand side for derived parameters \n");
			return(1);
		}
		derived[i].form=(int *)malloc(sizeof(int)*(n+2));
		for(k=0;k<n;k++) {
			derived[i].form[k]=f[k];
		}
	}
	evaluate_derived();
	return 0;
}


/* This evaluates all derived quantities in order of definition called before
 * any integration or numerical computation and after changing parameters
 * and constants
 */
void evaluate_derived(void) {
	int i;
	for(i=0;i<nderived;i++) {
		derived[i].value=evaluate(derived[i].form);
		constants[derived[i].index]=derived[i].value;
	}
}
