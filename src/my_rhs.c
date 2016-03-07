#include "my_rhs.h"

#include <stdlib.h>

#include "dae_fun.h"
#include "extra.h"
#include "form_ode.h"
#include "main.h"
#include "markov.h"
#include "parserslow.h"
#include "simplenet.h"
#include "xpplim.h"


/* --- Functions --- */
void extra(double *y__y, double t, int nod, int neq)  {
	int i;
	if(nod>=neq) {
		return;
	}
	SETVAR(0,t);
	for(i=0;i<nod;i++) {
		SETVAR(i+1,y__y[i]);
	}
	for(i=nod+FIX_VAR;i<nod+FIX_VAR+NMarkov;i++) {
		SETVAR(i+1,y__y[i-FIX_VAR]);
	}
	for(i=nod;i<nod+FIX_VAR;i++) {
		SETVAR(i+1,evaluate(my_ode[i]));
	}
	/* I dont think this is generally needed  */
	/* do_in_out();   */
	for(i=nod+NMarkov;i<neq;i++) {
		y__y[i]=evaluate(my_ode[i+FIX_VAR-NMarkov]);
	}
}


void fix_only(void) {
	int i;
	for(i=NODE;i<NODE+FIX_VAR;i++) {
		SETVAR(i+1,evaluate(my_ode[i]));
	}
}


int main(int argc, char **argv) {
	do_main(argc,argv);
	exit(0);
}


int my_rhs(double t, double *y, double *ydot, int neq) {
	int i;
	SETVAR(0,t);
	for(i=0;i<NODE;i++) {
		SETVAR(i+1,y[i]);
	}
	for(i=NODE;i<NODE+FIX_VAR;i++) {
		SETVAR(i+1,evaluate(my_ode[i]));
	}
	eval_all_nets();
	do_daes();
	do_in_out();
	for(i=0;i<NODE;i++) {
		ydot[i]=evaluate(my_ode[i]);
	}
	return(1);
}


void rhs_only(double *y,double *ydot) {
	int i;
	for(i=0;i<NODE;i++) {
		ydot[i]=evaluate(my_ode[i]);
	}
}


void set_fix_rhs(double t, double *y) {
	int i;
	SETVAR(0,t);
	for(i=0;i<NODE;i++) {
		SETVAR(i+1,y[i]);
	}
	for(i=0;i<NMarkov;i++) {
		SETVAR(i+1+NODE+FIX_VAR,y[i+NODE]);
	}
	for(i=NODE;i<NODE+FIX_VAR;i++) {
		SETVAR(i+1,evaluate(my_ode[i]));
	}
	eval_all_nets();
	do_in_out();
}


void update_based_on_current(void) {
	int i;
	for(i=NODE;i<NODE+FIX_VAR;i++) {
		SETVAR(i+1,evaluate(my_ode[i]));
	}
	eval_all_nets();
	do_in_out();
}
