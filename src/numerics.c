/* The input is primitive and eventually, I want to make it so
 * that it uses nice windows for input.
 * For now, I just will let it remain command driven
 */
#include "numerics.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "adj2.h"
#include "browse.h"
#include "color.h"
#include "delay_handle.h"
#include "flags.h"
#include "form_ode.h"
#include "gear.h"
#include "ggets.h"
#include "graf_par.h"
#include "integrate.h"
#include "load_eqn.h"
#include "main.h"
#include "many_pops.h"
#include "markov.h"
#include "menu.h"
#include "menudrive.h"
#include "menus.h"
#include "odesol2.h"
#include "parserslow.h"
#include "pop_list.h"
#include "pp_shoot.h"
#include "storage.h"
#include "struct.h"
#include "volterra2.h"


/*   This is numerics.c
 *   The input is primitive and eventually, I want to make it so
	that it uses nice windows for input.
	For now, I just will let it remain command driven
*/

/* --- Forward declarations --- */
static void check_pos(int *j);
static void get_method(void);
static void ruelle(void);
int (*solver)();


/* --- Data --- */
int METHOD;
int cv_bandflag=0,cv_bandupper=1,cv_bandlower=1;


/* --- Functions --- */
void check_delay(void) {
	if(DELAY>0.0) {
		free_delay();
		if(alloc_delay(DELAY)) {
			INFLAG=0; /*  Make sure no last ics allowed */
		}
	} else {
		free_delay();
	}
}


void do_meth(void) {
	if(NKernel>0) {
		METHOD=VOLTERRA;
	}
	switch(METHOD) {
	case 0:
		solver=discrete;
		DELTA_T=1;
		break;
	case 1:
		solver=euler;
		break;
	case 2:
		solver=mod_euler;
		break;
	case 3:
		solver=rung_kut;
		break;
	case 4:
		solver=adams;
		break;
	case 5:
		NJMP=1;
		break;
	case 6:
		solver=volterra;
		break;
	case SYMPLECT:
		solver=symplect3;
		break;
	case BACKEUL:
		solver=bak_euler;
		break;
	case RKQS:
	case STIFF:
	case CVODE:
	case DP5:
	case DP83:
	case RB23:
		NJMP=1; break;
	default: solver=rung_kut;
	}
}


void get_num_par(int ch) {
	double temp;
	int tmp;
	switch(ch) {
	case 'a':
		make_adj();
		break;
	case 't':
		/* total */
		new_float("total :",&TEND);
		FOREVER=0;
		if(TEND<0) {
			FOREVER=1;
			TEND=-TEND;
		}
		break;
	case 's':
		/* start */
		new_float("start time :",&T0);
		break;
	case 'r':
		/* transient */
		new_float("transient :",&TRANS);
		break;
	case 'd':
		/* DT */
		temp=DELTA_T;
		new_float("Delta t :",&DELTA_T);
		if(DELTA_T==0.0) {
			DELTA_T=temp;
		}
		if(DELAY>0.0) {
			free_delay();
			if(alloc_delay(DELAY)) {
				INFLAG=0; /*  Make sure no last ics allowed */
			}
		} else {
			free_delay();
		}
		if(NKernel>0) {
			INFLAG=0;
			MyStart=1;
			alloc_kernels(1);
		}

		break;
	case 'n':
		/* ncline */
		new_int("ncline mesh :",&NMESH);
		check_pos(&NMESH);
		break;
	case 'v':
		new_int("Maximum iterates :",&BVP_MAXIT);
		check_pos(&BVP_MAXIT);
		new_float("Tolerance :",&BVP_TOL);
		new_float("Epsilon :",&BVP_EPS);
		reset_bvp();
		break;
	case 'i':
		/* sing pt */
		new_int("Maximum iterates :",&EVEC_ITER);
		check_pos(&EVEC_ITER);
		new_float("Newton tolerance :",&EVEC_ERR);
		new_float("Jacobian epsilon :",&NEWT_ERR);
		if(NFlags>0) {
			new_float("SMIN :",&STOL);
		}
		break;
	case 'o':
		/* noutput */
		new_int("n_out :",&NJMP);
		check_pos(&NJMP);
		break;
	case 'b':
		/* bounds */
		new_float("Bounds :",&BOUND);
		BOUND=fabs(BOUND);
		break;
	case 'm':
		/* method */
		get_method();
		if(METHOD==VOLTERRA && NKernel==0) {
			err_msg("Volterra only for integral eqns");
			METHOD=4;
		}
		if(NKernel>0) {
			METHOD=VOLTERRA;
		}
		if(METHOD==GEAR || METHOD==RKQS || METHOD==STIFF) {
			new_float("Tolerance :",&TOLER);
			new_float("minimum step :",&HMIN);
			new_float("maximum step :",&HMAX);
		}
		if(METHOD==CVODE || METHOD==DP5 || METHOD==DP83 || METHOD==RB23)  {
			new_float("Relative tol:",&TOLER);
			new_float("Abs. Toler:",&ATOLER);
		}

		if(METHOD==BACKEUL || METHOD==VOLTERRA) {
			new_float("Tolerance :",&EulTol);
			new_int("MaxIter :",&MaxEulIter);
		}
		if(METHOD==VOLTERRA) {
			tmp=MaxPoints;
			new_int("MaxPoints:",&tmp);
			new_int("AutoEval(1=yes) :",&AutoEvaluate);
			allocate_volterra(tmp,1);
		}
		if(METHOD==CVODE || METHOD==RB23) {
			new_int("Banded system(0/1)?",&cv_bandflag);
			if(cv_bandflag==1) {
				new_int("Lower band:",&cv_bandlower);
				new_int("Upper band:",&cv_bandupper);
			}
		}
		if(METHOD==SYMPLECT) {
			if((NODE%2)!=0) {
				err_msg("Symplectic is only for even dimensions");
				METHOD=4;
			}
		}
		break;
	case 'e':
		/* delay */
		if(NDELAYS==0) {
			break;
		}
		new_float("Maximal delay :",&DELAY);
		new_float("real guess :", &AlphaMax);
		new_float("imag guess :", &OmegaMax);
		new_int("DelayGrid :",&DelayGrid);
		if(DELAY>0.0) {
			free_delay();
			if(alloc_delay(DELAY)) {
				INFLAG=0; /*  Make sure no last ics allowed */
			}
		} else {
			free_delay();
		}
		break;
	case 'c':
		/* color */
		if(COLOR==0) {
			break;
		}
		set_col_par();
		break;
	case 'h':
		do_stochast();
		break;
	case 'f':
		/* FFT */
		break;
	case 'p':
		/*Poincare map */
		get_pmap_pars();
		break;
	case 'u':
		/* ruelle */
		ruelle();
		break;
	case 'k':
		/*lookup table */
		new_lookup();
		break;
	case 27:
		do_meth();
		TEND=fabs(TEND);
		alloc_meth();
		help();
		break;
	}  /* End num switch */
}


void set_delay(void) {
	if(NDELAYS==0) {
		return;
	}
	if(DELAY>0.0) {
		free_delay();
		if(alloc_delay(DELAY)) {
			INFLAG=0;
		}
	}
}


void set_total(double total) {
	int n;
	n=(total/fabs(DELTA_T))+1;
	TEND=n*fabs(DELTA_T);
}


/* --- Static functions --- */
static void check_pos(int *j) {
	if(*j<=0) {
		*j=1;
	}
}


static void get_method(void) {
	char ch;
	int i;
	int nmeth;

	Window temp=main_win;
	static char *n[]={"(D)iscrete","(E)uler","(M)od. Euler",
					  "(R)unge-Kutta","(A)dams","(G)ear","(V)olterra","(B)ackEul",
					  "(Q)ualst.RK4","(S)tiff","(C)Vode","DoPri(5)","DoPri(8)3",
					  "Rosen(2)3","sYmplectic"};
	static char key[]="demragvbqsc582y";

#ifdef CVODE_YES
	nmeth=15;
#else
	nmeth=15;
#endif
	ch = (char)pop_up_list(&temp,"Method",n,key,nmeth,15,METHOD,10,DCURY+8,
						   meth_hint,info_pop,info_message);
	for(i=0;i<nmeth;i++)
		if(ch==key[i]) {
			METHOD=i;
		}
	if(i>(nmeth-1)) {
		i=nmeth-1;
	}
}


static void ruelle(void) {
	new_int("x-axis shift ",&(MyGraph->xshft));
	new_int("y-axis shift ",&(MyGraph->yshft));
	new_int("z-axis shift",&(MyGraph->zshft));
	if(MyGraph->xshft<0) {
		MyGraph->xshft=0;
	}
	if(MyGraph->yshft<0) {
		MyGraph->yshft=0;
	}
	if(MyGraph->zshft<0) {
		MyGraph->zshft=0;
	}
}
