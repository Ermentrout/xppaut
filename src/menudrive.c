/* the menu for XPP commands
 * this calls any command it also has lots of the direct X Gui stuff
 */
#include "menudrive.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>
#include <X11/Xlib.h>

#include "adj2.h"
#include "auto_nox.h"
#include "browse.h"
#include "calc.h"
#include "edit_rhs.h"
#include "extra.h"
#include "ggets.h"
#include "graf_par.h"
#include "histogram.h"
#include "init_conds.h"
#include "integrate.h"
#include "kinescope.h"
#include "load_eqn.h"
#include "lunch-new.h"
#include "many_pops.h"
#include "main.h"
#include "markov.h"
#include "menudrive.h"
#include "menus.h"
#include "mykeydef.h"
#include "nullcline.h"
#include "numerics.h"
#include "parserslow.h"
#include "pop_list.h"
#include "pp_shoot.h"
#include "tabular.h"
#include "torus.h"
#include "tutor.h"
#include "txtread.h"
#include "util/timeutil.h"

#include "bitmap/alert.bitmap"


/* --- Types --- */
typedef struct {
	Window w;
	char text[DEFAULT_STRING_LENGTH];
	int here;
} MSGBOXSTRUCT;


/* --- Forward declarations --- */
static void do_file_com(int com);
static void get_pmap_pars_com(int l);
static void quick_num(int com);
static void set_col_par_com(int i);

/* --- Data --- */
static int status;
static MSGBOXSTRUCT MsgBox;

int screen;

/* --- Functions --- */
void do_tutorial(void) {
	printf("Running tutorial!\n");
	int tut=0;
	while (1)	{
		char ans = (char)two_choice("Next","Done",tutorial[tut],"nd",DisplayWidth/2,DisplayHeight/2,
									RootWindow(display,screen),"Did you know you can...");
		if(ans=='d') {
			break;
		}
		tut++;
		tut = tut % N_TUTORIAL;
	}
}


void edit_xpprc(void) {
	pid_t child_pid;
	char rc[DEFAULT_STRING_LENGTH];
	char editor[DEFAULT_STRING_LENGTH];
	int child_status;

	char* ed = getenv("XPPEDITOR");
	if ((ed == NULL) || (strlen(ed)==0)) {
		err_msg("Environment variable XPPEDITOR needs to be set.");
		return;
	} else {
		sprintf(editor, "%s",ed);
	}

	child_pid = fork();
	if (child_pid == 0)	{
		sprintf(rc,"%s/.xpprc",getenv("HOME"));
		char *const args[] = {editor,rc,NULL};
		execvp(editor,args);
		wait(&child_status);
		return;
	} else {
		if (child_pid == -1) {
			err_msg("Unable to fork process for editor.");
		}
		return;
	}
}


void xpp_hlp(void) {
	char cmd[DEFAULT_STRING_LENGTH];

	if(getenv("XPPHELP")==NULL) {
		err_msg("Environment variable XPPHELP undefined.");
		return;
	}
	if (getenv("XPPBROWSER")==NULL) {
		err_msg("Environment variable XPPBROWSER undefined.");
		return;
	}
	sprintf(cmd,"file:///%s",getenv("XPPHELP"));
	if(fork()==0) {
		execlp(getenv("XPPBROWSER"),getenv("XPPHELP"),cmd,(char *)0);
		perror("Unable to open browser. Check your XPPBROWSER and XPPHELP environement variables.");
		exit(1);
	} else {
		wait(&status);
	}
}


void MessageBox(char *m) {
	int wid=strlen(m)*DCURX+20;
	int hgt=4*DCURY;
	MsgBox.w=make_plain_window(RootWindow(display,screen),
							   DisplayWidth/2,DisplayHeight/2, wid,hgt,4);

	make_icon((char*)alert_bits,alert_width,alert_height,MsgBox.w);
	MsgBox.here=1;
	set_window_title(MsgBox.w,"Yo!");
	strcpy(MsgBox.text,m);
	ping();
}


void RedrawMessageBox(Window w) {
	if(w==MsgBox.w) {
		Ftext(10,2*DCURY,MsgBox.text,MsgBox.w);
	}
}


void KillMessageBox(void) {
	if(MsgBox.here==0) {
		return;
	}
	MsgBox.here=0;
	waitasec(ClickTime);
	XDestroyWindow(display,MsgBox.w);
}


int TwoChoice(char *c1,char *c2, char *q,char *key) {
	return two_choice(c1,c2,q,key,DisplayWidth/2,DisplayHeight/2,
					  RootWindow(display,screen),NULL);
}


int GetMouseXY(int *x,int *y) {
	return get_mouse_xy(x,y,draw_win);
}


void FlushDisplay(void) {
	XFlush(display);
}


void clear_draw_window(void) {
	clr_scrn();
	hi_lite(draw_win);
}

void drw_all_scrns(void) {
	int i;
	int ic=current_pop;
	if(SimulPlotFlag==0) {
		redraw_all();
		return;
	}
	for(i=0;i<num_pops;i++) {
		make_active(ActiveWinList[i],1);
		redraw_all();
	}
	make_active(ic,1);
	hi_lite(draw_win);
}

void clr_all_scrns(void) {
	int i;
	int ic=current_pop;
	if(SimulPlotFlag==0) {
		clr_scrn();
		hi_lite(draw_win);
		return;
	}
	for(i=0;i<num_pops;i++) {
		make_active(ActiveWinList[i],1);
		clr_scrn();
	}
	make_active(ic,1);
	hi_lite(draw_win);
}


void run_the_commands(int com) {
	if(com<0) {
		return;
	}
	if(com<=MAX_M_I) {
		do_init_data(com);
		return;
	}
	if(com==M_C) {
		cont_integ();
		return;
	}
	if(com>=M_SG&&com<=M_SC) {
		find_equilib_com(com-M_SG);
		return;
	}
	if(com>=M_NFF&&com<=M_NFA) {
		froz_cline_stuff_com(com-M_NFF);
		return;
	}
	if(com>=M_NN&&com<=M_NS) {
		new_clines_com(com-M_NN);
		return;
	}
	if(com>=M_DD&&com<=M_DS) {
		direct_field_com(com-M_DD);
		if((com-M_DD)==1) {
			return;
		}
		create_new_cline();
		redraw_the_graph();
		return;
	}
	if(com>=M_WW&&com<=M_WS) {
		window_zoom_com(com-M_WW);
		return;
	}
	if(com>=M_AA&&com<=M_AC) {
		do_torus_com(com-M_AA);
		return;
	}
	if(com>=M_KC&&com<=M_KM) {
		do_movie_com(com-M_KC);
		return;
	}
	if(com>=M_GA &&com<=M_GC) {
		add_a_curve_com(com-M_GA);
		return;
	}
	if(com>=M_GFF&&com<=M_GFO) {
		freeze_com(com-M_GFF);
		return;
	}
	if(com>=M_GCN&&com<=M_GCU) {
		change_cmap_com(com-M_GCN);
		redraw_dfield();
		return;
	}
	if(com==M_GFKK||com==M_GFKN) {
		key_frz_com(com-M_GFKN);
		return;
	}
	if(com==M_UKE||com==M_UKV) {
		new_lookup_com(com-M_UKE);
		return;
	}
	if(com==M_R) {
		drw_all_scrns();
		return;
	}
	if(com==M_EE) {
		clr_all_scrns();
		DF_FLAG=0;
		return;
	}
	if(com==M_X) {
		xi_vs_t();
		return;
	}
	if(com==M_3) {
		get_3d_par_com();
		return;
	}
	if(com==M_P) {
		new_parameter();
		return;
	}
	if(com>=M_MC&&com<=M_MS) {
		do_windows_com(com-M_MC);
		return;
	}
	/* CLONE */
	if(com>=M_FP&&com<=M_FL) {
		do_file_com(com);
		return;
	}
	if(com>=M_TT&&com<=M_TS) {
		do_gr_objs_com(com-M_TT);
		return;
	}
	if(com>=M_TEM&&com<=M_TED) {
		edit_object_com(com-M_TEM);
		return;
	}
	if(com>=M_BR &&com<=M_BH) {
		find_bvp_com(com-M_BR);
		return;
	}
	if(com>=M_V2&&com<=M_VT) {
		change_view_com(com-M_V2);
	}
	if(com>=M_UAN&&com<=M_UAR) {
		make_adj_com(com-M_UAN);
	}
	if(com>=M_UCN&&com<=M_UCA) {
		set_col_par_com(com-M_UCN);
	}
	if(com>=M_UPN&&com<=M_UPP) {
		get_pmap_pars_com(com-M_UPN);
	}
	if(com>=M_UHN&&com<=M_UH2) {
		do_stochast_com(com-M_UHN);
	}
	if(com>=M_UT && com<=M_UC) {
		quick_num(com-M_UT);
	}
}


void do_stochast(void) {
	static char *n[]={"New seed","Compute","Data","Mean","Variance","Histogram",
					  "Old hist","Fourier","Power","fIt data","Stat","Liapunov",
					  "stAutocor","Xcorrel etc","spEc.dns","2D-hist"};
	static char key[]="ncdmvhofpislaxe2";
	Window temp=main_win;
	char ch;
	int i;
	ch=(char)pop_up_list(&temp,"Stochastic",n,key,16,10,0,10,2*DCURY+8,
						 stoch_hint,info_pop,info_message);
	for(i=0;i<16;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<16) {
		run_the_commands(M_UHN+i);
	}
}


void get_pmap_pars(void) {
	static char *map[]={"(N)one","(S)ection","(M)ax/min","(P)eriod"};
	static char mkey[]="nsmp";
	char ch;
	Window temp=main_win;
	int i;

	ch=(char)pop_up_list(&temp,"Poincare map",map,mkey,4,13,POIMAP,10,6*DCURY+8,
						 map_hint,info_pop,info_message);

	for(i=0;i<4;i++) {
		if(ch==mkey[i]) {
			break;
		}
	}
	if(i>=0&&i<4) {
		run_the_commands(M_UPN+i);
	}
}


void set_col_par(void) {
	char ch;
	Window tempw=main_win;
	static char *n[]={"(N)o color","(V)elocity","(A)nother quantity"};
	static char key[]="nva";
	int i;
	ch=(char)pop_up_list(&tempw,"Color code",n,key,
						 3,11,0,10,12*DCURY+8,color_hint,
						 info_pop,info_message);
	for(i=0;i<3;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<3) {
		run_the_commands(i+M_UCN);
	}
}


void make_adj(void) {
	Window temp=main_win;
	static char *n[]={"(N)ew adj","(M)ake H","(A)djoint","(O)rbit","(H)fun","(P)arameters","(R)ange"};
	static char key[]="nmaohpr";
	char ch;
	int i;
	ch=(char)pop_up_list(&temp,"Adjoint",n,key,7,10,0,10,11*DCURY+8,adj_hint,
						 info_pop,info_message);
	for(i=0;i<7;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<7) {
		run_the_commands(M_UAN+i);
	}
}


void do_gr_objs(void) {
	char ch;
	int i;
	static char *list[]={"(T)ext","(A)rrow","(P)ointer","(M)arker",
						 "(E)dit","(D)elete all","marker(S)"};
	static char key[]="tapmeds";
	static char title[]="Text,etc";
	static char *elist[]={"(M)ove","(C)hange","(D)elete"};
	static char ekey[]="mcd";
	static char etitle[]="Edit";
	Window temp=main_win;
	ch=(char)pop_up_list(&temp,title,list,key,7,10,0,10,10*DCURY+8,
						 text_hint,info_pop,info_message);
	if(ch==ESC) {
		return;
	}
	if(ch=='e') {
		ch=(char)pop_up_list(&temp,etitle,elist,ekey,3,9,0,10,10*DCURY+8,

							 edit_hint,info_pop,info_message);
		if(ch==ESC) {
			return;
		}
		for(i=0;i<3;i++) {
			if(ch==ekey[i]) {
				break;
			}
		}
		if(i>=0&&i<3) {
			run_the_commands(M_TEM+i);
		}
		return;
	}
	for(i=0;i<7;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<7) {
		run_the_commands(M_TT+i);
	}
}


void new_lookup(void) {
	static char *n[]={"(E)dit","(V)iew"};
	static char key[]="ev";
	char ch;
	Window temp=main_win;
	if(NTable==0) {
		return;
	}
	ch=(char)pop_up_list(&temp,"Tables",n,key,2,12,1,10,11*DCURY+8,
						 tab_hint,info_pop,info_message);
	if(ch==key[0]) {
		run_the_commands(M_UKE);
	}
	if(ch==key[1]) {
		run_the_commands(M_UKV);
	}
}

void find_bvp(void) {
	static char *n[]={"(R)ange","(N)o show","(S)how","(P)eriodic"};
	static char key[]="rnsp";
	char ch;
	int i;
	Window temp=main_win;
	ch=(char)pop_up_list(&temp,"Bndry Value Prob",n,key,4,16,1,10,6*DCURY+8,
						 bvp_hint,info_pop,info_message);
	if(ch==ESC) {
		return;
	}
	for(i=0;i<4;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<4) {
		run_the_commands(M_BR+i);
	}
}


void change_view(void) {
	Window temp=main_win;
	static char *n[]={"2D" ,"3D","Array","Toon"};
	static char key[]="23at";
	char ch;
	int i;

	ch=(char)pop_up_list(&temp,"Axes",n,key,4,5,0,10,13*DCURY+8,
						 view_hint,info_pop,info_message);
	for(i=0;i<4;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<4) {
		run_the_commands(M_V2+i);
	}
}


void do_windows(void) {
	int i;
	char ch;
	static char *list[]={"(C)reate","(K)ill all","(D)estroy","(B)ottom",
						 "(A)uto","(M)anual","(S)imPlot On"};
	static char *list2[]={"(C)reate","(K)ill all","(D)estroy","(B)ottom",
						  "(A)uto","(M)anual","(S)imPlot Off"};
	static char key[]="ckdbams";
	static char title[]="Make window";
	Window temp=main_win;
	if(SimulPlotFlag==0) {
		ch=(char)pop_up_list(&temp,title,list,key,7,11,0,10,14*DCURY+8,
							 half_hint,info_pop,info_message);
	} else {
		ch=(char)pop_up_list(&temp,title,list2,key,7,11,0,10,14*DCURY+8,
							 half_hint,info_pop,info_message);
	}
	for(i=0;i<7;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<7) {
		run_the_commands(M_MC+i);
	}
}


void add_a_curve(void) {
	int com=-1;
	static char *na[]={"(A)dd curve","(D)elete last","(R)emove all",
					   "(E)dit curve", "(P)ostscript","S(V)G","(F)reeze","a(X)es opts",
					   "exp(O)rt data", "(C)olormap"};
	static char *nc[]={"(N)ormal","(P)eriodic","(H)ot","(C)ool","(B)lue-red",
					   "(G)ray","c(U)behelix"};
	static char *nf[]={"(F)reeze","(D)elete","(E)dit","(R)emove all","(K)ey",
					   "(B)if.Diag","(C)lr. BD","(O)n freeze"};
	static char *nf2[]={"(F)reeze","(D)elete","(E)dit","(R)emove all","(K)ey",
						"(B)if.Diag","(C)lr. BD","(O)ff freeze"};
	static char *nk[]={"(N)o key","(K)ey"};
	static char keya[]="adrepvfxoc";
	static char keyc[]="nphcbgu";
	static char keyf[]="fderkbco";
	static char keyk[]="nk";
	Window temp=main_win;
	char ch;
	int i,j;
	ch=(char)pop_up_list(&temp,"Curves",na,keya,10,15,0,10,8*DCURY+8,
						 graf_hint,info_pop,info_message);
	for(i=0;i<10;i++) {
		if(ch==keya[i]) {
			break;
		}
	}
	if(i==6) {
		if(AutoFreezeFlag==0) {
			ch=(char)pop_up_list(&temp,"Freeze",nf,keyf,8,15,0,10,8*DCURY+8,
								 frz_hint,info_pop,info_message);
		} else {
			ch=(char)pop_up_list(&temp,"Freeze",nf2,keyf,8,15,0,10,8*DCURY+8,
								 frz_hint,info_pop,info_message);
		}
		for(j=0;j<8;j++) {
			if(ch==keyf[j]) {
				break;
			}
		}
		if(j==4) {
			ch=(char)pop_up_list(&temp,"Key",nk,keyk,2,9,0,10,8*DCURY+8,
								 no_hint,info_pop,info_message);
			if(ch==keyk[0]) {
				com=M_GFKN;
			}
			if(ch==keyk[1]) {
				com=M_GFKK;
			}
		} else {
			if(j>=0&&j<8) {
				com=M_GFF+j;
			}
		}
	} else {
		if(i==9) {
			ch=(char)pop_up_list(&temp,"Colormap",nc,keyc,7,15,0,10,8*DCURY+8,
								 cmap_hint,info_pop,info_message);
			for(j=0;j<7;j++) {
				if(ch==keyc[j]) {
					break;
				}
			}
			if(j>=0&&j<7) {
				com=M_GCN+j;
			}
		} else {
			if(i>=0&&i<10) {
				com=M_GA+i;
			}
		}
	}
	run_the_commands(com);
}


void do_movie(void) {
	int i;
	char ch;
	int nkc=6;
	static char *list[]={"(C)apture","(R)eset","(P)layback","(A)utoplay","(S)ave","(M)ake AniGif","(X)tra"};
	static char key[]="crpasmx";
	Window temp=main_win;
	ch=(char)pop_up_list(&temp,"Kinescope",list,key,nkc,11,0,10,8*DCURY+8,
						 kin_hint,info_pop,info_message);
	for(i=0;i<nkc;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<nkc) {
		run_the_commands(i+M_KC);
	}
}


void do_torus(void) {
	Window temp=main_win;
	static char *n[]={"(A)ll","(N)one","(C)hoose"};
	static char key[]="anc";
	char ch;
	int i;
	ch=(char)pop_up_list(&temp,"Torus",n,key,3,9,1-TORUS,10,4*DCURY+8,
						 phas_hint,info_pop,info_message);
	for(i=0;i<3;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<3) {
		run_the_commands(M_AA+i);
	}
}


void window_zoom(void) {
	static char *n[]={"(W)indow","(Z)oom In", "Zoom (O)ut", "(F)it","(D)efault","(S)croll"};
	static char key[]="wzofds";
	char ch;
	int i;
	Window temp=main_win;
	ch=(char)pop_up_list(&temp,"Window",n,key,6,13,0,10,13*DCURY+8,
						 wind_hint,info_pop,info_message);
	for(i=0;i<6;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<6) {
		run_the_commands(M_WW+i);
	}
}


void direct_field(void) {
	int i;
	static char *n[]={"(D)irect Field","(F)low","(N)o dir. fld.","(C)olorize",
					  "(S)caled Dir.Fld"};
	static char key[]="dfncs";
	char ch;
	Window temp=main_win;
	ch=(char)pop_up_list(&temp,"Two-D Fun",n,key,5,18,0,10,6*DCURY+8,
						 flow_hint,info_pop,info_message);
	for(i=0;i<5;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<5) {
		run_the_commands(M_DD+i);
	}
}


void new_clines(void) {
	int i;
	Window temp=main_win;
	static char *n[]={"(N)ew","(R)estore","(A)uto","(M)anual","(F)reeze","(S)ave"};
	static char key[]="nramfs";
	char ch;
	ch=(char)pop_up_list(&temp,"Nullclines",n,key,6,10,0,10,6*DCURY+8,
						 null_hint,info_pop,info_message);
	for(i=0;i<6;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<6) {
		run_the_commands(M_NN+i);
	}
}


void froz_cline_stuff(void) {
	Window temp=main_win;
	static char *n[]={"(F)reeze","(D)elete all","(R)ange","(A)nimate"};
	static char key[]="fdra";
	char ch;
	int i;
	ch=(char)pop_up_list(&temp,"Freeze cline",n,key,4,10,0,10,6*DCURY+8,
						 null_freeze,info_pop,info_message);
	for(i=0;i<4;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>=0&&i<4) {
		run_the_commands(M_NFF+i);
	}
}


void find_equilibrium(void) {
	int i;
	static char *n[]={"(G)o","(M)ouse","(R)ange","monte(C)ar"};
	static char key[]="gmrc";
	char ch;
	Window temp=main_win;
	ch=(char)pop_up_list(&temp,"Equilibria",n,key,4,12,1,10,6*DCURY+8,
						 sing_hint,info_pop,info_message);
	if(ch==ESC) {
		return;
	}
	for(i=0;i<4;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	if(i>-1&&i<4) {
		run_the_commands(i+M_SG);
	}
}


void ini_data_menu(void) {
	int i;
	Window temp=main_win;
	static char *n[]={"(R)ange","(2)par range","(L)ast","(O)ld","(G)o","(M)ouse","(S)hift","(N)ew",
					  "s(H)oot","(F)ile","form(U)la","m(I)ce","DAE guess",
					  "(B)ackward"};
	static char key[]="r2logmsnhfuidb";
	char ch;
	ch= (char)pop_up_list(&temp,"Integrate",n,key,14,13,3,10,3*DCURY+8,
						  ic_hint,info_pop,info_message);
	if(ch==ESC) {
		return;
	}
	for(i=0;i<14;i++) {
		if(ch==key[i]) {
			break;
		}
	}
	run_the_commands(i);
}


void new_param(void) {
	run_the_commands(M_P);
}


void clear_screens(void) {
	run_the_commands(M_EE);
}


void x_vs_t(void) {
	run_the_commands(M_X);
}


void redraw_them_all(void) {
	run_the_commands(M_R);
}

void get_3d_par(void) {
	run_the_commands(M_3);
}

/* --- Static functions --- */
static void do_file_com(int com) {
	switch(com) {
	case M_FT:
		do_transpose();
		break;
	case M_FG:
		get_intern_set();
		break;
	case M_FI:
		TipsFlag=1-TipsFlag;
		break;
	case M_FP:
		make_txtview();
		break;
	case M_FW:
		do_lunch(WRITEM);
		break;
	case M_FS:
		file_inf();
		break;
	case M_FA:
#ifdef AUTO
		do_auto_win();
#endif
		break;
	case M_FC:
		q_calc();
		break;
	case M_FR:
		do_lunch(READEM);
		break;
	case M_FB:
		tfBell=1-tfBell;
		break;
	case M_FH:
		break;
	case M_FX:
		edit_xpprc();
		break;
	case M_FU:
		do_tutorial();
		break;
	case M_FQ:
		if(yes_no_box())bye_bye();
		break;
	case M_FER:
		edit_rhs();
		break;
	case M_FEF:
		edit_functions();
		break;
	case M_FES:
		save_as();
		break;
	case M_FEL:
		load_new_dll();
		break;
	case M_FL:
		clone_ode();
		break;
	}
}


static void get_pmap_pars_com(int l) {
	static char mkey[]="nsmp";
	char ch;
	static char *n[]={"*0Variable","Section","Direction (+1,-1,0)","Stop on sect(y/n)"};
	char values[4][MAX_LEN_SBOX];
	static char *yn[]={"N","Y"};
	int status;
	char n1[15];
	int i1=POIVAR;
	ch=mkey[l];
	POIMAP=0;
	if(ch=='s') {
		POIMAP=1;
	}
	if(ch=='m') {
		POIMAP=2;
	}
	if(ch=='p') {
		POIMAP=3;
	}
	if(POIMAP==0) {
		return;
	}
	ind_to_sym(i1,n1);
	sprintf(values[0],"%s",n1);
	sprintf(values[1],"%.16g",POIPLN);
	sprintf(values[2],"%d",POISGN);
	sprintf(values[3],"%s",yn[SOS]);
	status=do_string_box(4,4,1,"Poincare map",n,values,45);
	if(status!=0) {
		find_variable(values[0],&i1);
		if(i1<0) {
			POIMAP=0;
			err_msg("No such section");
			return;
		}
		POIVAR=i1;
		POISGN=atoi(values[2]);
		if(values[3][0]=='Y' || values[3][0]=='y') {
			SOS=1;
		} else {
			SOS=0;
		}
		POIPLN=atof(values[1]);
	}
}


static void quick_num(int com) {
	char key[]="tsrdnviobec";
	if(com>=0 && com<11) {
		get_num_par(key[com]);
	}
}


static void set_col_par_com(int i) {
	int j,ivar;
	double temp[2];
	float maxder=0.0,minder=0.0,sum=0.0;
	char ch,name[20];
	MyGraph->ColorFlag=i;
	/* set color to black/white */
	if(MyGraph->ColorFlag==0) {
		return;
	}
	if(MyGraph->ColorFlag==2) {
		ind_to_sym(MyGraph->ColorValue,name);
		new_string("Color via:",name);
		find_variable(name,&ivar);
		if(ivar>=0) {
			MyGraph->ColorValue=ivar;
		} else {
			err_msg("No such quantity!");
			MyGraph->ColorFlag=0;
			return;
		}
	}

	/*   This will be uncommented    ..... */
	ch=TwoChoice("(O)ptimize","(C)hoose","Color","oc");

	if(ch=='c') {
		temp[0]=MyGraph->min_scale;
		temp[1]=MyGraph->min_scale+MyGraph->color_scale;
		new_float("Min :",&temp[0]);
		new_float("Max :",&temp[1]);
		if(temp[1]>temp[0]				&&
		   ((MyGraph->ColorFlag==2)	||
			(MyGraph->ColorFlag==1		&&
			 temp[0]>=0.0))) {
			MyGraph->min_scale=temp[0];
			MyGraph->color_scale=(temp[1]-temp[0]);
		} else {
			err_msg("Min>=Max or Min<0 error");
		}
		return;
	}
	if(MyGraph->ColorFlag==1) {
		if(storind<2) {
			return;
		}
		maxder=0.0;
		minder=1.e20;
		for(i=1;i<my_browser.maxrow;i++) {
			sum=0.0;
			for(j=0;j<NODE;j++) {
				sum+=(float)fabs((double)(my_browser.data[1+j][i]-my_browser.data[1+j][i-1]));
			}
			if(sum<minder) {
				minder=sum;
			}
			if(sum>maxder) {
				maxder=sum;
			}
		}
		if(minder>=0.0 && maxder>minder) {
			MyGraph->color_scale=(maxder-minder)/(fabs(DELTA_T*NJMP));
			MyGraph->min_scale=minder/(fabs(DELTA_T*NJMP));
		}
	} else {
		get_max(MyGraph->ColorValue,&temp[0],&temp[1]);
		MyGraph->min_scale=temp[0];
		MyGraph->color_scale=(temp[1]-temp[0]);
		if(MyGraph->color_scale==0.0) {
			MyGraph->color_scale=1.0;
		}
	}
}
