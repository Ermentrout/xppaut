AM_CFLAGS = -std=c99 -Wall -pedantic
AM_CPPFLAGS = \
	-D_XOPEN_SOURCE=600 \
	-DNOERRNO \
	-DNON_UNIX_STDIO \
	-DAUTO \
	-DCVODE_YES \
	-DHAVEDLL \
	-DMYSTR1=$(MAJORVER) \
	-DMYSTR2=$(MINORVER) \
	-I$(top_srcdir)/cvodesrc

bin_PROGRAMS = \
	xppaut

lib_LIBRARIES = \
	libxpp.a

check_PROGRAMS = \
	libexample.so

noinst_HEADERS = \
	abort.h \
	adj2.h \
	aniparse.h \
	arrayplot.h \
	array_print.h \
	autevd.h \
	autlim.h \
	auto.h \
	auto_c.h \
	auto_def2.h \
	auto_f2c.h \
	auto_nox.h \
	auto_types.h \
	auto_x11.h \
	axes2.h \
	browse.h \
	calc.h \
	choice_box.h \
	color.h \
	comline.h \
	config.h \
	dae_fun.h \
	delay_handle.h \
	del_stab.h \
	derived.h \
	diagram.h \
	dialog_box.h \
	do_fit.h \
	dormpri.h \
	edit_rhs.h \
	eig_list.h \
	extra.h \
	fftn.h \
	flags.h \
	form_ode.h \
	gear.h \
	getvar.h \
	ggets.h \
	graf_par.h \
	graphics.h \
	help_defs.h \
	histogram.h \
	init_conds.h \
	integrate.h \
	kbs.h \
	kinescope.h \
	load_eqn.h \
	lunch-new.h \
	macdirent.h \
	macsysdirent.h \
	main.h \
	many_pops.h \
	markov.h \
	menudrive.h \
	menu.h \
	menus.h \
	myfonts.h \
	mykeydef.h \
	my_pars.h \
	my_ps.h \
	my_rhs.h \
	my_svg.h \
	newhome.h \
	newpars.h \
	nullcline.h \
	numerics.h \
	odesol2.h \
	parserslow.h \
	phsplan.h \
	pop_list.h \
	pp_shoot.h \
	read_dir.h \
	rubber.h \
	run_auto.h \
	scrngif.h \
	shoot.h \
	simplenet.h \
	stiff.h \
	storage.h \
	struct.h \
	tabular.h \
	toons.h \
	torus.h \
	tutor.h \
	txtread.h \
	userbut.h \
	volterra2.h \
	volterra.h \
	xAuto.h \
	xpplim.h \
	$(bitmap_headers)

bitmap_headers = \
	bitmap/alert.bitmap \
	bitmap/aniwin.bitmap \
	bitmap/array.bitmap \
	bitmap/auto.bitmap \
	bitmap/bc.bitmap \
	bitmap/browse.bitmap \
	bitmap/delay.bitmap \
	bitmap/eqns.bitmap \
	bitmap/equilib.bitmap \
	bitmap/filebrowse.bitmap \
	bitmap/graph.bitmap \
	bitmap/home.bitmap \
	bitmap/ic.bitmap \
	bitmap/info.bitmap \
	bitmap/linedn.bitmap \
	bitmap/lines.bitmap \
	bitmap/lineup.bitmap \
	bitmap/pagedn.bitmap \
	bitmap/pageup.bitmap \
	bitmap/param.bitmap \
	bitmap/pp.bitmap \
	bitmap/start.bitmap \
	bitmap/txtview.bitmap

xppaut_LDADD = \
	$(top_builddir)/cvodesrc/libcvode.a

xppaut_SOURCES = \
	abort.c \
	adj2.c \
	aniparse.c \
	arrayplot.c \
	array_print.c \
	axes2.c \
	browse.c \
	cabs.c \
	calc.c \
	choice_box.c \
	color.c \
	comline.c \
	conpar2.c \
	dae_fun.c \
	delay_handle.c \
	del_stab.c \
	derived.c \
	dialog_box.c \
	d_imag.c \
	d_lg10.c \
	do_fit.c \
	dormpri.c \
	d_sign.c \
	edit_rhs.c \
	eig_list.c \
	eispack.c \
	extra.c \
	fftn.c \
	flags.c \
	form_ode.c \
	gear.c \
	ggets.c \
	gogoauto.c \
	graf_par.c \
	graphics.c \
	histogram.c \
	i_dnnt.c \
	i_nint.c \
	init_conds.c \
	integrate.c \
	kinescope.c \
	load_eqn.c \
	lunch-new.c \
	main.c \
	many_pops.c \
	markov.c \
	menu.c \
	menudrive.c \
	my_ps.c \
	my_rhs.c \
	my_svg.c \
	nullcline.c \
	numerics.c \
	odesol2.c \
	parserslow2.c \
	pop_list.c \
	pow_dd.c \
	pow_di.c \
	pow_ii.c \
	pp_shoot.c \
	read_dir.c \
	r_lg10.c \
	rubber.c \
	scrngif.c \
	setubv2.c \
	simplenet.c \
	stiff.c \
	storage.c \
	tabular.c \
	torus.c \
	txtread.c \
	userbut.c \
	volterra2.c \
	worker2.c \
	z_abs.c \
	z_exp.c \
	z_log.c \
	$(cvode_sources) \
	$(auto_sources)

libxpp_a_SOURCES = \
	abort.c \
	adj2.c \
	aniparse.c \
	arrayplot.c \
	array_print.c \
	axes2.c \
	browse.c \
	cabs.c \
	calc.c \
	choice_box.c \
	color.c \
	comline.c \
	conpar2.c \
	dae_fun.c \
	delay_handle.c \
	del_stab.c \
	derived.c \
	dialog_box.c \
	d_imag.c \
	d_lg10.c \
	do_fit.c \
	dormpri.c \
	d_sign.c \
	edit_rhs.c \
	eig_list.c \
	eispack.c \
	extra.c \
	fftn.c \
	flags.c \
	form_ode.c \
	gear.c \
	ggets.c \
	gogoauto.c \
	graf_par.c \
	graphics.c \
	histogram.c \
	i_dnnt.c \
	i_nint.c \
	init_conds.c \
	integrate.c \
	kinescope.c \
	load_eqn.c \
	lunch-new.c \
	main.c \
	many_pops.c \
	markov.c \
	menu.c \
	menudrive.c \
	my_ps.c \
	my_rhs.c \
	my_svg.c \
	nullcline.c \
	numerics.c \
	odesol2.c \
	parserslow2.c \
	pop_list.c \
	pow_dd.c \
	pow_di.c \
	pow_ii.c \
	pp_shoot.c \
	read_dir.c \
	r_lg10.c \
	rubber.c \
	scrngif.c \
	setubv2.c \
	simplenet.c \
	stiff.c \
	storage.c \
	tabular.c \
	torus.c \
	txtread.c \
	userbut.c \
	volterra2.c \
	worker2.c \
	z_abs.c \
	z_exp.c \
	z_log.c \
	$(cvode_sources) \
	$(auto_sources)

auto_sources = \
	autlib1.c \
	autlib2.c \
	autlib3.c \
	autlib4.c \
	autlib5.c \
	autevd.c \
	autpp.c \
	diagram.c \
	auto_nox.c \
	auto_x11.c

libexample_so_CFLAGS = -fPIC
libexample_so_LDFLAGS = -shared

libexample_so_SOURCES = \
	funexample.c
