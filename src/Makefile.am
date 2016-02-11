AM_CFLAGS = -std=c99 -Wall -pedantic
AM_CPPFLAGS = \
	-D_XOPEN_SOURCE=600 \
	-DNOERRNO \
	-DNON_UNIX_STDIO \
	-DAUTO \
	-DCVODE_YES \
	-DHAVEDLL \
	-DMYSTR1=$(MAJORVER) \
	-DMYSTR2=$(MINORVER)

bin_PROGRAMS = \
	xppaut

lib_LIBRARIES = \
	libxpp.a

check_PROGRAMS = \
	libexample.so

noinst_HEADERS = \
	abort.h \
	comline.h \
	extra.h \
	kinescope.h \
	newhome.h \
	stiff.h \
	adj2.h \
	cvband.h \
	f2c.h \
	lio.h \
	newpars.h \
	storage.h \
	aniparse_avi.h \
	cvdense.h \
	fftn.h \
	llnlmath.h \
	nullcline.h \
	struct.h \
	aniparse.h \
	cvdiag.h \
	fio.h \
	llnltyps.h \
	numerics.h \
	sue.h \
	arrayplot.h \
	cvode.h \
	flags.h \
	load_eqn.h \
	odesol2.h \
	tabular.h \
	array_print.h \
	cvspgmr.h \
	fmt.h \
	local.h \
	toons.h \
	autevd.h \
	dae_fun.h \
	form_ode.h \
	lread.h \
	open.h \
	torus.h \
	autlib1.h \
	delay_handle.h \
	fp.h \
	lunch-new.h \
	parser.h \
	txtread.h \
	autlib2.h \
	del_stab.h \
	gear.h \
	macdirent.h \
	parserslow.h \
	uio.h \
	autlim.h \
	dense.h \
	getvar.h \
	macsysdirent.h \
	phsplan.h \
	userbut.h \
	auto_define.h \
	derived.h \
	ggets.h \
	main.h \
	pop_list.h \
	util.h \
	auto.h \
	dfe.h \
	graf_par.h \
	many_pops.h \
	pp_shoot.h \
	vector.h \
	auto_nox.h \
	diagram.h \
	graphics.h \
	markov.h \
	rdfmt.h \
	volterra2.h \
	auto_x11.h \
	dialog_box.h \
	help_defs.h \
	menudrive.h \
	read_dir.h \
	volterra.h \
	axes2.h \
	do_fit.h \
	histogram.h \
	menu.h \
	rubber.h \
	wref.h \
	band.h \
	dormpri.h \
	homsup.h \
	menus.h \
	run_auto.h \
	wrtfmt.h \
	browse.h \
	due.h \
	iio.h \
	myfonts.h \
	scrngif.h \
	wsfe.h \
	calc.h \
	edit_rhs.h \
	init_conds.h \
	mykeydef.h \
	sfe.h \
	wsne.h \
	choice_box.h \
	eig_list.h \
	integrate.h \
	my_pars.h \
	shoot.h \
	xpplim.h \
	close.h \
	endfile.h \
	iterativ.h \
	my_ps.h \
	simplenet.h \
	color.h \
	err.h \
	kbs.h \
	my_rhs.h \
	spgmr.h \
	cv2.h \
	my_svg.h \
	tutor.h \
	$(bitmap_headers)

bitmap_headers = \
	bitmap/aniwin.bitmap \
	bitmap/alert.bitmap \
	bitmap/bc.bitmap \
	bitmap/eqns.bitmap \
	bitmap/info.bitmap \
	bitmap/ic.bitmap \
	bitmap/pp.bitmap \
	bitmap/array.bitmap \
	bitmap/browse.bitmap \
	bitmap/equilib.bitmap \
	bitmap/lines.bitmap \
	bitmap/txtview.bitmap \
	bitmap/auto.bitmap \
	bitmap/delay.bitmap \
	bitmap/graph.bitmap \
	bitmap/param.bitmap \
	bitmap/filebrowse.bitmap \
	bitmap/home.bitmap \
	bitmap/linedn.bitmap \
	bitmap/lineup.bitmap \
	bitmap/pagedn.bitmap \
	bitmap/pageup.bitmap \
	bitmap/start.bitmap

xppaut_SOURCES = \
	main.c \
	ggets.c \
	menu.c \
	rubber.c \
	derived.c \
	many_pops.c \
	pop_list.c \
	graphics.c \
	dialog_box.c \
	numerics.c \
	choice_box.c \
	color.c \
	init_conds.c \
	browse.c \
	kinescope.c \
	axes2.c \
	abort.c \
	parserslow2.c \
	storage.c \
	load_eqn.c \
	form_ode.c \
	odesol2.c \
	gear.c \
	eig_list.c \
	integrate.c \
	delay_handle.c \
	graf_par.c \
	dormpri.c \
	my_ps.c \
	my_svg.c \
	nullcline.c \
	torus.c \
	pp_shoot.c \
	lunch-new.c \
	calc.c \
	adj2.c \
	my_rhs.c \
	read_dir.c\
	volterra2.c \
	tabular.c \
	markov.c \
	histogram.c \
	comline.c \
	edit_rhs.c \
	do_fit.c \
	flags.c \
	del_stab.c \
	stiff.c \
	arrayplot.c \
	array_print.c \
	aniparse.c \
	simplenet.c \
	dae_fun.c \
	fftn.c \
	extra.c \
	scrngif.c \
	nagroutines.c \
	homsup.c \
	txtread.c \
	menudrive.c \
	userbut.c \
	Version.c \
	backspace.c \
	dfe.c \
	due.c \
	iio.c \
	inquire.c \
	rewind.c \
	rsfe.c \
	rdfmt.c \
	sue.c \
	uio.c \
	wsfe.c \
	sfe.c \
	fmt.c \
	lio.c \
	lread.c \
	open.c \
	close.c \
	util.c \
	endfile.c \
	wrtfmt.c \
	wref.c \
	err.c \
	fmtlib.c \
	rsne.c \
	wsne.c \
	fc77.c \
	f2cstart.c \
	$(cvode_sources) \
	$(auto_sources)

libxpp_a_SOURCES = \
	main.c \
	ggets.c \
	menu.c \
	rubber.c \
	derived.c \
	many_pops.c \
	pop_list.c \
	graphics.c \
	dialog_box.c \
	numerics.c \
	choice_box.c \
	color.c \
	init_conds.c \
	browse.c \
	kinescope.c \
	axes2.c \
	abort.c \
	parserslow2.c \
	storage.c \
	load_eqn.c \
	form_ode.c \
	odesol2.c \
	gear.c \
	eig_list.c \
	integrate.c \
	delay_handle.c \
	graf_par.c \
	dormpri.c \
	my_ps.c \
	my_svg.c \
	nullcline.c \
	torus.c \
	pp_shoot.c \
	lunch-new.c \
	calc.c \
	adj2.c \
	read_dir.c \
	volterra2.c \
	tabular.c \
	markov.c \
	histogram.c \
	comline.c \
	edit_rhs.c \
	do_fit.c \
	flags.c \
	del_stab.c \
	stiff.c \
	arrayplot.c \
	array_print.c \
	aniparse.c \
	simplenet.c \
	dae_fun.c \
	fftn.c \
	extra.c \
	scrngif.c \
	nagroutines.c \
	homsup.c \
	txtread.c \
	menudrive.c \
	userbut.c \
	$(cvode_sources) \
	$(auto_sources)

auto_sources = \
	autlib1.c \
	autlib2.c \
	autlib3.c \
	autevd.c \
	run_auto.c \
	autpp.c \
	diagram.c \
	auto_nox.c \
	auto_x11.c \
	flowkm_small.c

cvode_sources = \
	band.c \
	cv2.c \
	cvband.c \
	cvdense.c \
	cvdiag.c \
	cvode.c \
	cvspgmr.c \
	dense.c \
	iterativ.c \
	llnlmath.c \
	spgmr.c \
	vector.c

libexample_so_CFLAGS = -fPIC
libexample_so_LDFLAGS = -shared

libexample_so_SOURCES = \
	funexample.c
