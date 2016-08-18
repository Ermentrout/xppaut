AM_CFLAGS = -std=c99 -Wall -pedantic
AM_CPPFLAGS = \
		-D_XOPEN_SOURCE=600 \
		-DAUTO \
		-DMYSTR1=$(MAJORVER) \
		-DMYSTR2=$(MINORVER) \
		-Iauto \
		-Isolver \
		-Isolver/cvodesrc

bin_PROGRAMS = \
		xppaut

lib_LIBRARIES = \
		libxpp.a

noinst_HEADERS = \
		abort.h \
		adj2.h \
		aniparse.h \
		arrayplot.h \
		array_print.h \
		autevd.h \
		autlim.h \
		auto_nox.h \
		auto_x11.h \
		axes2.h \
		browse.h \
		calc.h \
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
		edit_rhs.h \
		eig_list.h \
		extra.h \
		fftn.h \
		find_fixedPoint.h \
		flags.h \
		form_ode.h \
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
		main.h \
		many_pops.h \
		markov.h \
		menudrive.h \
		menu.h \
		menus.h \
		myfonts.h \
		mykeydef.h \
		my_ps.h \
		my_rhs.h \
		my_svg.h \
		newhome.h \
		newpars.h \
		nullcline.h \
		numerics.h \
		parserslow.h \
		pop_list.h \
		pp_shoot.h \
		read_dir.h \
		rubber.h \
		scrngif.h \
		simplenet.h \
		storage.h \
		struct.h \
		tabular.h \
		toons.h \
		torus.h \
		tutor.h \
		txtread.h \
		userbut.h \
		xAuto.h \
		xpplim.h \
		$(bitmap_headers) \
		$(solver_headers) \
		$(util_headers)

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

solver_headers = \
		solver/adams.h \
		solver/backeuler.h \
		solver/cv2.h \
		solver/discrete.h \
		solver/dormpri.h \
		solver/euler.h \
		solver/gear.h \
		solver/heun.h \
		solver/rosenbrock.h \
		solver/runge_kutta.h \
		solver/stiff.h \
		solver/symplect.h \
		solver/volterra2.h

util_headers = \
		util/matrixalgebra.h \
		util/strutil.h \
		util/timeutil.h

xppaut_LDADD = \
		solver/cvodesrc/libcvode.a \
		auto/libauto.a

xppaut_SOURCES = \
		abort.c \
		adj2.c \
		aniparse.c \
		arrayplot.c \
		array_print.c \
		autevd.c \
		auto_nox.c \
		auto_x11.c \
		autpp.c \
		axes2.c \
		browse.c \
		calc.c \
		color.c \
		comline.c \
		dae_fun.c \
		delay_handle.c \
		del_stab.c \
		derived.c \
		dialog_box.c \
		diagram.c \
		do_fit.c \
		edit_rhs.c \
		eig_list.c \
		extra.c \
		fftn.c \
		find_fixedPoint.c \
		flags.c \
		form_ode.c \
		ggets.c \
		gogoauto.c \
		graf_par.c \
		graphics.c \
		histogram.c \
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
		menus.c \
		my_ps.c \
		my_rhs.c \
		my_svg.c \
		nullcline.c \
		numerics.c \
		parserslow.c \
		pop_list.c \
		pp_shoot.c \
		read_dir.c \
		rubber.c \
		scrngif.c \
		simplenet.c \
		storage.c \
		tabular.c \
		torus.c \
		txtread.c \
		userbut.c \
		$(solver_sources) \
		$(util_sources)

libxpp_a_SOURCES = \
		$(xppaut_SOURCES)

solver_sources = \
		solver/adams.c \
		solver/backeuler.c \
		solver/cv2.c \
		solver/discrete.c \
		solver/dormpri.c \
		solver/euler.c \
		solver/gear.c \
		solver/heun.c \
		solver/rosenbrock.c \
		solver/runge_kutta.c \
		solver/stiff.c \
		solver/symplect.c \
		solver/volterra2.c

util_sources = \
		util/matrixalgebra.c \
		util/strutil.c \
		util/timeutil.c
