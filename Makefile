# Copyright (C) 1990-2016 Bard Ermentrout & Daniel Dougherty & Robert McDougal
# Edited for Debian GNU/Linux.
#BINDIR = /usr/X11R6/bin
DESTDIR =
BINDIR = /usr/local/bin
DOCDIR = /usr/share/doc/xppaut
MANDIR = /usr/local/man/man1
# End Debian Edit
#################################
#
VERSION=8.0
MAJORVER=8.0
MINORVER=0

ODES=ode/*.ode ode/*.ani canonical/*.* tstauto/*.ode
DOC=xpp_doc.pdf xpp_sum.pdf install.pdf tree.pdf 
HELP=help/*.html
# Standard C compiler
#CC= cc
# Use Gnu compiler
#CC= gcc
# Use c99 to compile according to newer ISO C standards (IEEE Std 1003.1-2001)
CC= gcc -std=c99 -pedantic -D_XOPEN_SOURCE=600 -Wall  
#AUTLIBS=  -lsundials_cvode -lX11 -lm 
#These are the 32bit compat libraries.
#AUTLIBS= -lf2c -lsundials_cvode -lX11 -lm 
#OTHERLIBS= libcvode.a libf2cm.a 
#OTHERLIBS= 
#
# Free BSD 10.0 - to'th Marmaduke Woodman
#CC= clang -std=c99 -pedantic -D_XOPEN_SOURCE=600 -Wall
#CFLAGS= -g -pedantic -O -DNOERRNO -DNON_UNIX_STDIO -DAUTO -DCVODE_YES  -DHAVEDLL -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER)  -I/usr/local/include
#LDFLAGS= -m64 -L/usr/local/lib/
#LIBS=  -lX11 -l
################################## 
# Standard Linux distributions   #
##################################
#CFLAGS=   -g -O -DAUTO -DCVODE_YES -DHAVEDLL -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER)  -I/usr/X11R6/include
#CFLAGS=   -g -O -m32 -DNON_UNIX_STDIO -DAUTO -DCVODE_YES  -DHAVEDLL -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER)  -I/usr/X11R6/include
#CFLAGS=   -g -O -m64 -DNOERRNO -DNON_UNIX_STDIO -DAUTO -DCVODE_YES  -DHAVEDLL -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER)  -I/usr/include/X11

CFLAGS= -g -pedantic -O2     -DNOERRNO -DNON_UNIX_STDIO -DAUTO -DCVODE_YES  -DHAVEDLL -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER)  -I/opt/X11/include 
#LDFLAGS=  -m64 -L/usr/lib -L/usr/lib64
#LDFLAGS=  -m32 -L/usr/X11R6/lib 
LDFLAGS=  -L/usr/X11/lib 
LIBS=  -lX11 -lm -ldl   
# NOTE: Recent (RedHat 8) versions of GCC seem to no longer have
# the integer errno, so compile with the -DNOERRNO option as well
#
# some errors with ctype on newer machines HP ITANIUM, eg can be fixed
# with the -DWCTYPE 
#
#   64 Bit machines must use a different parser code.
#   In the OBJECTS section, replace parser2.o with parserslow2.o
# it is not really that much slower so don't have a cow
#
#################################
# MACOSX                        #
#################################
#CFLAGS=   -g -O -DMACOSX -DAUTO -DCVODE_YES  -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER) -I/usr/X11R6/include
# LIBS=  -lX11 -lm
# LDFLAGS=  -L/usr/X11R6/lib
#
#################################
# CYGWIN                        #
#################################
# CFLAGS=   -O -DBGR -DNON_UNIX_STDIO -DNORAND48 -DSTRUPR -DNOERRNO -DAUTO -DCVODE_YES -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER)  -I/usr/X11R6/include
#LIBS= -lX11 -lm
#LDFLAGS= -L/usr/X11R6/lib 
#
#################################
# SPARC                         #
#################################
# CFLAGS=  -O -DAUTO  -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER)  -DCVODE_YES -I/usr/openwin/include 
# LDFLAGS= -L/usr/openwin/lib
# LIBS= -lm -lX11 
#
################################
# SUNPro compiler              #
################################
# CFLAGS=  -O -DSUNPRO -DAUTO -DCVODE_YES  -I/usr/openwin/include
# LDFLAGS= -L/usr/openwin/lib
# LIBS= -lm -lX11 
#
################################
# SGI                          #
################################
#CFLAGS=   -g -O -DAUTO -DCVODE_YES  -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER) -I/usr/X11R6/include
#  Old SGIs without gcc
# 
#CFLAGS= -cckr -O -DAUTO -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER) -DCVODE_YES
# LDFLAGS=  -L/usr/X11R6/lib
# LIBS= -lm -lX11 
#
###############################
# HP                          #
###############################
#CFLAGS=   -g -O -DAUTO -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER) -DCVODE_YES  -DHAVEDLL -I/usr/X11R6/include
# LDFLAGS=  -L/usr/X11R6/lib
# LIBS= -lm -lX11 -ldl
#
###############################
# DEC OSF                     #
###############################
#
#CFLAGS= -O -DAUTO -DCVODE_YES -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER) -Olimit 1000
#
###############################################################  
#             You  can stop messing with it now, the rest is
#             probably OkeyDokey
###############################################################
HEADERS = abort.h browse.h do_fit.h iterativ.h my_rhs.h shoot.h \
	  adj2.h calc.h  dormpri.h kbs.h  my_svg.h simplenet.h \
	  aniparse.h choice_box.h edit_rhs.h kinescope.h myfonts.h \
	  aniparse_avi.h color.h  eig_list.h llnlmath.h mykeydef.h stiff.h \
	  autevd.h cvband.h fftn.h  local.h  nullcline.h tabular.h \
	  autlim.h cvdense.h flags.h  lunch-new.h numerics.h toons.h \
	  auto.h  cvdiag.h form_ode.h macdirent.h odesol2.h torus.h \
	  auto_c.h cvode.h  gear.h  macsysdirent.h parser.h tutor.h \
	  auto_def2.h cvspgmr.h getvar.h main.h  parserslow.h txtread.h \
	  auto_define.h dae_fun.h ggets.h  many_pops.h phsplan.h userbut.h \
	  auto_f2c.h del_stab.h graf_par.h markov.h pop_list.h vector.h \
	  auto_nox.h delay_handle.h graphics.h menu.h  pp_shoot.h volterra.h \
	  auto_types.h dense.h  help_defs.h menudrive.h read_dir.h volterra2.h \
	  auto_x11.h derived.h histogram.h menus.h  rubber.h xAuto.h \
	  axes2.h  diagram.h init_conds.h my_pars.h run_auto.h xpplim.h \
	  band.h  dialog_box.h integrate.h my_ps.h  scrngif.h spgmr.h \
	  arrayplot.h cv2.h  extra.h load_eqn.h newpars.h struct.h \
	  array_print.h comline.h endfile.h llnltyps.h newhome.h storage.h 	

BITMAPS = aniwin.bitmap alert.bitmap bc.bitmap eqns.bitmap info.bitmap ic.bitmap pp.bitmap\
	  array.bitmap browse.bitmap equilib.bitmap lines.bitmap txtview.bitmap\
	  auto.bitmap delay.bitmap graph.bitmap param.bitmap filebrowse.bitmap \
	  lineup.bitmap linedn.bitmap pageup.bitmap pagedn.bitmap home.bitmap start.bitmap  

SOURCES = abort.c  cabs.c  delay_handle.c graf_par.c my_ps.c  simplenet.c \
	  adj2.c  calc.c  dense.c  graphics.c my_rhs.c spgmr.c \
	  aniparse.c choice_box.c derived.c histogram.c my_svg.c stiff.c \
	  array_print.c color.c  diagram.c i_dnnt.c nullcline.c storage.c \
	  arrayplot.c comline.c dialog_box.c i_nint.c numerics.c tabular.c \
	  autevd.c conpar2.c do_fit.c init_conds.c odesol2.c torus.c \
	  autlib1.c cv2.c  dormpri.c integrate.c parserslow2.c txtread.c \
	  autlib2.c cvband.c edit_rhs.c iterativ.c pop_list.c userbut.c \
	  autlib3.c cvdense.c eig_list.c kinescope.c pow_dd.c vector.c \
	  autlib4.c cvdiag.c eispack.c llnlmath.c pow_di.c volterra2.c \
	  autlib5.c cvode.c  extra.c  load_eqn.c pow_ii.c worker2.c \
	  auto_nox.c cvspgmr.c fftn.c  lunch-new.c pp_shoot.c z_abs.c \
	  auto_x11.c d_imag.c flags.c  main.c  r_lg10.c z_exp.c \
	  autpp.c  d_lg10.c form_ode.c many_pops.c read_dir.c z_log.c \
	  axes2.c  d_sign.c gear.c  markov.c rubber.c \
	  band.c  dae_fun.c ggets.c  menu.c  scrngif.c \
	  browse.c del_stab.c gogoauto.c menudrive.c setubv2.c \


#
# WARNING: For 64 bit machines replace parser2.o with parserslow2.o 
#          in the OBJECTS section!
#
OBJECTS = abort.o  cabs.o  delay_handle.o graf_par.o my_ps.o  simplenet.o \
	  adj2.o  calc.o  dense.o  graphics.o my_rhs.o spgmr.o \
	  aniparse.o choice_box.o derived.o histogram.o my_svg.o stiff.o \
	  array_print.o color.o  diagram.o i_dnnt.o nullcline.o storage.o \
	  arrayplot.o comline.o dialog_box.o i_nint.o numerics.o tabular.o \
	  autevd.o conpar2.o do_fit.o init_conds.o odesol2.o torus.o \
	  autlib1.o cv2.o  dormpri.o integrate.o parserslow2.o txtread.o \
	  autlib2.o cvband.o edit_rhs.o iterativ.o pop_list.o userbut.o \
	  autlib3.o cvdense.o eig_list.o kinescope.o pow_dd.o vector.o \
	  autlib4.o cvdiag.o eispack.o llnlmath.o pow_di.o volterra2.o \
	  autlib5.o cvode.o  extra.o  load_eqn.o pow_ii.o worker2.o \
	  auto_nox.o cvspgmr.o fftn.o  lunch-new.o pp_shoot.o z_abs.o \
	  auto_x11.o d_imag.o flags.o  main.o  r_lg10.o z_exp.o \
	  autpp.o  d_lg10.o form_ode.o many_pops.o read_dir.o z_log.o \
	  axes2.o  d_sign.o gear.o  markov.o rubber.o \
	  band.o  dae_fun.o ggets.o  menu.o  scrngif.o \
	  browse.o del_stab.o gogoauto.o menudrive.o setubv2.o \

#If you are using sundials cvode libraries then you may comment out following above
#	cvode.o  cvdense.o  dense.o  cvband.o \
#        band.o  cvdiag.o  cvspgmr.o  spgmr.o  \
#        iterativ.o  vector.o  llnlmath.o cv2.o
######################################################################
#
#
xppaut: $(OBJECTS)
#
###########################################################
##  shared library method - delete my_fun.o from OBJECTS
#	$(CC) -DAUTO -o xppaut $(OBJECTS) $(AUTOOBJ) $(LDFLAGS) libcvode.a libf2cm.a libmyfun.so $(LIBS) 	
####  You made your own libraries but installed locally (2)
#	$(CC) -DAUTO -o xppaut $(OBJECTS) $(AUTOOBJ) $(LDFLAGS) libcvode.a libI77.a libF77.a $(LIBS) 	
##
####  You have previously installed the f2c libraries
#$(CC) -DAUTO -o xppaut $(OBJECTS) $(AUTOOBJ) $(LDFLAGS) $(LIBS) libcvode.a $(AUTLIBS) 
#  
####  You have previously installed the f2c libraries and sundials
#	 $(CC) -DAUTO -o xppaut $(OBJECTS) $(AUTOOBJ) $(LDFLAGS) $(LIBS) $(AUTLIBS) 
####  You have previously installed the sundials
	 #$(CC) -DAUTO -o xppaut $(OBJECTS) $(AUTOOBJ) $(LDFLAGS) $(LIBS) $(AUTLIBS) ./libf2cm.a
	 $(CC) -DAUTO -o xppaut $(OBJECTS) $(LDFLAGS) $(LIBS) 
tarfile:
	tar zcvf xppaut$(VERSION).tgz $(SOURCES) $(HEADERS) $(BITMAPS) default.opt \
	 xpp_doc.tex README Makefile.old Makefile Makefile.64 Makefile.nice  Makefile.lib Makefile.avi Makefile.old Makefile.w7b \
	ode/*.*  tstauto/*.ode cuda/README cuda/*.* canonical/*.* xpp_doc.pdf xpp_sum.tex xpp_sum.pdf  xpp_doc.tex tree.tex tree.pdf  \
	xppaut.1 install*.html\
        help/*.html \
	help/odes/*.ode help/odes/*.c install.pdf install.tex LICENSE HISTORY \
	Makefile.s2x sbml2xpp.c *.xbm *.jar changes.txt \
	xppaut-stylesheet.css	xppaut-stylesheet2.css animsvgwww


#   
# 
###################################################################
################################################################
#
# Edited for Debian GNU/Linux.
install: xppaut 
# Make necessary installation directories
	mkdir -p $(DESTDIR)$(BINDIR)
	mkdir -p $(DESTDIR)$(DOCDIR)/html
	mkdir -p $(DESTDIR)$(DOCDIR)/examples
	mkdir -p $(DESTDIR)$(DOCDIR)/xbm
	mkdir -p $(DESTDIR)$(MANDIR)
# Put everything home
	strip xppaut
	install -m 755 xppaut $(DESTDIR)$(BINDIR)
	cp -r ode* $(DESTDIR)$(DOCDIR)/examples
	cp -r help/* $(DESTDIR)$(DOCDIR)/html
	cp README *.pdf $(DESTDIR)$(DOCDIR)
	cp *.xbm $(DESTDIR)$(DOCDIR)/xbm 
	cp xppaut.1 $(DESTDIR)$(MANDIR)
# End Debian Ed
uninstall: 
# Remove everything you created
	rm $(DESTDIR)$(BINDIR)/xppaut
	rm -r $(DESTDIR)$(DOCDIR)
	rm -r $(DESTDIR)$(MANDIR)/xppaut.1
# End Debian Ed
##################################################
#    Make s stand alone library -must link your rhs
###################################################
####################################################
#  tar file
####################################################
##############################################
#  pack up a binary
##############################################
binary:
	tar zvcf binary.tgz xppaut *.pdf *.xbm xppaut.1 canonical/* $(ODES) $(DOC) $(HELP) README HISTORY LICENSE
##############################################
#  clean
##############################################
clean:
	rm -f *.o xppaut 
#######################################################
#  Documentation
#######################################################
xppdoc:     
	 latex xpp_doc
	 latex xpp_doc
	 latex xpp_doc		
	 dvips -o xpp_doc.ps  xpp_doc
	 ps2pdf xpp_doc.ps
	 latex xpp_sum
	 latex xpp_sum
	 dvips -o xpp_sum.ps  xpp_sum
	  ps2pdf xpp_sum.ps
	  latex install.tex
	 latex install.tex
	 dvips -o install.ps  install
	  ps2pdf install.ps
         
