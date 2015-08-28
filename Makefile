# Copyright (C) 1990-2013 Bard Ermentrout & Daniel Dougherty & Robert McDougal
# Edited for Debian GNU/Linux.
#BINDIR = /usr/X11R6/bin
DESTDIR =
BINDIR = /usr/local/bin
DOCDIR = /usr/share/doc/xppaut
MANDIR = /usr/local/man/man1
# End Debian Edit
#################################
#
VERSION=7.0
MAJORVER=7.0
MINORVER=0

ODES=ode/*.ode ode/*.ani
DOC=xpp_doc.pdf xpp_sum.pdf install.pdf
HELP=help/*.html
# Standard C compiler
#CC= cc
# Use Gnu compiler
#CC= gcc
# Use c99 to compile according to newer ISO C standards (IEEE Std 1003.1-2001)
CC= gcc -std=c99 -pedantic -D_XOPEN_SOURCE=600 -Wall  
#AUTLIBS=  -lsundials_cvode -lX11 -lm 
#These are the 32bit compat libraries.
AUTLIBS=  -lm 
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

CFLAGS= -g -pedantic -O -m32 -DNOERRNO -DNON_UNIX_STDIO -DAUTO -DCVODE_YES  -DHAVEDLL -DMYSTR1=$(MAJORVER) -DMYSTR2=$(MINORVER)  -I/opt/X11/include 
#LDFLAGS=  -m64 -L/usr/lib -L/usr/lib64
LDFLAGS=  -m32 -L/usr/X11R6/lib 
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
HEADERS = abort.h comline.h extra.h kinescope.h newhome.h stiff.h\
	  adj2.h cvband.h f2c.h lio.h newpars.h storage.h\
	  aniparse_avi.h cvdense.h fftn.h llnlmath.h nullcline.h struct.h\
	  aniparse.h cvdiag.h fio.h llnltyps.h numerics.h sue.h\
	  arrayplot.h cvode.h flags.h load_eqn.h odesol2.h tabular.h\
	  array_print.h cvspgmr.h fmt.h local.h toons.h\
	  autevd.h dae_fun.h form_ode.h lread.h open.h torus.h\
	  autlib1.h delay_handle.h fp.h lunch-new.h parser.h txtread.h\
	  autlib2.h del_stab.h gear.h macdirent.h parserslow.h uio.h\
	  autlim.h dense.h getvar.h macsysdirent.h phsplan.h userbut.h\
	  auto_define.h derived.h ggets.h main.h pop_list.h util.h\
	  auto.h dfe.h graf_par.h many_pops.h pp_shoot.h vector.h\
	  auto_nox.h diagram.h graphics.h markov.h rdfmt.h volterra2.h\
	  auto_x11.h dialog_box.h help_defs.h menudrive.h read_dir.h volterra.h\
	  axes2.h do_fit.h histogram.h menu.h rubber.h wref.h\
	  band.h dormpri.h homsup.h menus.h run_auto.h wrtfmt.h\
	  browse.h due.h iio.h myfonts.h scrngif.h wsfe.h\
	  calc.h edit_rhs.h init_conds.h mykeydef.h sfe.h wsne.h\
	  choice_box.h eig_list.h integrate.h my_pars.h shoot.h xpplim.h\
	  close.h endfile.h iterativ.h my_ps.h my_svg.h simplenet.h\
	  color.h err.h kbs.h my_rhs.h spgmr.h cv2.h tutor.h

BITMAPS = aniwin.bitmap alert.bitmap bc.bitmap eqns.bitmap info.bitmap ic.bitmap pp.bitmap\
	  array.bitmap browse.bitmap equilib.bitmap lines.bitmap txtview.bitmap\
	  auto.bitmap delay.bitmap graph.bitmap param.bitmap filebrowse.bitmap \
	  lineup.bitmap linedn.bitmap pageup.bitmap pagedn.bitmap home.bitmap start.bitmap  

SOURCES = main.c ggets.c menu.c rubber.c derived.c init_condold.c \
	  many_pops.c pop_list.c graphics.c dialog_box.c \
	  numerics.c choice_box.c color.c init_conds.c \
	  browse.c kinescope.c  axes2.c abort.c \
	   parser2.c storage.c load_eqn.c lunch-new.c \
	  form_ode.c odesol2.c gear.c eig_list.c \
	  integrate.c delay_handle.c graf_par.c\
	  my_ps.c my_svg.c nullcline.c torus.c pp_shoot.c\
	  lunch.c calc.c adj2.c  my_rhs.c dormpri.c\
	  volterra2.c tabular.c markov.c histogram.c \
	  autlib1.c autlib2.c autlib3.c autevd.c run_auto.c autpp.c \
	  diagram.c auto.c flowkm.c comline.c edit_rhs.c do_fit.c \
	  flags.c del_stab.c stiff.c arrayplot.c array_print.c \
	  aniparse.c simplenet.c dae_fun.c read_dir.c  parserslow2.c \
	  kinescope_old.c fftn.c  extra.c funexample.c scrngif.c \
	  kinescope_avi.c aniparse_avi.c  nagroutines.c flowkm_small.c \
	  homsup.c txtread.c menudrive.c rtsafe.c vector.c userbut.c \
	  lbf_drive.c auto_nox.c auto_x11.c cli.c rdfmt_old.c\
	  Version.c backspace.c dfe.c due.c iio.c inquire.c rewind.c rsfe.c \
	  rdfmt.c sue.c uio.c wsfe.c sfe.c fmt.c lio.c lread.c open.c \
	  close.c util.c endfile.c wrtfmt.c wref.c err.c fmtlib.c rsne.c wsne.c \
	  fc77.c f2cstart.c \
	  cvode.c  cvdense.c  dense.c  cvband.c \
	  band.c  cvdiag.c  cvspgmr.c  spgmr.c  \
	  iterativ.c  vector.c  llnlmath.c cv2.c \
	  mainold.c load_eqnold.c ggetsold.c comlineold.c 
	  

          
#
# WARNING: For 64 bit machines replace parser2.o with parserslow2.o 
#          in the OBJECTS section!
#
OBJECTS = main.o ggets.o menu.o  rubber.o derived.o\
	many_pops.o  pop_list.o  graphics.o dialog_box.o \
	numerics.o choice_box.o color.o init_conds.o \
        browse.o kinescope.o axes2.o abort.o \
        parserslow2.o storage.o load_eqn.o\
	form_ode.o odesol2.o gear.o eig_list.o\
        integrate.o delay_handle.o graf_par.o dormpri.o\
	my_ps.o my_svg.o nullcline.o torus.o pp_shoot.o \
	lunch-new.o calc.o adj2.o  my_rhs.o read_dir.o\
        volterra2.o tabular.o markov.o histogram.o \
	comline.o edit_rhs.o do_fit.o flags.o del_stab.o stiff.o \
        arrayplot.o array_print.o aniparse.o simplenet.o dae_fun.o \
        fftn.o extra.o scrngif.o nagroutines.o homsup.o txtread.o \
        menudrive.o userbut.o \
	Version.o backspace.o dfe.o due.o iio.o inquire.o rewind.o rsfe.o \
	rdfmt.o sue.o uio.o wsfe.o sfe.o fmt.o lio.o lread.o open.o \
	close.o util.o endfile.o wrtfmt.o wref.o err.o fmtlib.o rsne.o wsne.o \
        fc77.o f2cstart.o \
	cvode.o  cvdense.o  dense.o  cvband.o \
        band.o  cvdiag.o  cvspgmr.o  spgmr.o  \
        iterativ.o  vector.o  llnlmath.o cv2.o
	
#If you are using sundials cvode libraries then you may comment out following above
#	cvode.o  cvdense.o  dense.o  cvband.o \
#        band.o  cvdiag.o  cvspgmr.o  spgmr.o  \
#        iterativ.o  vector.o  llnlmath.o cv2.o
 
LIB_OBJECTS = main.o ggets.o menu.o  rubber.o derived.o\
	many_pops.o  pop_list.o  graphics.o dialog_box.o \
	numerics.o choice_box.o color.o init_conds.o \
        browse.o kinescope.o axes2.o abort.o \
        parserslow2.o storage.o load_eqn.o\
	form_ode.o odesol2.o gear.o eig_list.o\
        integrate.o delay_handle.o graf_par.o dormpri.o\
	my_ps.o my_svg.o nullcline.o torus.o pp_shoot.o \
	lunch-new.o calc.o adj2.o read_dir.o\
        volterra2.o tabular.o markov.o histogram.o \
	comline.o edit_rhs.o do_fit.o flags.o del_stab.o stiff.o \
        arrayplot.o array_print.o aniparse.o simplenet.o dae_fun.o \
         fftn.o extra.o scrngif.o nagroutines.o homsup.o txtread.o \
        menudrive.o userbut.o 
AUTOOBJ = autlib1.o autlib2.o autlib3.o autevd.o run_auto.o autpp.o \
	diagram.o auto_nox.o auto_x11.o flowkm_small.o 
           
######################################################################
#
#
xppaut: $(OBJECTS) $(AUTOOBJ)
#
	#cd libI77 ; make
	#cd cvodesrc ; make
###########################################################
# Okay, here we go
#############################################
	#$(CC) -DAUTO -o xppaut $(OBJECTS) $(AUTOOBJ) $(LDFLAGS) $(LIBS) ./libf2cm.a ./libcvode.a
	#$(CC) -DAUTO -o xppaut $(OBJECTS) $(AUTOOBJ) $(LDFLAGS) $(LIBS)  ./libcvode.a
########################################################
########################################################
# Probably never need this but here it is
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
	 
	 $(CC) -DAUTO -o xppaut $(OBJECTS) $(AUTOOBJ) $(LDFLAGS) $(LIBS) $(AUTLIBS)

#   
# 
###################################################################
################################################################
mkmyfun:
	make -f Makefile.lib
mkI77:
	cd libI77 ; make
mkcvode:
	cd cvodesrc ; make
mkmkavi:
	cd mkavi ; make
#
# Edited for Debian GNU/Linux.
install: xppaut 
# Make necessary installation directories
	mkdir -p $(DESTDIR)$(BINDIR)
	mkdir -p $(DESTDIR)$(DOCDIR)/html
	mkdir -p $(DESTDIR)$(DOCDIR)/examples
	mkdir -p $(DESTDIR)$(DOCDIR)/xbm
	#mkdir -p /usr/X11R6/man/man1
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
xpplib: $(LIB_OBJECTS) $(AUTOOBJ)
	(ar rcv libxpp.a $(LIB_OBJECTS) $(AUTOOBJ) cvodesrc/*.o  ; ranlib libxpp.a)
	#(ar rcv libxpp.a $(LIB_OBJECTS) $(AUTOOBJ) cvodesrc/*.o libI77/*.o ; ranlib libxpp.a)
####################################################
#  tar file
####################################################
tarfile:
	tar cvf xppaut$(VERSION).tar $(SOURCES) $(AUTOSRC) $(HEADERS) $(BITMAPS) default.opt \
	 xpp_doc.tex README Makefile Makefile.64 Makefile.nice  Makefile.lib Makefile.avi Makefile.old Makefile.w7b \
	ode/*.*  canonical/*.* xpp_doc.pdf xpp_sum.tex xpp_sum.pdf  xpp_doc.tex tree.tex tree.pdf nullcline_bw.c  \
	xppaut.1\
        mkavi/*.cc mkavi/*.h mkavi/Makefile mkavi/drive.c help/*.html \
	help/odes/*.ode help/odes/*.c install.pdf install.tex LICENSE HISTORY \
	Makefile.s2x sbml2xpp.c *.xbm *.jar changes.txt \
xppaut-stylesheet.css	xppaut-stylesheet2.css animsvgwww
	gzip xppaut$(VERSION).tar 
otarfile:
	tar cvf xppaut$(VERSION).tar $(SOURCES) $(AUTOSRC) $(HEADERS) $(BITMAPS) default.opt \
	 xpp_doc.tex README Makefile Makefile.nice  Makefile.lib Makefile.avi Makefile.old \
	ode/*.* xpp_doc.ps xpp_doc.pdf xpp_sum.tex xpp_sum.pdf xpp_sum.ps nullcline_bw.c  \
        libI77/*.c libI77/*.h libI77/Makefile \
	cvodesrc/*.c cvodesrc/*.h cvodesrc/Makefile xppaut.1\
        mkavi/*.cc mkavi/*.h mkavi/Makefile mkavi/drive.c help/*.html \
	help/odes/*.ode help/odes/*.c install.pdf install.tex LICENSE HISTORY \
	Makefile.s2x sbml2xpp.c *.xbm *.jar
	gzip xppaut$(VERSION).tar 
##############################################
#  pack up a binary
##############################################
binary:
	strip xppaut;tar zvcf binary.tgz xppaut *.pdf *.xbm xppaut.1 $(ODES) $(DOC) $(HELP) README HISTORY LICENSE
##############################################
#  clean
##############################################
clean:
	#rm -f *.o *.a  libI77/*.o libI77/*.a cvodesrc/*.o cvodesrc/*.a xppaut
	#rm -f *.o *.a   cvodesrc/*.o cvodesrc/*.a xppaut
	rm -f *.o *.a xppaut 
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
         
