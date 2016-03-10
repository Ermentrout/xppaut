AM_CFLAGS = -std=c99 -Wall -pedantic -Ilibf2c
AM_CPPFLAGS = \
	-D_XOPEN_SOURCE=600

lib_LIBRARIES = \
	libauto.a

noinst_HEADERS = \
	fcon.h \
	auto.h \
	auto_c.h \
	auto_def2.h \
	auto_f2c.h \
	auto_mpi.h \
	auto_types.h

libauto_a_SOURCES = \
	autlib1.c \
	autlib2.c \
	autlib3.c \
	autlib4.c \
	autlib5.c \
	conpar.c \
	eispack.c \
	setubv.c \
	worker.c \
	$(libf2c_SOURCES)

libf2c_SOURCES = \
	libf2c/cabs.c \
	libf2c/d_lg10.c \
	libf2c/d_imag.c \
	libf2c/d_sign.c \
	libf2c/i_dnnt.c \
	libf2c/i_nint.c \
	libf2c/pow_dd.c \
	libf2c/pow_di.c \
	libf2c/pow_ii.c \
	libf2c/r_lg10.c \
	libf2c/z_abs.c \
	libf2c/z_exp.c \
	libf2c/z_log.c
