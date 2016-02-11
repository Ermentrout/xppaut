#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([xppaut], [7.0.0], [xppaut-l_at_mailman.ucsd.edu])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([foreign])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([X11], [XOpenDisplay])
AC_CHECK_LIB([dl], [dlopen])
AC_CHECK_LIB([m], [floor])

# Checks for header files.
AC_PATH_X
AC_CHECK_HEADERS([ \
    fcntl.h \
    limits.h \
    memory.h \
    stdlib.h \
    string.h \
    strings.h \
    sys/time.h \
    unistd.h \
    wctype.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([ \
    atexit \
    floor \
    getcwd \
    gettimeofday \
    memmove \
    memset \
    pow \
    sqrt \
    strcasecmp \
    strcspn \
    strdup \
    strncasecmp])

MAJORVER=${PACKAGE_VERSION%.*}
AC_SUBST([MAJORVER])
MINORVER=${PACKAGE_VERSION##*.}
AC_SUBST([MINORVER])

AC_CONFIG_FILES([
    Makefile
    betty/Makefile
    bin/Makefile
    canonical/Makefile
    doc/Makefile
    help/Makefile
    mkavi/Makefile
    ode/Makefile
    share/Makefile
    src/Makefile])
AC_OUTPUT
