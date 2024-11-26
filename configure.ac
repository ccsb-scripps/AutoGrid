#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
# 4.2.5 November 2 2012
# 4.2.5.1 December 2012
# 4.2.6 July 2014
#AC_INIT([autogrid], [4.2.6], [autodock@scripps.edu])
# the m4_esyscmd following creates a string such as VERSION='4.2.5.x.2014-12-25'
# indicating when the "autoreconf" was most recently run.
# This date is also inserted into PACKAGE_VERSION and PACKAGE_STRING
# AC_INIT([autogrid],[m4_translit(m4_esyscmd(date +4.2.7.x.%Y-%m-%d),m4_newline)],[autodock@scripps.edu])
AC_INIT([autogrid], [4.2.8], [autodock@scripps.edu])
AC_CONFIG_SRCDIR([atom_parameter_manager.cpp])
#AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE(subdir-objects)

# Checks for programs.
AC_PROG_CXX

AC_CANONICAL_HOST

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([limits.h stddef.h stdlib.h string.h sys/param.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
#AC_CHECK_FUNCS([gethostname pow round sqrt strncasecmp])
AC_CHECK_FUNCS([gethostname getpid GetProcessId fchmod setlinebuf sysconf times getrusage pow sqrt strncasecmp])



case $host_os in
  *cygwin* ) CYGWIN=yes;;
         * ) CYGWIN=no;;
esac
AM_CONDITIONAL([WANT_BIGSTACK], [test "$CYGWIN" = yes]) 

# AC_CONFIG_LINKS(parse_param_line.cpp:../autodock/parse_param_line.cc)
# AC_CONFIG_LINKS(read_parameter_library.cpp:../autodock/read_parameter_library.cc)
# AC_CONFIG_LINKS(stop.cpp:../autodock/stop.cc)

# M Pique: Open MP detection taken from Joachim Wuttke at
# http://stackoverflow.com/questions/9037434/enable-openmp-checks-and-flags-in-autoconf-and-automake
# who suggests also:   AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")
AC_LANG(C++)
AC_PROG_CXX
AC_OPENMP  # requires autoconf 2.62 or newer, sets OPENMP_CXXFLAGS if successful

AM_CONDITIONAL([WANT_OMP], [test "x$OPENMP_CXXFLAGS" != "x"])
# See "automake" manual, version 1.10, section 10.5 (Python)
#AM_PATH_PYTHON([2.1],, [:])
#AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
