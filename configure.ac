#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([autogrid], [4.2.3], [autodock@scripps.edu])
AC_CONFIG_SRCDIR([atom_parameter_manager.cpp])
#AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE()

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
AC_CHECK_FUNCS([gethostname pow round sqrt strncasecmp])

case $host_os in
  *cygwin* ) CYGWIN=yes;;
         * ) CYGWIN=no;;
esac
AM_CONDITIONAL([WANT_BIGSTACK], [test "$CYGWIN" = yes]) 

AC_CONFIG_LINKS(read_parameter_library.cpp:../autodock/read_parameter_library.cc)

# See "automake" manual, version 1.10, section 10.5 (Python)
#AM_PATH_PYTHON([2.1],, [:])
#AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
