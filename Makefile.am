## Process this file with automake to produce Makefile.in

# $Id: Makefile.am,v 1.35 2020/05/28 21:33:47 mp Exp $
# 
# AutoGrid 
# 
# Copyright (C) 1989-2009, TSRI
# All Rights Reserved.
# 
# AutoGrid is a Trade Mark of The Scripps Research Institute.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# Makefile.am to build AutoGrid 4
#

# Autotools

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = Tests

dist-hook:
	rm -rf `find $(distdir)/Tests -name CVS -o -name '*.glg'` 

bin_PROGRAMS = autogrid4 

autogrid4_SOURCES = main.cpp \
       bondmanager.cpp \
       bhtree.cpp \
       check_size.cpp \
       setflags.cpp \
       ad4_shared/timesys.cc \
       ad4_shared/timesyshms.cc \
       ad4_shared/printhms.cc \
       prHMSfixed.cpp \
       ad4_shared/printdate.cc \
       strindex.cpp \
       banner.cpp \
       gpfparser.cpp \
       parsetypes.cpp \
       atom_parameter_manager.cpp \
       ad4_shared/read_parameter_library.cc \
       ad4_shared/parse_param_line.cc \
       distdepdiel.cpp \
       ad4_shared/memalloc.cc \
       ad4_shared/mingw_getrusage.cc \
       calc_vina_potential.cpp \
       ad4_shared/threadlog.cc \
       ad4_shared/targetfile.cc \
       ad4_shared/stop.cc \
       atom_parameter_manager.h \
       ad4_shared/autoglobal.h \
       autogrid.h \
       bhtree.h \
       bondmanager.h \
       ad4_shared/constants.h \
       default_parameters.h \
       distdepdiel.h \
       gpftoken.h \
       grid.h \
       ad4_shared/memalloc.h \
       ad4_shared/mingw_sysconf.h \
       ad4_shared/mingw_sys_times.h \
       ad4_shared/mystdlib.h \
       parameters.h \
       parm.h \
       parse_param_line.h \
       ad4_shared/partokens.h \
       ad4_shared/printdate.h \
       ad4_shared/printhms.h \
       prototypes.h \
       ranlib.h \
       ad4_shared/read_parameter_library.h \
       ad4_shared/targetfile.h \
       times.h \
       ad4_shared/timesys.h \
       ad4_shared/timesyshms.h \
       ad4_shared/typedefs.h \
       util.h \
       calc_vina_potential.h \
       ad4_shared/threadlog.h \
       gridmap.h \
       version.h 


# Define the AutoGrid 4 source code files:
# See GNU Automake documentation, section 10.4 re BUILT_SOURCES;
# "default_parameters.h" is generated using the rules below
# which invoke "paramdat2h.csh" to convert the AD4.1_bound.dat
# parameter library file into a C++ header file.
BUILT_SOURCES = default_parameters.h


AM_CFLAGS = -Wall -O3 # All warnings, Agressive optimization

AM_CPPFLAGS = -Wall -O3 -I $(srcdir) -I$(srcdir)/ad4_shared # All warnings, Agressive optimization
AM_CPPFLAGS += -Wno-write-strings # increase compatibility with pre 4.x compilers for string handling
#  M Pique - add USE_8A_NBCUTOFF
AM_CPPFLAGS += -DUSE_8A_NBCUTOFF
#  M Pique - add USE_DOUBLE to define Real as double instead of float, see typedefs.h
AM_CPPFLAGS += -DUSE_DOUBLE
#  MS 2015 Pique - add USE_BHTREE to speed up calculations
# AM_CPPFLAGS += -DUSE_BHTREE

# autogrid4_CXXFLAGS = -Wall -O3 -ftree-vectorize # All warnings, Agressive optimization, GCC vectorization in Mac OS X

if WANT_OMP  
# M Pique - I'm not sure which of these are really necessary
# could this be replaced by :  AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")
autogrid4_CFLAGS = $(AM_CFLAGS)  -fopenmp  $(OPENMP_CXXFLAGS)
autogrid4_CPPFLAGS = $(AM_CPPFLAGS)  $(OPENMP_CXXFLAGS)
autogrid4_LDFLAGS = $(AM_LDFLAGS) -fopenmp  $(OPENMP_CXXFLAGS)

autogrid4_CXXFLAGS = $(AM_CXXFLAGS)  -fopenmp  $(OPENMP_CXXFLAGS)
#autogrid4_LDADD = $(autogrid4_LDADD)
else
autogrid4_CFLAGS = $(AM_CFLAGS)  
autogrid4_CPPFLAGS = $(AM_CPPFLAGS)  
autogrid4_LDFLAGS = $(AM_LDFLAGS) 

autogrid4_CXXFLAGS = $(AM_CXXFLAGS)
endif

# LINKOPT = $(CSTD) $(OPT) -Wl,--stack=0x2000000 # Cygwin, 32MB stacksize

if WANT_BIGSTACK
    autogrid4_LDFLAGS += -Wl,--stack,8388608 # Cygwin, 32MB stacksize
endif


# DBUG = -DNDEBUG # No debugging and no assert code.
# DBUG = # Use assert code.
# DBUG = -g # dbx.
# DBUG = -g -DDEBUG # dbx + DEBUG-specific code.
# DBUG = -g3 # dbx + optimization.
# DBUG = -g3 -DDEBUG # dbx + optimization, + DEBUG-specific code.
# DBUG = -DDEBUG # Just DEBUG-specific code.
# DBUG = -DDEBUG2 # Just DEBUG2-specific code for tracking prop.selection.
# DBUG = -DDEBUG3 # Just DEBUG3-specific code for print age of individuals.

check : test

test : autogrid4 $(srcdir)/Tests/test_autogrid4.py
#if HAVE_PYTHON
	@echo " "
	@echo "Running the Python test_autogrid4.py unit-tests on autogrid4"
	@echo " "
# for technical reasons, autogrid must be run in the directory where the
# test files reside because it writes its maps into its current directory
# and expects to find AD4.1_bound.dat in that same directory.
# We create a run-time Testruns directory for this purpose.
# The name is different to allow "configure;make;make check" in the source dir.
	outputdir=`pwd` ;\
	umask 0 ;\
	test -d $$outputdir/Testruns && chmod -fR a+w $$outputdir/Testruns;\
	test -d $$outputdir/Testruns || mkdir $$outputdir/Testruns ;\
	cp $(srcdir)/Tests/AD4.1_bound.dat   $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.py   $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.gpf  $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.pdbqt  $$outputdir/Testruns;\
	cp $(srcdir)/Tests/*.map*  $$outputdir/Testruns;\
	chmod -fR a+w $$outputdir/Testruns;\
	cd $$outputdir/Testruns;\
	python ./test_autogrid4.py -d . -e ../autogrid4$(EXEEXT) -o .
#else
#	@echo " "
#	@echo "Python version 2.1 or higher needed for 'make test', but not found."
#	@echo "Skipping AutoGrid4 unit tests."
#	@echo " "
#endif


main.o : $(srcdir)/ad4_shared/autocomm.h

default_parameters.h : $(srcdir)/ad4_shared/paramdat2h.csh  $(srcdir)/ad4_shared/AD4_parameters.dat $(srcdir)/ad4_shared/AD4.1_bound.dat
	rm -f $@
	csh $(srcdir)/ad4_shared/paramdat2h.csh $(srcdir)/ad4_shared/AD4_parameters.dat  $(srcdir)/ad4_shared/AD4.1_bound.dat > $@

# 2007-04-30 MP & GMM: this line is necessary to be able to compile read_parameter_library.o
read_parameter_library.o : default_parameters.h

# see Autoconf manual 4.7.3 (p. 25) support for VPATH
# enabling multi-platform builds
srcdir = @srcdir@
VPATH  = @srcdir@

# see Autoconf manual for how to extend 'make clean' to directories:
clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -f default_parameters.h
	-rm -rf Testruns
#
# EOF
#
