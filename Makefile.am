## Process this file with automake to produce Makefile.in

# $Id: Makefile.am,v 1.21 2009/09/09 23:02:11 rhuey Exp $
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

EXTRA_DIST = Tests

dist-hook:
	rm -rf `find $(distdir)/Tests -name CVS -o -name '*.glg'` 

bin_PROGRAMS = autogrid4 

autogrid4_SOURCES = main.cpp \
       check_size.cpp \
       setflags.cpp \
       timesys.cpp \
       timesyshms.cpp \
       printhms.cpp \
       prHMSfixed.cpp \
       printdate.cpp \
       strindex.cpp \
       banner.cpp \
       gpfparser.cpp \
       parsetypes.cpp \
       atom_parameter_manager.cpp \
       read_parameter_library.cpp \
       parse_param_line.cpp \
       distdepdiel.cpp \
            atom_parameter_manager.h \
            autoglobal.h \
            autogrid.h \
            constants.h \
            default_parameters.h \
            distdepdiel.h \
            gpftoken.h \
            grid.h \
            parameters.h \
            parm.h \
            parse_param_line.h \
            partokens.h \
            printdate.h \
            printhms.h \
            prototypes.h \
            ranlib.h \
            read_parameter_library.h \
            structs.h \
            times.h \
            timesys.h \
            timesyshms.h \
            typedefs.h \
            util.h \
            version.h


# Define the AutoGrid 4 source code files:
# See GNU Automake documentation, section 10.4 re BUILT_SOURCES;
# "default_parameters.h" is generated using the rules below
# which invoke "paramdat2h.csh" to convert the AD4.1_bound.dat
# parameter library file into a C++ header file.
BUILT_SOURCES = default_parameters.h

CLEANFILES = default_parameters.h

AM_CFLAGS = -Wall -O3 # All warnings, Agressive optimization

AM_CPPFLAGS = -Wall -O3 -I$(srcdir)/../autodock # All warnings, Agressive optimization
AM_CPPFLAGS += -Wno-write-strings # increase compatibility with pre 4.x compilers for string handling

# autogrid4_CXXFLAGS = -Wall -O3 -ftree-vectorize # All warnings, Agressive optimization, GCC vectorization in Mac OS X

autogrid4_CPPFLAGS=$(AM_CPPFLAGS)

# LINKOPT = $(CSTD) $(OPT) -Wl,--stack=0x2000000 # Cygwin, 32MB stacksize

if WANT_BIGSTACK
    autogrid4_LDFLAGS = -Wl,--stack,8388608 # Cygwin, 32MB stacksize
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
# for technical reasons, autogrid must be run in the Tests
# directory because it writes its maps into its current directory
# and expects to find AD4.1_bound.dat in that same directory.
	outputdir=`pwd` ;\
	umask 0 ;\
	test -d $$outputdir/Tests || mkdir $$outputdir/Tests ;\
	cp $(srcdir)/Tests/*.py   $$outputdir/Tests;\
	cp $(srcdir)/Tests/*.gpf  $$outputdir/Tests;\
	cp $(srcdir)/Tests/*.dat  $$outputdir/Tests;\
	cp $(srcdir)/Tests/*.pdbqt  $$outputdir/Tests;\
	cp $(srcdir)/Tests/*.map*  $$outputdir/Tests;\
	cd Tests ;  python ./test_autogrid4.py -d . -e ../autogrid4$(EXEEXT) -o .
#else
#	@echo " "
#	@echo "Python version 2.1 or higher needed for 'make test', but not found."
#	@echo "Skipping AutoGrid4 unit tests."
#	@echo " "
#endif


main.o : $(srcdir)/../autodock/autocomm.h

default_parameters.h : $(srcdir)/../autodock/paramdat2h.csh  $(srcdir)/../autodock/AD4_parameters.dat $(srcdir)/../autodock/AD4.1_bound.dat
	rm -f $@
	csh $(srcdir)/../autodock/paramdat2h.csh $(srcdir)/../autodock/AD4_parameters.dat  $(srcdir)/../autodock/AD4.1_bound.dat > $@

# 2007-04-30 MP & GMM: this line is necessary to be able to compile read_parameter_library.o
read_parameter_library.o : default_parameters.h

# see Autoconf manual 4.7.3 (p. 25) support for VPATH
# enabling multi-platform builds
srcdir = @srcdir@
VPATH  = @srcdir@

#
# EOF
#
