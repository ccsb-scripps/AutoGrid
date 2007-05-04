## Process this file with automake to produce Makefile.in

# $Id: Makefile.am,v 1.3 2007/05/04 01:38:39 garrett Exp $
# 
# AutoGrid 
# 
# Copyright (C) 1989-2007,  Garrett M. Morris,   TSRI
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

EXTRA_DIST = AD4_parameters.dat paramdat2h.csh Tests

dist-hook:
	rm -rf `find $(distdir)/Tests -name CVS -o -name '*.glg'` 

bin_PROGRAMS = autogrid4

# Define the AutoGrid 4 source code files:

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

AM_CFLAGS = -Wall -O3 # All warnings, Agressive optimization

AM_CXXFLAGS = -Wall -O3 # All warnings, Agressive optimization

AM_CPPFLAGS = -I$(srcdir)/../autodock-4.0.0

autogrid4_CFLAGS = -Wall -O3 # All warnings, Agressive optimization
# autogrid4_CFLAGS = -Wall -O3 -ftree-vectorize # All warnings, Agressive optimization, GCC vectorization in Mac OS X

autogrid4_CXXFLAGS = -Wall -O3 # All warnings, Agressive optimization
# autogrid4_CXXFLAGS = -Wall -O3 -ftree-vectorize # All warnings, Agressive optimization, GCC vectorization in Mac OS X


# LINKOPT = $(CSTD) $(OPT) -Wl,--stack=0x2000000 # Cygwin, 32MB stacksize

# DBUG = -DNDEBUG # No debugging and no assert code.
# DBUG = # Use assert code.
# DBUG = -g # dbx.
# DBUG = -g -DDEBUG # dbx + DEBUG-specific code.
# DBUG = -g3 # dbx + optimization.
# DBUG = -g3 -DDEBUG # dbx + optimization, + DEBUG-specific code.
# DBUG = -DDEBUG # Just DEBUG-specific code.
# DBUG = -DDEBUG2 # Just DEBUG2-specific code for tracking prop.selection.
# DBUG = -DDEBUG3 # Just DEBUG3-specific code for print age of individuals.

check : autogrid4 $(srcdir)/Tests/test_autogrid4.py
#if HAVE_PYTHON
	@echo " "
	@echo "Running the Python test_autogrid4.py unit-tests on autogrid4"
	@echo " "
	cd $(srcdir)/Tests ; chmod u+w . ; python ./test_autogrid4.py
#else
#	@echo " "
#	@echo "Python version 2.1 or higher needed for 'make check', but not found."
#	@echo "Skipping AutoGrid4 unit tests."
#	@echo " "
#endif


main.o : $(srcdir)/../autodock-4.0.0/autocomm.h

# 2007-04-30 MP & GMM: this line is necessary to be able to compile read_parameter_library.o
read_parameter_library.o : default_parameters.h

default_parameters.h : AD4_parameters.dat paramdat2h.csh
	rm -f $@
	csh ./paramdat2h.csh > tmp-paramdat
	mv -f tmp-paramdat $@

# see Autoconf manual 4.7.3 (p. 25) support for VPATH
# enabling multi-platform builds
srcdir = @srcdir@
VPATH  = @srcdir@

#
# EOF
#
