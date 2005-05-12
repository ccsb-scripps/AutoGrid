#
# Makefile to build AutoGrid from Object files.
#
# (c) 1994-2004, TSRI
# Garrett M. Morris, Ruth Huey, David S. Goodsell
#

EXE = .

OBJS = main.o \
       check_size.o \
       setflags.o \
       timesys.o \
       timesyshms.o \
       printhms.o \
       prHMSfixed.o \
       printdate.o \
       strindex.o \
       banner.o \
       gpfparser.o \
       parsetypes.o \
       atom_parameter_manager.o 

LNS = \
      main.ln \
      check_size.ln \
      setflags.ln \
      timesyshms.ln \
      timesys.ln \
      printhms.ln \
      prHMSfixed.ln \
      printdate.ln \
      banner.ln \
      gpfparser.ln \
      parsetypes.ln \
      strindex.ln \
      atom_parameter_manager.ln


CC = g++
# CC = cc # SGI, Sun, MacOS X
# CC = cxx # Alpha.
# CC = gcc # use this if you have the Gnu compiler, as on Linux, MkLinux, LinuxPPC systems.

LIB = -lm # SGI, Sun, Linux, MacOS X
# LIB = -lm -lc # Alpha, Convex.
# LIB = -lm -lg++ # HP, Gnu.

CFLAGS = $(DBUG) $(OPT) -I../autodock $(WARN) $(PROF) # SGI, HP, Alpha, Sun, Convex, Linux, MacOS X

OLIMIT = # SGI, Convex, Linux, MacOS X, -g debugging
# OLIMIT = -Olimit 2048 # Sun, Hewlett-Packard, DEC Alpha

OPTLEVEL = -O3 # Agressive optimization.
# OPTLEVEL = -O2 # High optimization.
# OPTLEVEL = -O1 # Do optimizations that can be done quickly;default.
# OPTLEVEL = -O0 # Do not optimize.

OPT_SGI_IPNUM = # Alpha, HP, Sun, Convex, SGI, Linux, MacOS X.
# OPT_SGI_IPNUM = -Ofast=ip19 # SGI, 'uname -a' says 'IP19'
# OPT_SGI_IPNUM = -Ofast=ip21 # SGI, 'uname -a' says 'IP21'
# OPT_SGI_IPNUM = -Ofast=ip25 # SGI, 'uname -a' says 'IP25' PowerChallenge is R10000, IP25
# OPT_SGI_IPNUM = -Ofast=ip27 # SGI, 'uname -a' says 'IP27'  Atlas Cluster is R12000, IP27
# OPT_SGI_IPNUM = -Ofast=ip30 # SGI, 'uname -a' says 'IP30'

OPT_SGI_R000 = # Alpha, HP, Sun, Convex, SGI, Linux, MacOS X.
# OPT_SGI_R000 = -r4000 -mips2 # SGI, 'hinv' says MIPS Processor is R4000
# OPT_SGI_R000 = -r8000 -mips4 # SGI, 'hinv' says MIPS Processor is R8000
# OPT_SGI_R000 = -r10000 -mips4 # SGI, 'hinv' says MIPS Processor is R10000
# OPT_SGI_R000 = -r12000 -mips4 # SGI, 'hinv' says MIPS Processor is R12000
# OPT_SGI_R000 = -r14000 -mips4 # SGI, 'hinv' says MIPS Processor is R14000

OPT_ARCH_SPECIFIC = # Alpha, HP, Sun, Convex, SGI, Linux, MacOS X.
# OPT_ARCH_SPECIFIC = -n32 $(OPT_SGI_IPNUM) $(OPT_SGI_R000) -IPA $(LNO_OPT) # SGI, new 32-bit
# OPT_ARCH_SPECIFIC = -n32 $(OPT_SGI_IPNUM) $(OPT_SGI_R000) -IPA $(LNO_OPT) -DUSE_INT_AS_LONG # SGI (long is 8bytes).
# OPT_ARCH_SPECIFIC = $(OPT_SGI_IPNUM) $(OPT_SGI_R000) $(LNO_OPT) # SGI, not new 32-bit
# OPT_ARCH_SPECIFIC = -ifo # DEC Alpha
# OPT_ARCH_SPECIFIC = -O +O3 +Obb2048 # Hewlett-Packard
# OPT_ARCH_SPECIFIC = -ur # Convex
# OPT_ARCH_SPECIFIC = -g # debugging
# OPT_ARCH_SPECIFIC = # Mac OS X

OPT = $(OPTLEVEL) $(OPT_ARCH_SPECIFIC) # All platforms

LINKOPT = $(OPT) # SGI, HP, Alpha, Sun, Convex, Linux.
# LINKOPT = -O2 -r4000 -IPA $(LNO_OPT) # SGI/IRIX5: R4000.

LINK = $(LINKOPT) # Linking flags.
# LINK = $(LINKOPT) -cord # Procedure rearranger on SGI.

LNO_OPT = # SGI, no special optimization at link time
# LNO_OPT = -LNO:auto_dist=ON:gather_scatter=2 # SGI

LINT    = lint # C code checking

LINTFLAGS = $(LIB) -lansi -c # DEC-Alpha # MA
# LINTFLAGS = -DHPPA -D_HPUX_SOURCE -Aa $(LIB) -c # HP/Hewlett-Packard
# LINTFLAGS = -Aa $(LIB) -c # HP/Hewlett-Packard Alt
# LINTFLAGS = $(LIB) -n -lansic # Sun Sparc


DBUG = -DNDEBUG # No debugging and no assert code.
# DBUG = # Use assert code.
# DBUG = -g # dbx.
# DBUG = -g -DDEBUG # dbx + DEBUG-specific code.
# DBUG = -g3 # dbx + optimization.
# DBUG = -g3 -DDEBUG # dbx + optimization, + DEBUG-specific code.
# DBUG = -DDEBUG # Just DEBUG-specific code.
# DBUG = -DDEBUG2 # Just DEBUG2-specific code for tracking prop.selection.
# DBUG = -DDEBUG3 # Just DEBUG3-specific code for print age of individuals.

PROF = # No profiling.
# PROF = -p # Profiling.

WARN = # Default warning level.
# WARN = -woff all # For no warnings.
# WARN = -fullwarn -ansiE -ansiW # For full warnings during compilation.


autogrid4 : $(OBJS)
	@echo $(EXE) on `date` using' >>> ' `hostname` ' <<<' >> LATEST_MAKE
	@echo 'Flags: '$(CC) $(LINK) $(LIB) >> LATEST_MAKE
	$(CC) $(LINK) -o $@ $(OBJS) $(LIB)

convertmap : convertmap.o
	$(CC) $(OPT) -o $@ convertmap.o $(LIB)


main.o : main.c autogrid.h autoglobal.h  banner.c gpfparser.c gpftoken.h parsetypes.c

check_size.o : check_size.c autogrid.h

setflags.o : setflags.c autogrid.h

timesyshms.o : timesyshms.c autogrid.h

timesys.o : timesys.c autogrid.h

printhms.o : printhms.c autogrid.h

prHMSfixed.o : prHMSfixed.c autogrid.h

printdate.o : printdate.c autogrid.h

banner.o : banner.c autogrid.h

gpfparser.o : gpfparser.c autogrid.h gpftoken.h

parsetypes.o : parsetypes.c autogrid.h

strindex.o : strindex.c autogrid.h

convertmap.o : convertmap.c

atom_parameter_manager.o : atom_parameter_manager.c parameters.h


#
# lcheck dependencies...
#

main.ln : main.c autogrid.h autoglobal.h
	$(LINT) $(LINTFLAGS) $?

check_size.ln : check_size.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

setflags.ln : setflags.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

timesyshms.ln : timesyshms.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

timesys.ln : timesys.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

printhms.ln : printhms.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

prHMSfixed.ln : prHMSfixed.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

printdate.ln : printdate.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

banner.ln : banner.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

gpfparser.ln : gpfparser.c autogrid.h gpftoken.h
	$(LINT) $(LINTFLAGS) $?

parsetypes.ln : parsetypes.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

strindex.ln : strindex.c autogrid.h
	$(LINT) $(LINTFLAGS) $?


#
# Remove objects, lint files, cores, etc.
#
clean:
	/bin/rm -rf *.o *.s *.ln a.out core autogrid4 convertmap

#
# DEC Alpha
#
openalpha :
	/bin/cp -f obj.alpha/*.o .
	/bin/cp -f obj.alpha/*.ln .

closealpha :
	/bin/mv -f *.o obj.alpha
	/bin/mv -f *.ln obj.alpha
#
# Convex
#
openc2 :
	/bin/cp -f obj.c2/*.o .
	/bin/cp -f obj.c2/*.ln .

closec2 :
	/bin/mv -f *.o obj.c2
	/bin/mv -f *.ln obj.c2

#
# Hewlett-Packard Precision Architecture
#
openhppa :
	/bin/cp -f obj.hppa/*.o .
	/bin/cp -f obj.hppa/*.ln .

closehppa :
	/bin/mv -f *.o obj.hppa
	/bin/mv -f *.ln obj.hppa

#
# Silicon Graphics
#
opensgi4D :
	/bin/cp -f obj.sgi4D/*.o .
	/bin/cp -f obj.sgi4D/*.ln .

closesgi4D :
	/bin/mv -f *.o obj.sgi4D
	/bin/mv -f *.ln obj.sgi4D

#
# Sun
#
opensun4 :
	/bin/cp -f obj.sun4/*.o .
	/bin/cp -f obj.sun4/*.ln .

closesun4 :
	/bin/mv -f *.o obj.sun4
	/bin/mv -f *.ln obj.sun4

#
# Architecture
#
openxxx :
	/bin/cp -f obj.xxx/*.o .
	/bin/cp -f obj.xxx/*.ln .

closexxx :
	/bin/mv -f *.o obj.xxx
	/bin/mv -f *.ln obj.xxx
