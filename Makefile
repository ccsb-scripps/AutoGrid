#
# Makefile to build AutoGrid from Object files.
#
# (c) 1994-2001, TSRI
# Garrett M. Morris
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
       get_atom_type.o \
	   parm-ii.o

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
      get_atom_type.ln \
      strindex.ln \
      parm-ii.ln \


CC = cc # SGI, Sun, MacOS X
# CC = cxx # Alpha.
# CC = gcc # use this if you have the Gnu compiler, as on Linux, MkLinux, LinuxPPC systems.

LIB = -lm # SGI, Sun, Linux, MacOS X
# LIB = -lm -lc # Alpha, Convex.
# LIB = -lm -lg++ # HP, Gnu.

CSTD = $(DBUG) $(PROF) $(WARN) # SGI, Sun, Linux, MacOS X
# CSTD = $(DBUG) $(PROF) $(WARN) -std # Convex.
# CSTD = -std -verbose $(PROF) $(DBUG) $(WARN) # Alpha. Not sarah
# CSTD = -std arm -verbose $(PROF) $(DBUG) $(WARN) # Alpha. sarah
# CSTD = -DHPPA -D_HPUX_SOURCE -ansi $(PROF) $(DBUG) $(WARN) # HP.

CFLAGS = $(OPT) # SGI, HP, Alpha, Sun, Convex, MacOS X: Optimize the object files, too.

CFLAGS = $(PROF) $(DBUG) # SGI, Linux, MacOS X
# CFLAGS = -std -verbose $(PROF) $(DBUG) # DEC Alpha, Convex
# CFLAGS = -Aa -D_HPUX_SOURCE $(PROF) $(DBUG) # Hewlett-Packard

OLIMIT = # SGI, Convex, Linux, MacOS X, -g debugging
# OLIMIT = -Olimit 2048 # Sun, Hewlett-Packard, DEC Alpha

OPTLEVEL = -O3 # Agressive optimization.
# OPTLEVEL = -O2 # High optimization.
# OPTLEVEL = -O1 # Do optimizations that can be done quickly;default.
# OPTLEVEL = -O0 # Do not optimize.

# OPT_SGI_IPNUM = # Alpha, HP, Sun, Convex, SGI, Linux, MacOS X.
# OPT_SGI_IPNUM = -Ofast=ip19 # SGI, 'uname -a' says 'IP19'
# OPT_SGI_IPNUM = -Ofast=ip21 # SGI, 'uname -a' says 'IP21'
# OPT_SGI_IPNUM = -Ofast=ip25 # SGI, 'uname -a' says 'IP25' PowerChallenge is R10000, IP25
OPT_SGI_IPNUM = -Ofast=ip27 # SGI, 'uname -a' says 'IP27'  Atlas Cluster is R12000, IP27
# OPT_SGI_IPNUM = -Ofast=ip30 # SGI, 'uname -a' says 'IP30'

# OPT_SGI_R000 = # Alpha, HP, Sun, Convex, SGI, Linux, MacOS X.
# OPT_SGI_R000 = -r4000 -mips2 # SGI, 'hinv' says MIPS Processor is R4000
# OPT_SGI_R000 = -r8000 -mips4 # SGI, 'hinv' says MIPS Processor is R8000
# OPT_SGI_R000 = -r10000 -mips4 # SGI, 'hinv' says MIPS Processor is R10000
# OPT_SGI_R000 = -r12000 -mips4 # SGI, 'hinv' says MIPS Processor is R12000
OPT_SGI_R000 = -r14000 -mips4 # SGI, 'hinv' says MIPS Processor is R14000

# OPT_ARCH_SPECIFIC = # Alpha, HP, Sun, Convex, SGI, Linux, MacOS X.
# OPT_ARCH_SPECIFIC = -n32 $(OPT_SGI_IPNUM) $(OPT_SGI_R000) -IPA $(LNO_OPT) # SGI, new 32-bit
# OPT_ARCH_SPECIFIC = -n32 $(OPT_SGI_IPNUM) $(OPT_SGI_R000) -IPA $(LNO_OPT) -DUSE_INT_AS_LONG # SGI (long is 8bytes).
# OPT_ARCH_SPECIFIC = $(OPT_SGI_IPNUM) $(OPT_SGI_R000) $(LNO_OPT) # SGI, not new 32-bit
# OPT_ARCH_SPECIFIC = -ifo # DEC Alpha
# OPT_ARCH_SPECIFIC = -O +O3 +Obb2048 # Hewlett-Packard
# OPT_ARCH_SPECIFIC = -ur # Convex
# OPT_ARCH_SPECIFIC = -g # debugging
OPT_ARCH_SPECIFIC = # Mac OS X

OPT = $(CSTD) $(OPTLEVEL) $(OPT_ARCH_SPECIFIC) # All platforms

LINKOPT = $(OPT) # SGI, HP, Alpha, Sun, Convex, Linux.
# LINKOPT = $(CSTD) -O2 -r4000 -IPA $(LNO_OPT) # SGI/IRIX5: R4000.

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


autogrid3 : $(OBJS)
	@echo $(EXE) on `date` using' >>> ' `hostname` ' <<<' >> LATEST_MAKE
	@echo 'Flags: '$(CC) $(LINK) $(LIB) >> LATEST_MAKE
	$(CC) $(LINK) -o $@ $(OBJS) $(LIB)

convertmap : convertmap.o
	$(CC) $(OPT) -o $@ convertmap.o $(LIB)


main.o : main.c autogrid.h autoglobal.h parm-ii.c parm-ii.h banner.c gpfparser.c gpftoken.h
	$(CC) $(OPT) $(OLIMIT) -c main.c

check_size.o : check_size.c autogrid.h
	$(CC) $(OPT) -c check_size.c

setflags.o : setflags.c autogrid.h
	$(CC) $(OPT) -c setflags.c

timesyshms.o : timesyshms.c autogrid.h
	$(CC) $(OPT) -c timesyshms.c

timesys.o : timesys.c autogrid.h
	$(CC) $(OPT) -c timesys.c

printhms.o : printhms.c autogrid.h
	$(CC) $(OPT) -c printhms.c

prHMSfixed.o : prHMSfixed.c autogrid.h
	$(CC) $(OPT) -c prHMSfixed.c

printdate.o : printdate.c autogrid.h
	$(CC) $(OPT) -c printdate.c

banner.o : banner.c autogrid.h
	$(CC) $(OPT) -c banner.c

gpfparser.o : gpfparser.c autogrid.h gpftoken.h
	$(CC) $(OPT) -c gpfparser.c

get_atom_type.o : get_atom_type.c autogrid.h gpftoken.h
	$(CC) $(OPT) -c get_atom_type.c

strindex.o : strindex.c autogrid.h
	$(CC) $(OPT) -c strindex.c

convertmap.o : convertmap.c
	$(CC) $(OPT) -c convertmap.c

parm-ii.o : parm-ii.c parm-ii.h
	$(CC) $(OPT) -c parm-ii.c

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

get_atom_type.ln : get_atom_type.c autogrid.h gpftoken.h
	$(LINT) $(LINTFLAGS) $?

strindex.ln : strindex.c autogrid.h
	$(LINT) $(LINTFLAGS) $?

parm-ii.ln : parm-ii.c parm-ii.h
	$(LINT) $(LINTFLAGS) $?

#
# Remove objects, lint files, cores, etc.
#
clean:
	/bin/rm -rf *.o *.s *.ln a.out core autogrid3 convertmap

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
