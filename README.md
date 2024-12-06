# Autogrid4 

Autogrid4 is a support software for docking programs such as AutoDock4 and Autodock-GPU. 
Its function is to precalculate the grids used by the docking software. 

AutoDock is a suite of automated docking tools. It is designed to
predict how small molecules, such as substrates or drug candidates,
bind to a receptor of known 3D structure.

AutoDock actually consists of two main programs: AutoDock performs
the docking of the ligand to a set of grids describing the target
protein; AutoGrid pre-calculates these grids.

In addition to using them for docking, the atomic affinity grids
can be visualised. This can help, for example, to guide organic
synthetic chemists design better binders.

We have also developed a graphical user interface called AutoDockTools,
or ADT for short, which amongst other things helps to set up which
bonds will treated as rotatable in the ligand and to analyze dockings.

AutoDock has applications in:

- X-ray crystallography; 
- structure-based drug design; 
- lead optimization; 
- virtual screening (HTS); 
- combinatorial library design; 
- protein-protein docking; 
- chemical mechanism studies.

## Web site and email help for AutoDock Suite

http://autodock.scripps.edu

autodock@scripps.edu

# Installation 

The compilation of Autogrid is now standalone and does not require Autodock. It should still be used in conjunction with Autodock or Autodock-GPU. 

There are two ways to install Autogrid, using meson and using autotools. 

## Autotools (original method)

This method should work on Linux, MacOS, and Windows (Cygwin). 

First, clone the repository and run autoreconf: 

```bash
git clone https://github.com/ccsb-scripps/AutoGrid 
cd AutoGrid 
autoreconf -i # creates all the configuration files necessary 
```
It is recommended you create a directory with the name of your operating system in which the software will be compiled, e.g Linux64. 

```bash
mkdir Linux64 
cd Linux64
../configure 
```

Then make: 

```bash
make
# make check # optional, may be currently broken
make install # optional
```

## Meson (newer method)

Recently the meson build tool has been added to compile Autogrid. It has only been tested on Linux and MacOS. 
It may or may not work on Windows. 

Meson and Ninja are required for this option. 

Installation steps: 

```bash 
git clone https://github.com/ccsb-scripps/AutoGrid
cd AutoGrid 
meson setup builddir # the builddir can be any name you want. 
cd builddir
meson compile 
```

# Documentation

For documentation of AutoDock and AutoGrid see: https://autodock.scripps.edu/documentation/documentation/


# Old installation instructions

These are the contents of the old README file. They're somewhat out of date and should not be used unless you have an older version of Autodock/Autogrid.  


### Installing AutoDock

The web site offers pre-compiled executables for several
popular computers and operating systems.  If you wish
to modify AutoDock or build it for other platforms, see
"Building AutoDock from Source Tar File" below.


#### Building AutoDock from source tar File

Download the "src" tar file from the web site.
Unzip the downloaded tar file into a directory of your choice;
this will create directories src/autodock and src/autogrid.

We suggest creating subdirectories within each of "autodock" and "autogrid"
to do the compilations in.  We suggest naming them after the computer
architecture and operating system, which we abbreviate as 'ARCHOSV'.

For ARCHOSV, substitute your
computer type and operating system version, such as "i86Linux2".
You can use any name you like; nothing in the build process depends
on the name. We suggest following these patterns:

- i86Linux2 : 32-bit Intel Linux
- x86_64Linux2 : 64-bit Intel Linux
- MacOSX : Mac OS X PPC or Intel
- i86Windows : Microsoft Windows Intel
- sun4SunOS5 : SPARC Solaris

The general process is:

First, build autodock:

```bash
 "cd" # to autodock source directory
 "autoreconf -i"  # [ optional, see below: "Modifying AutoDock and AutoGrid" ]
 "mkdir" # and "cd" to your architecture's build directory (such as x86_64Linux3)
 ../configure  #   (* see below)
 make
 ```

On Microsoft Windows when building for MINGW within Cygwin, 
you need to type here:   rm .deps/*; ../configure
to run "configure" a second time for reasons we do not understand.

```bash
make check  # (optional but recommended, requires "python")
make install # (optional, this will install autodock4 executable)
```

Second, in the src/autogrid directory, do the same steps  to configure, build, and install autogrid4.

```bash
 "cd" # to autogrid source directory
 "autoreconf -i"  # [ optional, see below: "Modifying AutoDock and AutoGrid" ]
 "mkdir" # and "cd" your architecture's build directory (such as x86_64Linux3)
 ../configure  # (* see below)
 make 
 make check 
 make install # (optional)
```

On Mac OS X, you can instead use "../configure-universalDarwin"
if you wish to compile a "universal binary" for PPC/Intel 32/64 bit
as appropriate for your OS X version.

### Troubleshooting Compilations

If "autoreconf -i" hangs or runs very slowly with messages:

```
autom4te: cannot lock autom4te.cache/requests with mode 2: Input/output error
```
this is caused by building on a remote-mounted NFS file system
that does not support certain kinds of file locking.  You can
either wait for the autoreconf to finish (if it does) or do the
autoreconf on the NFS host computer, where the files are local.

### Copying AutoDock and AutoGrid

Please refer to the file "COPYING" in this directory for
more information on copying AutoDock.


### Configuration notes


By default, AutoDock 4.2 uses the AD4.1_bound.dat parameter file,
so AD4.1_compact.dat and AD4.1_extended.dat would not be used,
but if "../autodock/paramdat2h.csh" is changed to use either of
these 4.1 scoring functions, then either of these other two ".dat"
files would become necessary:

```
../autodock/AD4.1_compact.dat
../autodock/AD4.1_extended.dat
```

If you want to build in a custom parameter file to be used
as the default, such as for a new atom type "K", add your
parameters to the end of AD4.1_bound.dat and save as (for example)
/home/myusername/AD4.1_bound_plusK.dat

In the autodock source directory, edit Makefile.am the same way (about line 646)
from
```
csh $(srcdir)/paramdat2h.csh $(srcdir)/AD4_parameters.dat $(srcdir)/AD4.1_bound.dat > $@
```

to

```
csh $(srcdir)/paramdat2h.csh $(srcdir)/AD4_parameters.dat /home/myusername/AD4.1_bound_plusK.dat > $@
```

run "autoreconf -i; cd $ARCHOSV; ../configure;  make"
where $ARCHOSV is your computer type, such as "i86Linux2"


In the autogrid source directory, edit Makefile.am and change the line (about line 159)

```
csh $(srcdir)/../autodock/paramdat2h.csh $(srcdir)/../autodock/AD4_parameters.dat  $(srcdir)/../autodock/AD4.1_bound.dat > $@
```

to

```
csh $(srcdir)/../autodock/paramdat2h.csh  $(srcdir)/../autodock/AD4_parameters.dat /home/myusername/AD4.1_bound_plusK.dat > $@
```

and again run "autoreconf; cd $ARCHOSV; ../configure;  make"


### Modifying AutoDock and AutoGrid


Please feel free to change the programs to suit your needs.
If you have suggestions for improvements, please email us.


If you add source files you will need to modify Makefile.am to include any
new source files.  If you add source files or change the dependencies
of one source file upon another, you will need to generate a new 'configure'
file (for AutoDock or for AutoGrid), using 'autoreconf' as above.

