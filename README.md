# Autogrid4 

Autogrid4 is a support software for docking programs such as AutoDock4 and Autodock-GPU. 
Its function is to precalculate the grids used by the docking software. 


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