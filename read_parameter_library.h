#ifndef _READ_PARAMETER_LIBRARY
#define _READ_PARAMETER_LIBRARY

/* $Id: read_parameter_library.h,v 1.1 2005/12/12 22:55:09 garrett Exp $ */

#include "autocomm.h"

void read_parameter_library(
        char FN_parameter_library[MAX_CHARS],
        int outlev
        );

void setup_parameter_library(
        int outlev
        );

void setup_distdepdiel( int outlev, 
                        EnergyTables *ptr_ad_energy_tables  // Holds vdw+Hb, desolvation & dielectric lookup tables
                      );


#endif
