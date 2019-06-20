/*

 $Id: bondmanager.h,v 1.1 2019/06/20 19:07:27 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifndef BONDMANAGER
#define BONDMANAGER

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "autogrid.h"

// bondmanager for AutoGrid - adapted from AutoDock nonbonds.cc - M Pique 2018
// functions: addbond, isbonded, isbonded1, printbonds

// addbond - returns 0 if OK, else non-zero for error
int
addbond(
              const int i, // from_atom
              const int j, // to_atom
	      /* not const */ int nbonds[],
              /* not const */ int bonded[][AG_MAX_NBONDS],
	      const int outlev,
	      FILE *logFile
	      );

bool isbonded(const int i, const int j, const int nbonds[], const int bonded[][AG_MAX_NBONDS], 
const int outlev, FILE *logFile);

bool isbonded1(const int i, const int j, const int nbonds[], const int bonded[][AG_MAX_NBONDS], 
const int outlev, FILE *logFile);

void printbonds(const int natom, const int nbonds[], const int bonded[][AG_MAX_NBONDS], 
const int outlev, FILE *logFile);

#endif
/* EOF */
