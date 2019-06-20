/*

 $Id: bondmanager.cpp,v 1.1 2019/06/20 19:07:27 mp Exp $

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include "autogrid.h"
#include "bondmanager.h"

// bondmanager for AutoGrid - adapted from AutoDock nonbonds.cc - M Pique 2018
// functions: addbond, isbonded

extern char *programname;

static bool debug=false;

// addbond - returns 0 if OK, else non-zero for error
int
addbond(
              const int i, // from_atom
              const int j, // to_atom
	      /* not const */ int nbonds[],
              /* not const */ int bonded[][AG_MAX_NBONDS],
	      const int outlev,
	      FILE *logFile
	      )
{
    int errorcode = 0;

	            // Remember:   int bonded[AG_MAX_ATOMS][AG_MAX_NBONDS];	

                    if (nbonds[i] >= AG_MAX_NBONDS || nbonds[j] >= AG_MAX_NBONDS) {
                    if (nbonds[i] >= AG_MAX_NBONDS) {
                        // bonded array for the i-th atom is full; return failure code here.
                        fprintf( logFile, "%s: WARNING!  Atom %d has too many bonds %d (%d is the maximum): ", programname, i+1, nbonds[i], AG_MAX_NBONDS );
                        for(int k=0;k<AG_MAX_NBONDS;k++) fprintf( logFile, "%d ", bonded[i][k]+1);
                        fprintf( logFile, "[%d...]\n", j+1);
			errorcode |= 1;
			}
                    if (nbonds[j] >= AG_MAX_NBONDS) {
                        // bonded array for the j-th atom is full; return failure code here.
                        fprintf( logFile, "%s: WARNING!  Atom %d has too many bonds %d (%d is the maximum): ", programname, j+1, nbonds[j], AG_MAX_NBONDS );
                        for(int k=0;k<AG_MAX_NBONDS;k++) fprintf( logFile, "%d ", bonded[j][k]+1);
                        fprintf( logFile, "[%d...]\n", i+1);
			errorcode |= 2;
			}

		    } else {
                // add a symmetric bond between i and j
				bonded[i][ nbonds[i] ] = j;
				bonded[j][ nbonds[j] ] = i;
                // increment the number of bonds to i and j
				nbonds[i] += 1;
				nbonds[j] += 1;
		   }

                if (debug) {
                    // print out details
                    (void)fprintf(logFile,"Adding a bond between %d and %d \n", i+1, j+1);
                    (void)fprintf(logFile,"  bonded[%d][ nbonds=%d ] = %d\n", i+1, nbonds[i], j+1);
                    (void)fprintf(logFile,"  bonded[%d][ nbonds=%d ] = %d\n", j+1, nbonds[j], i+1);
                    (void)fprintf(logFile,"  nbonds[%d]=%d\n", i+1, nbonds[i]);
                    (void)fprintf(logFile,"  nbonds[%d]=%d\n", j+1, nbonds[j]);
                }
				
	return errorcode;
} // end of get bonds

/*----------------------------------------------------------------------------*/

bool isbonded(const int i, const int j, const int nbonds[], const int bonded[][AG_MAX_NBONDS], 
const int outlev, FILE *logFile)
{

// return true iff atom i is bonded to atom j or atom j is bonded to atom i
    return isbonded1(i, j, nbonds, bonded, outlev, logFile) ||
      isbonded1(j, i, nbonds, bonded, outlev, logFile);
}

bool isbonded1(const int i, const int j, const int nbonds[], const int bonded[][AG_MAX_NBONDS], 
const int outlev, FILE *logFile)
{

// return true if atom i is bonded to atom j (or if j i
    if(debug) pr(logFile, "DEBUG: isbonded  atom %d to atom %d. \n", i+1, j+1);
    if(debug) pr(logFile, "atom %d  %d bonds: ", i+1, nbonds[i]);
    if ( i == j ) return true; // or is this an error?
    for (int ii = 0; ii<AG_MAX_NBONDS && ii<nbonds[i]; ii++) {  // loop over slots in i looking for j
        if(debug) pr(logFile, "   slot %d->atom %d ", ii, bonded[i][ii]+1);
	if (bonded[i][ii] == j) {
	    if(debug) pr(logFile, " BONDED \n");
	    return true;
	}
    }
    return false;
} // end of isbonded1

void printbonds(const int natom, const int nbonds[], const int bonded[][AG_MAX_NBONDS], 
const int outlev, FILE *logFile)
{

    //pr(logFile, "%s", message);
    for (int i = 0; i<natom; i++) {  // loop over atoms, "i", from 1 to natom
        if(debug) pr(logFile, "DEBUG:  atom %d  bonded to ", i+1);
        for (int j=0; j<AG_MAX_NBONDS;j++) {// loop over all the slots for this atom's "j" bonds
	    if(j==nbonds[i]) pr(logFile, "|");
	    pr(logFile, "  %d", bonded[i][j]+1);
	    }
        pr(logFile, "\n");
    } // i
} // end of printbonds

/* EOF */
