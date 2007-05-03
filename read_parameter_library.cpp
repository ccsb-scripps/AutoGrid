/*

 $Id: read_parameter_library.cpp,v 1.3 2007/05/03 20:46:06 garrett Exp $

 AutoGrid 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "autogrid.h"
#include "autocomm.h"
#include "constants.h"
#include "parameters.h"
// #include "openfile.h"
#include "parse_param_line.h"
#include "atom_parameter_manager.h"
#include "partokens.h"
#include "default_parameters.h"


extern FILE *logFile;
extern char *programname;
extern int debug;
extern Linear_FE_Model AD4;


void read_parameter_library(
        char FN_parameter_library[MAX_CHARS],
        int outlev
        )
{
    static ParameterEntry thisParameter;
    FILE *parameter_library_file;
    char parameter_library_line[MAX_CHARS];
    int nfields;
    int param_keyword = -1;
    int int_hbond_type = 0;

    pr(logFile, "Using read_parameter_library\n");

    // Open and read the parameter library
    //
    if ((parameter_library_file = ag_fopen(FN_parameter_library, "r")) == NULL) {
         fprintf(stderr,"Sorry, I can't find or open %s\n", FN_parameter_library);
         exit(-1);
    }

    while (fgets(parameter_library_line, sizeof(parameter_library_line), parameter_library_file) != NULL) {
        param_keyword = parse_param_line( parameter_library_line );
        if (debug > 0) {
            pr(logFile, "DEBUG: parameter_library_line = %sDEBUG: param_keyword          = %d\n", parameter_library_line, param_keyword);
        }

        switch (param_keyword) {
            case PAR_:
            case PAR_NULL:
            case PAR_COMMENT:
                break;

            case PAR_VDW:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_vdW);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the van der Waals term = \t%.4lf\n\n", AD4.coeff_vdW);
                break;

            case PAR_HBOND:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_hbond);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the H-bonding term     = \t%.4lf\n\n", AD4.coeff_hbond);
                break;

            case PAR_ESTAT:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_estat);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the electrostatic term = \t%.4lf\n\n", AD4.coeff_estat);
                break;

            case PAR_DESOLV:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_desolv);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the desolvation term   = \t%.4lf\n\n", AD4.coeff_desolv);
                break;

            case PAR_TORS:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_tors);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the torsional term     = \t%.4lf\n\n", AD4.coeff_tors);
                break;

            case PAR_ATOM_PAR:
                // Read in one line of atom parameters;
                // NB: scanf doesn't try to write missing fields
                nfields = sscanf(parameter_library_line, "%*s %s %lf %lf %lf %lf %lf %lf %d %d %d %d",
                                    thisParameter.autogrid_type,
                                    &thisParameter.Rij,
                                    &thisParameter.epsij,
                                    &thisParameter.vol,
                                    &thisParameter.solpar,
                                    &thisParameter.Rij_hb,
                                    &thisParameter.epsij_hb,
                                    &int_hbond_type,
                                    &thisParameter.rec_index,
                                    &thisParameter.map_index,
                                    &thisParameter.bond_index);
                if (nfields < 2) {
                    continue; // skip any parameter_library_line without enough info
                }

                if (int_hbond_type == 0) {
                    thisParameter.hbond = NON;
                } else if (int_hbond_type == 1) {
                    thisParameter.hbond = DS;
                } else if (int_hbond_type == 2) {
                    thisParameter.hbond = D1;
                } else if (int_hbond_type == 3) {
                    thisParameter.hbond = AS;
                } else if (int_hbond_type == 4) {
                    thisParameter.hbond = A1;
                } else if (int_hbond_type == 5) {
                    thisParameter.hbond = A2;
                } else {
                    thisParameter.hbond = NON;
                }

                thisParameter.epsij    *= AD4.coeff_vdW;
                thisParameter.epsij_hb *= AD4.coeff_hbond;

                apm_enter(thisParameter.autogrid_type, thisParameter);
                pr(logFile, "Parameters for the atom type named \"%s\" were read in from the parameter library as follows:\n", thisParameter.autogrid_type);

                if (outlev > 2) {
                    pr(logFile, "\tR-eqm = %5.2f Angstrom\n\tweighted epsilon = %5.3f\n\tAtomic fragmental volume = %5.3f\n\tAtomic solvation parameter = %5.3f\n\tH-bonding R-eqm = %5.3f\n\tweighted H-bonding epsilon = %5.3f\n\tH-bonding type = %d,  bond index = %d\n\n",
                            thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                            thisParameter.Rij_hb, thisParameter.epsij_hb, thisParameter.hbond, thisParameter.bond_index );
                } else {
                    pr(logFile, "\tR-eqm = %.2f Angstrom,  weighted epsilon = %.3f,  At.frag.vol. = %.3f,  At.solv.par. = %.3f, \n\tHb R-eqm = %.3f,  weighted Hb epsilon = %.3f,  Hb type = %d,  bond index = %d\n\n",
                            thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                            thisParameter.Rij_hb, thisParameter.epsij_hb, thisParameter.hbond, thisParameter.bond_index );
                }
                break;

            default:
                break;
        } // switch
    } // while there is another line of parameters to read in
}

void setup_parameter_library( int outlev )
{
    static ParameterEntry thisParameter;
    char parameter_library_line[MAX_CHARS];
    int nfields;
    int param_keyword = -1;
    int int_hbond_type = 0;
    register int counter = 0;

    pr(logFile, "Setting up parameter library with factory defaults.\n\n\n");

    // Default parameters
    //
    // These are set up in "default_parameters.h"
    // and stored in the param_string[MAX_LINES] array

    while ( param_string[counter] != NULL) {
        param_keyword = parse_param_line( param_string[counter] );

        (void)strcpy(parameter_library_line, param_string[counter]);
        counter++;
        if (debug > 0) {
            pr(logFile, "DEBUG: parameter_library_line = %sDEBUG: param_keyword          = %d\n", parameter_library_line, param_keyword);
        }

        switch (param_keyword) {
            case PAR_:
            case PAR_NULL:
            case PAR_COMMENT:
                break;

            case PAR_VDW:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_vdW);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the van der Waals term = \t%.4lf\n\n", AD4.coeff_vdW);
                break;

            case PAR_HBOND:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_hbond);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the H-bonding term     = \t%.4lf\n\n", AD4.coeff_hbond);
                break;

            case PAR_ESTAT:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_estat);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the electrostatic term = \t%.4lf\n\n", AD4.coeff_estat);
                break;

            case PAR_DESOLV:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_desolv);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the desolvation term   = \t%.4lf\n\n", AD4.coeff_desolv);
                break;

            case PAR_TORS:
                nfields = sscanf(parameter_library_line, "%*s %lf", &AD4.coeff_tors);
                if (nfields < 1) {
                    pr( logFile, "%s: WARNING:  Please supply a coefficient as a floating point number.\n\n", programname);
                    continue; // skip any parameter_library_line without enough info
                }
                pr( logFile, "Free energy coefficient for the torsional term     = \t%.4lf\n\n", AD4.coeff_tors);
                break;

            case PAR_ATOM_PAR:
                // Read in one line of atom parameters;
                // NB: scanf doesn't try to write missing fields
                nfields = sscanf(parameter_library_line, "%*s %s %lf %lf %lf %lf %lf %lf %d %d %d %d",
                                    thisParameter.autogrid_type,
                                    &thisParameter.Rij,
                                    &thisParameter.epsij,
                                    &thisParameter.vol,
                                    &thisParameter.solpar,
                                    &thisParameter.Rij_hb,
                                    &thisParameter.epsij_hb,
                                    &int_hbond_type,
                                    &thisParameter.rec_index,
                                    &thisParameter.map_index,
                                    &thisParameter.bond_index);
                if (nfields < 2) {
                    continue; // skip any parameter_library_line without enough info
                }

                if (int_hbond_type == 0) {
                    thisParameter.hbond = NON;
                } else if (int_hbond_type == 1) {
                    thisParameter.hbond = DS;
                } else if (int_hbond_type == 2) {
                    thisParameter.hbond = D1;
                } else if (int_hbond_type == 3) {
                    thisParameter.hbond = AS;
                } else if (int_hbond_type == 4) {
                    thisParameter.hbond = A1;
                } else if (int_hbond_type == 5) {
                    thisParameter.hbond = A2;
                } else {
                    thisParameter.hbond = NON;
                }

                thisParameter.epsij    *= AD4.coeff_vdW;
                thisParameter.epsij_hb *= AD4.coeff_hbond;

                apm_enter(thisParameter.autogrid_type, thisParameter);
                pr(logFile, "Parameters for the atom type named \"%s\" were read in from the parameter library as follows:\n", thisParameter.autogrid_type);

                if (outlev > 2) {
                    pr(logFile, "\tR-eqm = %5.2f Angstrom\n\tweighted epsilon = %5.3f\n\tAtomic fragmental volume = %5.3f\n\tAtomic solvation parameter = %5.3f\n\tH-bonding R-eqm = %5.3f\n\tweighted H-bonding epsilon = %5.3f\n\tH-bonding type = %d,  bond index = %d\n\n",
                            thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                            thisParameter.Rij_hb, thisParameter.epsij_hb, thisParameter.hbond, thisParameter.bond_index );
                } else {
                    pr(logFile, "\tR-eqm = %.2f Angstrom,  weighted epsilon = %.3f,  At.frag.vol. = %.3f,  At.solv.par. = %.3f, \n\tHb R-eqm = %.3f,  weighted Hb epsilon = %.3f,  Hb type = %d,  bond index = %d\n\n",
                            thisParameter.Rij, thisParameter.epsij, thisParameter.vol, thisParameter.solpar,
                            thisParameter.Rij_hb, thisParameter.epsij_hb, thisParameter.hbond, thisParameter.bond_index );
                }
                break;

            default:
                break;
        } // switch
    } // while there is another line of parameters to read in
}

/* EOF */
