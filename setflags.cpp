/*

 $Id: setflags.cpp,v 1.29 2020/05/07 21:25:29 mp Exp $

 AutoGrid 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "autogrid.h"
#include "constants.h"
#include <unistd.h>
#include <stdlib.h> // POSIX definitions of EXIT_SUCCESS and EXIT_FAILURE

extern FILE *GPF;
extern char *programname;
static char    AutoGridHelp[] = "\t-p parameter_filename\n\t\t\t-l log_filename\n\t\t\t-d (increment debug level)\n\t\t\t-h (display this message)\n\t\t\t--version (print version information, copyright, and license)\n";
extern char grid_param_fn[];
extern int  debug;

/*----------------------------------------------------------------------------*/

int setflags( int argc, char **argv, char *version, int use_bhtree, int use_omp, int maxthreads, FILE **logFile /* may be modified here */ )

/*----------------------------------------------------------------------------*/

/******************************************************************************/
/*      Name: setflags                                                        */
/*  Function: read flags from argv; return argindex of first non arg.         */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Matthew Morris, TSRI.                                   */
/*            (Adapted from code supplied by Bruce Duncan, TSRI.)             */
/*      Date: 06/11/92                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: argc,argv,version                                               */
/*   Returns: argindex                                                        */
/*   Globals: *GPF;                                                           */
/*            *programname;                                                   */
/*            grid_param_fn[];                                                */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/11/92 GMM     Modified for Autogrid flags:                              */
/*                  -p = Parameter filename;                                  */
/*                  -l = Log filename;                                        */
/*                  -o = Use old PDBq format (q in columns 55-61)             */
/* 04/01/93 GMM     Created for use in makefile.                              */
/******************************************************************************/

{
    int argindex;
/*----------------------------------------------------------------------------*/
/* Initialize                                                                 */
/*----------------------------------------------------------------------------*/
    argindex = 1;
    programname = strdup(argv[0]);
    GPF = stdin;
    *logFile = stderr;
/*----------------------------------------------------------------------------*/
/* Loop over arguments                                                        */
/*----------------------------------------------------------------------------*/
    while((argc > 1) && (argv[1][0] == '-')){
        if (argv[1][1] == '-') argv[1]++;

        switch(argv[1][1]){
#ifdef FOO
	// sample code
        case 'n':
            ncount = atoi(argv[2]);
            argv++;
            argc--;
            argindex++;
            break;
#endif
        case 'd':
            debug++;
            break;
        case 'u':
        case 'h':
	        fprintf(stdout, "usage: AutoGrid %s\n", AutoGridHelp);
	        exit(EXIT_SUCCESS);
            break;
        case 'l':
            if (argc < 3){
                fprintf(stderr, "\n%s: Sorry, -l requires a filename.\n\t%s\n", programname, AutoGridHelp);
                exit(EXIT_FAILURE);
            }
            if ( (*logFile = ad_fopen(argv[2], "w", *logFile)) == NULL ) {
                fprintf(stderr, "\n%s: Sorry, I can't create the log file \"%s\"\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                exit(EXIT_FAILURE);
            }
            argv++;
            argc--;
            argindex++;
            break;
        case 'p':
            if (argc < 3){
                fprintf(stderr, "\n%s: Sorry, -p requires a filename.\n\t%s\n", programname, AutoGridHelp);
                exit(EXIT_FAILURE);
            }
            strncpy(grid_param_fn, argv[2], PATH_MAX);
            grid_param_fn[PATH_MAX-1] = '\0';

            if ( (GPF = ad_fopen(argv[2], "r", *logFile)) == NULL ) {
                fprintf(stderr, "\n%s: Sorry, I can't find or open Grid Parameter File \"%s\"\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                exit(EXIT_FAILURE);
            }
            argv++;
            argc--;
            argindex++;
            break;
        case 'v':
            fprintf(stdout, "AutoGrid %-8s\n", version);
	    fprintf(stdout, "compilation options:\n");
            fprintf(stdout, "  Double-precision calculations (USE_DOUBLE): ");
#ifdef USE_DOUBLE
	    fprintf(stdout, " yes\n");
#else
	    fprintf(stdout, " no\n");
#endif
	    fprintf(stdout, "  Non-bond cutoff for internal energy calculation (SOFTNBC): %.2f Ang.\n", SOFTNBC);
	    fprintf(stdout, "  Hard Non-bond cutoff for vdW/Hb calculation (NBC): %.2f Ang.\n", NBC);
	    fprintf(stdout, "  Electrostatic/desolvation cutoff (DIELCUTOFF): %.2f Ang.\n", DIELCUTOFF);
            fprintf(stdout, "  Optimize internal energy scoring (USE_8A_NBCUTOFF): ");
#ifdef USE_8A_NBCUTOFF
	    fprintf(stdout, " yes\n");
#else
	    fprintf(stdout, " no\n");
#endif
            fprintf(stdout, "  Faster search for nearby atoms (USE_BHTREE): ");
	    fprintf(stdout, use_bhtree?" yes\n":" no\n");
            fprintf(stdout, "  Run calculations in parallel if possible (_OPENMP): ");
	    fprintf(stdout, use_omp?" yes\n":" no\n");
	    if ( use_omp) fprintf(stdout, "  Maximum number of parallel threads (MAXTHREADS): %d\n", maxthreads);
	    fprintf(stdout, "  Maximum number of receptor atoms (AG_MAX_ATOMS): %d\n", AG_MAX_ATOMS);
	    fprintf(stdout, "  Maximum number of receptor atom types (NUM_RECEPTOR_TYPES): %d\n", NUM_RECEPTOR_TYPES);
	    fprintf(stdout, "  Maximum number of bonds per receptor atom (AG_MAX_NBONDS): %d\n", AG_MAX_NBONDS);
	    fprintf(stdout, "  Maximum number of atom types (MAX_ATOM_TYPES): %d\n", MAX_ATOM_TYPES);
	    fprintf(stdout, "  Maximum number of maps (MAX_MAPS): %d\n", MAX_MAPS);
	    fprintf(stdout, "  Maximum dimension of map x, y, or z (MAX_GRID_PTS): %d\n", MAX_GRID_PTS);
	    /* print sizes of key types for this compilation */
	    fprintf(stdout, "  Size of int %d, long %d, long long %d, float %d, double %d, Real %d bytes.\n",
		(int)(sizeof(int)), (int)(sizeof(long)), (int)(sizeof(long long)),
		(int)(sizeof(float)), (int)(sizeof(double)), (int)(sizeof(Real)) );


            fprintf(stdout, "\n Copyright (C) 2009 The Scripps Research Institute.\n");
// GNU BEGIN   (see maintenance script update_license_de-GNU)
            fprintf(stdout, " License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>\n");
            fprintf(stdout, " This is free software: you are free to change and redistribute it.\n");
// GNU END   (see maintenance script update_license_de-GNU)
            fprintf(stdout, " There is NO WARRANTY, to the extent permitted by law.\n");
            exit(EXIT_SUCCESS);
            break;
        default:
            fprintf(stderr,"%s: unknown switch -%c\n",programname,argv[1][1]);
            exit(EXIT_FAILURE);
            break;
        }
        argindex++;
        argc--;
        argv++;
    }
    //no gpf specified and input is terminal, very likely an error
    if (GPF==stdin && isatty(fileno(stdin))){
	    fprintf(stdout, "usage: AutoGrid %s\n", AutoGridHelp);
	    exit(EXIT_FAILURE); 
    }
    return(argindex);
}
  
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
