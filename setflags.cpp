/*

 $Id: setflags.cpp,v 1.12 2009/06/30 00:26:24 rhuey Exp $

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
#include <unistd.h>

extern FILE *GPF;
extern FILE *logFile;
extern char *programname;
static char    AutoGridHelp[] = "\t-p parameter_filename\n\t\t\t-l log_filename\n\t\t\t-d (increment debug level)\n\t\t\t-h (display this message)\n\t\t\t--version (print version information, copyright, and license)\n";
extern char grid_param_fn[];
extern int  debug;

/*----------------------------------------------------------------------------*/

int setflags( int argc, char **argv )

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
/*    Inputs: argc,argv                                                       */
/*   Returns: argindex                                                        */
/*   Globals: *GPF;                                                           */
/*            *logFile;                                                       */
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
    programname = argv[0];
    GPF = stdin;
    logFile = stdout;
/*----------------------------------------------------------------------------*/
/* Loop over arguments                                                        */
/*----------------------------------------------------------------------------*/
    while((argc > 1) && (argv[1][0] == '-')){
        if (argv[1][1] == '-') argv[1]++;

        switch(argv[1][1]){
#ifdef FOO
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
	        exit(0);
            break;
        case 'l':
            if ( (logFile = ad_fopen(argv[2], "w")) == NULL ) {
                fprintf(stderr, "\n%s: Sorry, I can't create the log file \"%s\"\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                exit(-1);
            }
            argv++;
            argc--;
            argindex++;
            break;
        case 'p':
            strcpy(grid_param_fn, argv[2]);
            if ( (GPF = ad_fopen(argv[2], "r")) == NULL ) {
                fprintf(stderr, "\n%s: Sorry, I can't find or open Grid Parameter File \"%s\"\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                exit(-1);
            }
            argv++;
            argc--;
            argindex++;
            break;
        case 'v':
            fprintf(stdout, "AutoGrid %-8s\n", VERSION);
            fprintf(stdout, " Copyright (C) 2009 The Scripps Research Institute.\n");
            fprintf(stdout, " License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>\n");
            fprintf(stdout, " This is free software: you are free to change and redistribute it.\n");
            fprintf(stdout, " There is NO WARRANTY, to the extent permitted by law.\n");
            exit(0);
            break;
        default:
            fprintf(stderr,"%s: unknown switch -%c\n",programname,argv[1][1]);
            exit(1);
            break;
        }
        argindex++;
        argc--;
        argv++;
    }
    //no gpf specified and input is terminal, very likely an error
    if (GPF==stdin && isatty(fileno(stdin))){
	    fprintf(stdout, "usage: AutoGrid %s\n", AutoGridHelp);
	    exit(-1); 
    }
    return(argindex);
}
  
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
