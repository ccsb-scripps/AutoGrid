#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "autogrid.h"

extern FILE *GPF;
extern FILE *logFile;
extern char *programname;
extern char AutoGridHelp[];
extern char grid_param_fn[];
extern int  debug;
extern int  oldpdbq;

/*----------------------------------------------------------------------------*/

int setflags( int argc, char **argv )

/*----------------------------------------------------------------------------*/

/******************************************************************************/
/*      Name: setflags                                                        */
/*  Function: read flags from argv; return argindex of first non arg.         */
/* Copyright: (C) Garrett Matthew Morris, TSRI.                               */
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
        switch(argv[1][1]){
#ifdef FOO
        case 'n':
            ncount = atoi(argv[2]);
            argv++;
            argc--;
            argindex++;
            break;
#endif
        case 'o':
            oldpdbq = TRUE;
            break;
        case 'd':
            debug++;
            break;
        case 'u':
	    fprintf(stderr, "usage: %s %s\n", programname, AutoGridHelp);
	    exit(0);
            break;
        case 'l':
            if ( (logFile = fopen(argv[2], "w")) == NULL ) {
                fprintf(stderr, "\n%s: Sorry, I can't create the log file \"%s\"\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                exit(911);
            }
            argv++;
            argc--;
            argindex++;
            break;
        case 'p':
            strcpy(grid_param_fn, argv[2]);
            if ( (GPF = fopen(argv[2], "r")) == NULL ) {
                fprintf(stderr, "\n%s: Sorry, I can't find or open Grid Parameter File \"%s\"\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                exit(911);
            }
            argv++;
            argc--;
            argindex++;
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
    return(argindex);
}
  
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
