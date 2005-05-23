/*
  $Id: check_size.cpp,v 1.6 2005/05/23 22:16:19 gillet Exp $
*/

#include <iostream>
#include <math.h>
#include "autogrid.h"


extern char *programname;
extern FILE *logFile;

/*----------------------------------------------------------------------------*/
int check_size(int nelements, 
	       char axischar)

/*----------------------------------------------------------------------------*/

/******************************************************************************/
/*      Name: check_size                                                      */
/*  Function: Checks that number of grid elements is valid.                   */ 
/* Copyright: (C) 1994, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 13/07/92                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: nelements, axischar                                             */
/*   Returns: nelements                                                       */
/*   Globals: MAX_GRID_PTS                                                    */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 04/01/93 GMM     Created for use in makefile.                              */
/******************************************************************************/

{
    int oldnelements;

    if (nelements < 0) {
        fprintf(stderr, "\n%s: Error! Negative number of %c-grid elements!  Aborting.\n\n", programname, axischar);
        exit(-2);
    }
    if (nelements == 0) {
        fprintf(stderr, "\n%s: Warning! 0 %c-grid elements!\n\n", programname, axischar);
    }
    if (nelements>MAX_GRID_PTS) {
        fprintf(logFile, "%s: Warning! Maximum number of %c-grid elements allowed is %d. Using this value.\n", programname, axischar, MAX_GRID_PTS);
        nelements = MAX_GRID_PTS;
    }
    oldnelements = nelements;
    nelements = (int) ((nelements/2) * 2); // N.B.: integer divide truncates remainder.
    if (oldnelements != nelements)
        fprintf(logFile, "%s: Number of grid elements must be even; %c-elements changed to: %d\n", programname, axischar, nelements);

    return nelements;
}
 
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
