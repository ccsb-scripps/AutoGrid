/*

 $Id: banner.c,v 1.2 2003/02/12 19:32:29 lindy Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* banner.c */

#include <stdio.h>
#include "autogrid.h"

extern FILE *logFile;

void banner( double version_num )

{

/*----------------------------------------------------------------------------*/
/* Output banner...                                                           */
/*----------------------------------------------------------------------------*/

(void) fprintf(logFile,"\n       _______________________________________________________\n");
(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"__________//____________________________/////_________________/________\n");
(void) fprintf(logFile,"_________/__/____________/_____________/______________/_______/________\n");
(void) fprintf(logFile,"________/____/___________/_____________/_/////________________/________\n");
(void) fprintf(logFile,"________/____/__/_____/_/////___/////__/_____/_/_///__/__////_/________\n");
(void) fprintf(logFile,"_______/______/_/_____/__/_____/_____/_/_____/_//___/_/_/____//________\n");
(void) fprintf(logFile,"_______////////_/_____/__/_____/_____/_/_____/_/______/_/_____/________\n");
(void) fprintf(logFile,"_______/______/_/____//__/___/_/_____/_/_____/_/______/_/____//________\n");
(void) fprintf(logFile,"_______/______/__////_/___///___/////___/////__/______/__////_/________\n");
(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"       _______________________________________________________\n");

(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"                                ______\n");
(void) fprintf(logFile,"                               /      \\\n");
(void) fprintf(logFile,"                              /        \\\n");
(void) fprintf(logFile,"                             /          \\\n");
(void) fprintf(logFile,"                             \\    /\\    /\n");
(void) fprintf(logFile,"                              \\  /  \\  /\n");
(void) fprintf(logFile,"                               \\/ /\\ \\/\n");
(void) fprintf(logFile,"                                 /  \\\n");
(void) fprintf(logFile,"                                /____\\\n");
(void) fprintf(logFile,"\n");

(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"                  ______________________________________ \n");
(void) fprintf(logFile,"                 |                                      |\n");
(void) fprintf(logFile,"                 |            AutoGrid %3.2lf             |\n",version_num);
(void) fprintf(logFile,"                 |                                      |\n");
(void) fprintf(logFile,"                 |        Garrett M. Morris, TSRI       |\n");
(void) fprintf(logFile,"                 |        David S. Goodsell, TSRI       |\n");
(void) fprintf(logFile,"                 |        Arthur J. Olson, TSRI         |\n");
(void) fprintf(logFile,"                 |                                      |\n");
(void) fprintf(logFile,"                 |           (c) 2000, TSRI             |\n");
(void) fprintf(logFile,"                 |   The Scripps Research Institute     |\n");
(void) fprintf(logFile,"                 |______________________________________|\n");
(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"                  ______________________________________ \n");
(void) fprintf(logFile,"                 |                                      |\n");
(void) fprintf(logFile,"                 | Calculation of Non-Bond and Electro- |\n");
(void) fprintf(logFile,"                 | static Energy Grid Maps for AutoDock |\n");
(void) fprintf(logFile,"                 |   Including Solvation Free Energy    |\n");
(void) fprintf(logFile,"                 |       and Optional Disordered-H      |\n");
(void) fprintf(logFile,"                 |______________________________________|\n");
(void) fprintf(logFile,"\n");

}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
