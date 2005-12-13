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
(void) fprintf(logFile,"                 |            Ruth Huey, TSRI           |\n");
(void) fprintf(logFile,"                 |        David S. Goodsell, TSRI       |\n");
(void) fprintf(logFile,"                 |         Arthur J. Olson, TSRI        |\n");
(void) fprintf(logFile,"                 |                                      |\n");
(void) fprintf(logFile,"                 |        (c) 1989-2005, TSRI           |\n");
(void) fprintf(logFile,"                 |   The Scripps Research Institute     |\n");
(void) fprintf(logFile,"                 |______________________________________|\n");
(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"                  ______________________________________ \n");
(void) fprintf(logFile,"                 |                                      |\n");
(void) fprintf(logFile,"                 | Calculation of van der Waals, H-Bond,|\n");
(void) fprintf(logFile,"                 |   Electrostatic Potential Energy, &  |\n");
(void) fprintf(logFile,"                 |   Desolvation Free Energy Grid Maps  |\n");
(void) fprintf(logFile,"                 |             for AutoDock             |\n");
(void) fprintf(logFile,"                 |______________________________________|\n");
(void) fprintf(logFile,"\n");

}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
