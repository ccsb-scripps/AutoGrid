/*

 $Id: banner.cpp,v 1.12 2007/05/03 20:46:06 garrett Exp $

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
(void) fprintf(logFile,"________/____/___________/_____________/______________________/________\n");
(void) fprintf(logFile,"________/____/__/_____/_/////___/////__/__////_/_///__/__////_/________\n");
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
(void) fprintf(logFile,"                ______________________________________ \n");
(void) fprintf(logFile,"               |                                      |\n");
(void) fprintf(logFile,"               |            AutoGrid %3.2lf             |\n",version_num);
(void) fprintf(logFile,"               |                                      |\n");
(void) fprintf(logFile,"               |        Garrett M. Morris, TSRI       |\n");
(void) fprintf(logFile,"               |            Ruth Huey, TSRI           |\n");
(void) fprintf(logFile,"               |        David S. Goodsell, TSRI       |\n");
(void) fprintf(logFile,"               |         Arthur J. Olson, TSRI        |\n");
(void) fprintf(logFile,"               |                                      |\n");
(void) fprintf(logFile,"               |        (c) 1989-2005, TSRI           |\n");
(void) fprintf(logFile,"               |   The Scripps Research Institute     |\n");
(void) fprintf(logFile,"               |______________________________________|\n");
(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"                ______________________________________ \n");
(void) fprintf(logFile,"               |                                      |\n");
(void) fprintf(logFile,"               | Calculation of van der Waals, H-Bond,|\n");
(void) fprintf(logFile,"               |   Electrostatic Potential Energy, &  |\n");
(void) fprintf(logFile,"               |   Desolvation Free Energy Grid Maps  |\n");
(void) fprintf(logFile,"               |             for AutoDock             |\n");
(void) fprintf(logFile,"               |______________________________________|\n");
(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"\n");
(void) fprintf(logFile,"\n");

}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
