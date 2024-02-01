/*

 $Id: stop.cc,v 1.12 2020/05/21 15:35:51 mp Exp $

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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "stop.h"
#include "targetfile.h" 

using std::string;

extern char *programname;
extern FILE *logFile;


/*----------------------------------------------------------------------------*/
void stop(const char *const reason)
/*----------------------------------------------------------------------------*/
{
    
    /* remove any files created by extraction from a zip archive, see targetfile.cc
     * This call is harmless if no files were created 
     */
    (void) target_file_remove_tfiles("all");

    if (logFile == stdout) {
	fprintf( logFile, "%s: FATAL ERROR: %s\n", programname, reason);
	fprintf( logFile, "%s: Unsuccessful Completion.\n\n", programname);
	fflush(  logFile  );
    } else {
        string pn=programname;
        string r=reason;
	string message =  pn + ": FATAL ERROR: "+r+"\n";
	fprintf( logFile,  "%s", message.c_str() );
	fprintf( stderr, "%s", message.c_str() );

	message = pn + ": Unsuccessful Completion.\n\n";
	fprintf( logFile,  "%s", message.c_str() );
	fprintf( stderr, "%s", message.c_str() );

	fflush(logFile);
	fflush(stderr);
    }

    exit(EXIT_FAILURE); // POSIX, defined in stdlib.h
}
/*----------------------------------------------------------------------------*/
void  stop_if_nonnumeric(const char *const str, const char *const reason)
/*----------------------------------------------------------------------------*/
{
double junk1;
char junk2[200];
     if (1 != sscanf (str, "%lf%s", &junk1, junk2) ) {
	char msg[200];
	snprintf(msg, sizeof(msg), 
   	   "%s: Cannot interpret '%s' as a number.", reason, str);
	stop(msg);
	}
}
/* EOF */
