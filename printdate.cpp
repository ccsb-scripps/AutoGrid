#include <stdio.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/time.h>
#else
#include "times.h"
#endif

#include <time.h>

#include "autogrid.h"

/* #include "autodock.h" */

void printdate( FILE *fp,
                int flag)
{
	time_t	  time_now;
	struct tm *ts;

	time_now = time(&time_now);

	ts = localtime(&time_now);
	
	if (flag) {
	    fprintf(fp, "%d:%02d %02d\" %s, %02d/%02d/%02d\n", 
	    ( (ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour ), ts->tm_min, ts->tm_sec,
	    ( (ts->tm_hour >= 12) ? "p.m." : "a.m." ),
	    (ts->tm_mon + 1), ts->tm_mday, ts->tm_year );
	} else {
	    fprintf(fp, "%d:%02d %02d\" %s\n", 
	    ( (ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour ), ts->tm_min, ts->tm_sec,
	    ( (ts->tm_hour >= 12) ? "pm" : "am" ) );
	}
}
  
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
