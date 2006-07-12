/*

 $Id: timesys.cpp,v 1.7 2006/07/12 21:24:13 garrett Exp $

 */

#ifndef _WIN32
#   include <sys/times.h>
#   include <unistd.h>
#else
#   include "times.h"
#endif

#ifdef HAVE_CONFIG_H
#   include <config.h>
#endif

#include <stdio.h>
#include "timesys.h"

extern  FILE    *logFile;
extern	Real	idct;

/*----------------------------------------------------------------------------*/

void timesys( Clock       duration,
              struct tms  *start,
              struct tms  *end)

/*----------------------------------------------------------------------------*/

{
	fprintf( logFile, "Real= %.2f,  CPU= %.2f,  System= %.2f\n",     (Real)duration * idct,
                         (Real)(end->tms_utime  - start->tms_utime) * idct,
                         (Real)(end->tms_stime  - start->tms_stime) * idct );
}
/* EOF */
