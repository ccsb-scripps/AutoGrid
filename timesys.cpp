/*

 $Id: timesys.cpp,v 1.6 2005/09/27 22:58:40 garrett Exp $

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
extern	FloatOrDouble	idct;

/*----------------------------------------------------------------------------*/

void timesys( Clock       duration,
              struct tms  *start,
              struct tms  *end)

/*----------------------------------------------------------------------------*/

{
	fprintf( logFile, "Real= %.2f,  CPU= %.2f,  System= %.2f\n",     (FloatOrDouble)duration * idct,
                         (FloatOrDouble)(end->tms_utime  - start->tms_utime) * idct,
                         (FloatOrDouble)(end->tms_stime  - start->tms_stime) * idct );
}
/* EOF */
