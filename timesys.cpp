/* timesys.c */

#include <stdio.h>
#include <sys/types.h>

#ifndef _WIN32
#include <sys/times.h>
#include <unistd.h>
#else
#include "times.h"
#endif

#include <time.h>
#include "autogrid.h"

extern  FILE    *logFile;
extern	float	idct;

/*----------------------------------------------------------------------------*/

void timesys( Clock       duration,
              struct tms  *start,
              struct tms  *end)

/*----------------------------------------------------------------------------*/

{
	fprintf( logFile, "%.2f %.2f %.2f rus\n",   (float) duration * idct,
                         (float)(end->tms_utime  - start->tms_utime) * idct,
                         (float)(end->tms_stime  - start->tms_stime) * idct );
}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
