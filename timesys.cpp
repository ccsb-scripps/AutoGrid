/* timesys.c */

#include <stdio.h>
#include <sys/types.h>
#include <sys/times.h>
extern	float	idct;
#include <unistd.h>
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