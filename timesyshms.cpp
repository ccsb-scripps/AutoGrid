/* timesyshms.c */

#include <stdio.h>
#include <sys/types.h>
#include <sys/times.h>

#include <unistd.h>

void timesyshms( Clock     duration,
                 struct tms  *start,
                 struct tms  *end)

/*----------------------------------------------------------------------------*/

{
	fprintf( logFile, "elapsed (wall clock): " );
	printhms( (float)(duration * idct) );
	fprintf( logFile, ",  user: " );
	printhms( (float)((end->tms_utime  - start->tms_utime) * idct) );
	fprintf( logFile, ",  system: " );
	printhms( (float)((end->tms_stime  - start->tms_stime) * idct) );
	fprintf( logFile, "\n" );
}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
