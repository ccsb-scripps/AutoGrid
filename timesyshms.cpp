/*

 $Id: timesyshms.cpp,v 1.6 2005/09/27 22:58:40 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* timesyshms.cc */

#include <stdio.h>
#include <sys/types.h>

#ifndef _WIN32
#include <sys/times.h>
#include <unistd.h>
#else
#include "times.h"
#endif
#include "timesyshms.h"

#include <time.h>


extern  FILE    *logFile;
extern	FloatOrDouble	idct;

/*----------------------------------------------------------------------------*/

void timesyshms( Clock     duration,
                 struct tms  *start,
                 struct tms  *end)

/*----------------------------------------------------------------------------*/

{
    int   h, m;
    FloatOrDouble t, T, s;
#ifndef USE_DOUBLE
    const FloatOrDouble min = 60., hrs = 3600.;
#else
    const FloatOrDouble min = 60.L, hrs = 3600.L;
#endif

    (void)fprintf( logFile, "Real= " );
    t = (FloatOrDouble)duration * idct;
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
#ifndef USE_DOUBLE
            (void)fprintf(logFile,       "%.2fs",       s );
#else
            (void)fprintf(logFile,       "%.2lfs",       s );
#endif
        else
#ifndef USE_DOUBLE
            (void)fprintf(logFile,    "%dm %05.2fs",    m, s );
#else
            (void)fprintf(logFile,    "%dm %05.2lfs",    m, s );
#endif
    } else {
#ifndef USE_DOUBLE
            (void)fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
#else
            (void)fprintf(logFile, "%dh %02dm %05.2lfs", h, m, s );
#endif
    }

    (void)fprintf( logFile, ",  CPU= " );
    t =      (FloatOrDouble)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
#ifndef USE_DOUBLE
            (void)fprintf(logFile,       "%.2fs",       s );
#else
            (void)fprintf(logFile,       "%.2lfs",       s );
#endif
        else
#ifndef USE_DOUBLE
            (void)fprintf(logFile,    "%dm %05.2fs",    m, s );
#else
            (void)fprintf(logFile,    "%dm %05.2lfs",    m, s );
#endif
    } else {
#ifndef USE_DOUBLE
            (void)fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
#else
            (void)fprintf(logFile, "%dh %02dm %05.2lfs", h, m, s );
#endif
    }

    (void)fprintf( logFile, ",  System= " );
    t = (FloatOrDouble)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;
    if (h == 0) {
        if (m == 0)
#ifndef USE_DOUBLE
            (void)fprintf(logFile,       "%.2fs",       s );
#else
            (void)fprintf(logFile,       "%.2lfs",       s );
#endif
        else
#ifndef USE_DOUBLE
            (void)fprintf(logFile,    "%dm %05.2fs",    m, s );
#else
            (void)fprintf(logFile,    "%dm %05.2lfs",    m, s );
#endif
    } else {
#ifndef USE_DOUBLE
            (void)fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
#else
            (void)fprintf(logFile, "%dh %02dm %05.2lfs", h, m, s );
#endif
    }

    (void)fprintf( logFile, "\n" );
}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
