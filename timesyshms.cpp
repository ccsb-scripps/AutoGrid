/*

 $Id: timesyshms.cpp,v 1.7 2006/07/12 21:24:13 garrett Exp $

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
extern	Real	idct;

/*----------------------------------------------------------------------------*/

void timesyshms( Clock     duration,
                 struct tms  *start,
                 struct tms  *end)

/*----------------------------------------------------------------------------*/

{
    int   h, m;
    Real t, T, s;
#ifndef USE_DOUBLE
    const Real min = 60., hrs = 3600.;
#else
    const Real min = 60.L, hrs = 3600.L;
#endif

    (void)fprintf( logFile, "Real= " );
    t = (Real)duration * idct;
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
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
    t =      (Real)((end->tms_utime  - start->tms_utime) * idct);
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
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
    t = (Real)((end->tms_stime  - start->tms_stime) * idct);
    h = (int)(t/hrs);
    T = t - ((Real)h)*hrs;
    m = (int)(T/min);
    s = T - ((Real)m)*min;
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
