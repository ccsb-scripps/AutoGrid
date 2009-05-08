/*

 $Id: timesyshms.cpp,v 1.9 2009/05/08 23:17:35 rhuey Exp $

 AutoGrid 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
