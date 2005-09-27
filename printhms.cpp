/*

 $Id: printhms.cpp,v 1.4 2005/09/27 22:58:40 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* printhms.cc */

    #include <stdio.h>
    #include "printhms.h"

extern FILE *logFile;

void printhms( FloatOrDouble t )

{
    int   h,
          m;
    FloatOrDouble T, s;
    FloatOrDouble min = 60.,
	  hrs = 3600.;

    h = (int)(t/hrs);
    T = t - ((FloatOrDouble)h)*hrs;
    m = (int)(T/min);
    s = T - ((FloatOrDouble)m)*min;

    if (h == 0) {
        if (m == 0)
            fprintf(logFile,       "%.2fs",       s );
        else
            fprintf(logFile,    "%dm %05.2fs",    m, s );
    } else {
            fprintf(logFile, "%dh %02dm %05.2fs", h, m, s );
    }
}
/* EOF */
