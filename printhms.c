/*
  $Id: printhms.c,v 1.2 2004/11/23 22:07:45 lindy Exp $
*/

#include <stdio.h>
#include "autogrid.h"


extern FILE *logFile;

void printhms( float t )

{
    int   h, m;
    float T, s;
    float hrs = 3600., min = 60.;

    h = (int) (t/hrs);
    T = t - h*hrs;
    m = (int) (T/min);
    s = T - m*min;

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
