/*

 $Id: prHMSfixed.c,v 1.2 2003/02/12 19:32:29 lindy Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "autogrid.h"


extern FILE *logFile;

void prHMSfixed( float t )
{
    int   h, m;
    float T, s;
    float hrs = 3600., min = 60.;

    h = t/hrs;
    T = t - h*hrs;
    m = T/min;
    s = T - m*min;

    if (h == 0) {
        if (m == 0)
            fprintf(logFile,    "        %5.2fs",        s );
        else
            fprintf(logFile,   "    %2dm %05.2fs",    m, s );
    } else {
            fprintf(logFile, "%2dh %02dm %05.2fs", h, m, s );
    }
}
/* EOF */
