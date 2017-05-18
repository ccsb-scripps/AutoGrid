#include <stdio.h>
#include "autogrid.h"


void prHMSfixed( float t, FILE *logFile)
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
            fprintf(logFile,    "        %5.2fs",        s );
        else
            fprintf(logFile,   "    %2dm %05.2fs",    m, s );
    } else {
            fprintf(logFile, "%2dh %02dm %05.2fs", h, m, s );
    }
}
/* EOF */
