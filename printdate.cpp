/*

 $Id: printdate.cpp,v 1.6 2005/09/27 22:58:40 garrett Exp $

*/

#include <stdio.h>
#include <sys/types.h>

#ifndef _WIN32
#   include <sys/time.h>
#else
#   include "times.h"
#endif

#ifdef HAVE_CONFIG_H
#   include <config.h>
#endif

#include "printdate.h"

void printdate( FILE *fp, int flag )
{
    time_t tn; /* tn = "time_now" */
    char *StringTimeDate;
    struct tm *ts;

    tn = time( &tn );

    ts = localtime( &tn );
    
    if (flag==1) {
        fprintf(fp, "%d:%02d %02d\" %s, %02d/%02d/%4d\n", 
        ( (ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour ), ts->tm_min, ts->tm_sec,
        ( (ts->tm_hour >= 12) ? "p.m." : "a.m." ),
        (ts->tm_mon + 1), ts->tm_mday, 1900+ts->tm_year );
    } else if (flag==2) {
          StringTimeDate = ctime( &tn );
          fprintf(fp, "%s", StringTimeDate);
    } else {
        fprintf(fp, "%d:%02d %02d\" %s\n", 
        ( (ts->tm_hour >  12) ? (ts->tm_hour-12) : ts->tm_hour ), ts->tm_min, ts->tm_sec,
        ( (ts->tm_hour >= 12) ? "pm" : "am" ) );
    }
}
/* EOF */
