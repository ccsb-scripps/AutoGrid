/*
  $Id: timesys.h,v 1.1.6.1 2005/09/30 22:45:06 alther Exp $
*/

#ifndef TIMESYS
#define TIMESYS
#include <sys/types.h>

#ifdef _WIN32
   #include "times.h"
#else
   #include <sys/times.h>
#endif

#include <time.h>
#include "autocomm.h"
void  timesys( Clock  duration,
               struct tms *start,
               struct tms *end );
#endif
