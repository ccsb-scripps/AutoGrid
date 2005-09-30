/*
  $Id: timesyshms.h,v 1.1.6.1 2005/09/30 22:45:22 alther Exp $
*/

#ifndef TIMESYSHMS
#define TIMESYSHMS
#include <sys/types.h>

#ifdef _WIN32
   #include  "times.h"
#else
   #include <sys/times.h>
#endif

#include <time.h>
#include "autocomm.h"
#include "printhms.h"
void  timesyshms( Clock  duration,
                  struct tms *start,
                  struct tms *end );
#endif
