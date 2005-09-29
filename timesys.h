

#ifndef TIMESYS
#define TIMESYS
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include "autocomm.h"
void  timesys( Clock  duration,
               struct tms *start,
               struct tms *end );
#endif
