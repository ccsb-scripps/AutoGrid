
#ifndef TIMESYSHMS
#define TIMESYSHMS
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include "autocomm.h"
#include "printhms.h"
void  timesyshms( Clock  duration,
                  struct tms *start,
                  struct tms *end );
#endif
