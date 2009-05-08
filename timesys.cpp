/*

 $Id: timesys.cpp,v 1.9 2009/05/08 23:17:35 rhuey Exp $

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

#ifndef _WIN32
#   include <sys/times.h>
#   include <unistd.h>
#else
#   include "times.h"
#endif

#ifdef HAVE_CONFIG_H
#   include <config.h>
#endif

#include <stdio.h>
#include "timesys.h"

extern  FILE    *logFile;
extern	Real	idct;

/*----------------------------------------------------------------------------*/

void timesys( Clock       duration,
              struct tms  *start,
              struct tms  *end)

/*----------------------------------------------------------------------------*/

{
	fprintf( logFile, "Real= %.2f,  CPU= %.2f,  System= %.2f\n",     (Real)duration * idct,
                         (Real)(end->tms_utime  - start->tms_utime) * idct,
                         (Real)(end->tms_stime  - start->tms_stime) * idct );
}
/* EOF */
