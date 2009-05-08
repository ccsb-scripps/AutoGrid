/*

 $Id: times2.cpp,v 1.3 2009/05/08 23:17:35 rhuey Exp $

 AutoGrid 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 All Rights Reserved.

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

// times.c
// 
// emulate UNIX times() function for CPU time

#include <time.h>

clock_t times (struct tms *buffer)
{
	return clock();
}
//The times function stores the processor time information for the calling process in buffer. 
//The return value is the same as the value of clock(): the elapsed real time relative to an arbitrary base. The base is a constant within a particular process, and typically represents the time since system start-up. A value of (clock_t)(-1) is returned to indicate failure. 
//Portability Note: The clock function described in section Basic CPU Time Inquiry
