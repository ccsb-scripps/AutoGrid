/*

 $Id: targetfile.h,v 1.1 2020/05/21 15:26:38 mp Exp $

 AutoDock  

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

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

/* targetfile.h */

#ifndef _TARGETFILE
#define _TARGETFILE

int
target_file_open ( FILE * dpf_in, FILE ** dpf_out, char *dpf_out_filename, 
	const int outLev, FILE *logFile );


int
target_file_remove_tfiles ( const char *path );

int
target_file_capability ();
#endif /*_TARGETFILE*/
