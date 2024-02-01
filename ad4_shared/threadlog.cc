/* thread-safe log file utility functions - M Pique, 2014 */

/*

 $Id: threadlog.cc,v 1.7 2020/07/05 16:44:32 mp Exp $

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "threadlog.h"
#include "constants.h"
#include "stop.h"
/* include stdlib.h for "free" and unistd.h for "unlink"  */ 
/* tempnam is in <stdio.h>  */
#include <stdlib.h>
#include <unistd.h>

static char **tfilename /*[max_threads]*/;
static FILE **tfileptr /*[max_threads]*/;
static int max_threads;

void
threadLogAlloc(int j)
{
	
	free(tfileptr);
	free(tfilename);
	tfilename =  (char **) calloc(j, sizeof (char*));
	tfileptr =  (FILE **) calloc(j, sizeof (FILE*));
	if(tfilename==NULL || tfileptr==NULL) stop("threadLogAlloc failure");
	max_threads=j;
	}

FILE *
threadLogOpen(int j)
{
	if(j<0) {
		char msg[200];
		sprintf(msg,"threadLogOpen thread number j %d < 0", j);
		stop(msg);
	}
	if(j>=max_threads) {
		char msg[200];
		sprintf(msg,"threadLogOpen thread number too large j %d >= max_threads %d",
				j,max_threads);
		stop(msg);
	}
	// note that tempnam does its own malloc() and is thread-safe MPique
	// We look in environment variables to choose the directory, hoping this will work on Windows
	char *tempdir = NULL;
	if (tempdir==NULL) tempdir = getenv("TMP"); 
	if (tempdir==NULL) tempdir = getenv("TMPDIR"); 
	if (tempdir==NULL) tempdir = getenv("TEMP"); 
#ifdef TEMPDIR_DEBUG
	if (tempdir!=NULL) fprintf(stderr, "tempdir chosen as '%s'\n", tempdir);
	else fprintf(stderr, "tempdir is default value NULL\n");
#endif

	tfilename[j] = tempnam(tempdir, "autod");
#ifdef TEMPDIR_DEBUG
	fprintf(stderr, "tempnam is '%s'\n", tfilename[j]);
	fflush(stderr);
#endif

	tfileptr[j] = fopen(tfilename[j], "w");

	if(NULL==tfileptr[j]) stop("cannot allocate or open temp log file");
#ifdef TEMPDIR_DEBUG
	if (tempdir!=NULL) fprintf(tfileptr[j], "thread %d tempdir chosen as '%s'\n", j, tempdir);
	else fprintf(tfileptr[j], "thread %d tempdir is default value NULL\n",j);
#endif
	return tfileptr[j];
	}
void
threadLogClose(int j) {
	if(j<0) stop("threadLogClose thread number < 0");
	if(j>=max_threads) stop("threadLogClose thread number too large");
	if(NULL==tfileptr[j]) stop("closing non-open temp log file");
	fclose(tfileptr[j]);
	}
void
threadLogConcat(FILE * logFile, int j) {
	int c;
	if(j<0) stop("threadLogConcat thread number < 0");
	if(j>=max_threads) stop("threadLogConcat thread number too large");
	FILE * tmpfd = fopen(tfilename[j], "r");
	if(NULL==tmpfd) stop("cannot obtain new fd to concatenate log file");
	fflush(logFile);
	while( EOF != (c=getc(tmpfd)) ) putc(c, logFile);
	fflush(logFile);
	fclose(tmpfd);
	}
void
threadLogFree(int j) {
	if(j<0) stop("threadLogFree thread number < 0");
	if(j>=max_threads) stop("threadLogFree thread number too large");
	if(NULL==tfileptr[j]) stop("freeing non-active temp log file");
#pragma omp critical
{
	unlink(tfilename[j]);
	free(tfilename[j]);
	tfileptr[j]=NULL;
}
	}
void
threadLogFreeAll(void) {
	// for emergency cleanup of the temporary files on abnormal exit
   for(int j=0;j<max_threads;j++) {
	if(NULL==tfileptr[j]) continue;
	unlink(tfilename[j]);
	//free(tfilename[j]);
	//tfileptr[j]=NULL;
   }
}
