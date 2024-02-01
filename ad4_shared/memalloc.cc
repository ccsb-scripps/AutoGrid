/*
 "safe" malloc() and calloc()
 M Pique - adapted from DOT code 2019-06
*/


#include <stdlib.h>
#include <stdio.h>

#include "memalloc.h"
#include "stop.h"

void *
malloc_t(size_t bytes, const char *name)
{
void * p = malloc(bytes);

if(p==NULL) {
	char msg[1000];
	sprintf(msg,  "Unable to allocate memory for %s (%d bytes).", 
	  name, (int) bytes);
	stop(msg);
	}
return p;
}

void *
calloc_t(size_t nelem, size_t elsize, const char *name)
{
void * p = calloc(nelem, elsize);

if(p==NULL) {
	char msg[1000];
	sprintf(msg, "Unable to allocate memory for %s (%d elements of %d bytes each).",
	  name, (int) nelem, (int) elsize);
	stop(msg);
	}
return p;
}

