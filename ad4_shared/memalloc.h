/*
*/

#ifndef MEMALLOC_H


void *
malloc_t(size_t bytes, const char *name);

void *
calloc_t(size_t nelem, size_t elsize, const char *name);

#define MEMALLOC_H
#endif
