/* prototypes.h */

#ifndef _WIN32
#include <sys/times.h>
#endif

#include <sys/types.h>
#include <stdio.h>
#include "parameters.h"


void	banner( double version_num );
int	    setflags( int argc, char **argv );
ParameterEntry * apm_find( const char key[] );
void    apm_enter( const char key[], ParameterEntry value );
int	    check_size( int nelements, char axischar );
int	    gpfparser( char line[LINE_LEN] );
int	    main( int argc, char **argv );
int	    parsetypes(char * line, char *words[], int maxwords);
void	prHMSfixed( float t );
void	printdate( FILE *fp, int flag );
void	printhms( float t );
int	    strindex( char s[], char t[] );
void	timesys( Clock duration, struct tms *start, struct tms *end );
void	timesyshms( Clock duration, struct tms *start, struct tms *end );

/* EOF */
