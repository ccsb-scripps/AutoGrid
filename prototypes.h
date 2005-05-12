/* prototypes.h */
#include <sys/times.h>
#include <sys/types.h>
#include <stdio.h>
#include "parameters.h"

int	strindex( char s[], char t[] );
void	banner( double version_num );
int	check_size( int nelements, char axischar );
int	gpfparser( char line[LINE_LEN] );
int	parsetypes(char * line, char *words[], int maxwords);
int	main( int argc, char **argv );
void	printdate( FILE *fp, int flag );
void	printhms( float t );
void	prHMSfixed( float t );
int	    setflags( int argc, char **argv );
void	timesys( Clock duration, struct tms *start, struct tms *end );
void	timesyshms( Clock duration, struct tms *start, struct tms *end );
void    apm_enter(const char key[], ParameterEntry value);
ParameterEntry * apm_find(const char key[]);

/* EOF */
