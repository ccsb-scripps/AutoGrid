/* get_atom_type.cc */



    #include <stdio.h>
    #include "get_atom_type.h"


extern char *programname;
extern FILE *logFile;


int get_atom_type( char aname[4], 
		   char chtype[MAX_TYPES] )

{
    char ch, ch1;
    register int i; 
    int type = -1;

    /* DEBUG (void) fprintf (stderr, "DEBUG: aname='%s', chtype='%s'\n", aname, chtype); */

    ch1 = aname[0];
    /* DEBUG (void) fprintf (stderr, "DEBUG: get_atom_type: ch1= '%c'\n", ch1); */

    switch (ch1) {
        case ' ':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
            ch = aname[1];
            break;
        default:
            ch = ch1;
            break;
    }
    /* ch = (ch1 == ' ') ? aname[1] : ch1; */

    /* DEBUG (void) fprintf (stderr, "DEBUG: get_atom_type: ch= '%c'\n", ch); */
    for ( i=0; i<MAX_TYPES; i++) {
        /* DEBUG (void) fprintf (stderr, "DEBUG: get_atom_type: i= %d\n", i); */
        if (ch == chtype[i]) {
            type = i;
            /* DEBUG (void) fprintf (stderr, "DEBUG: get_atom_type: type= %d\n", type); */
            break; /* break out of this for i loop */
        }
    }    
    if (type == -1) {
        (void) fprintf(stderr, "%s: Atom type error, can't find type for \"%s\" in typelist \"%s\".\n", programname, aname, chtype );
        (void) fprintf(logFile, "%s: Atom type error, can't find type for \"%s\" in typelist \"%s\".\n", programname, aname, chtype );
        type = UNKNOWN;
    }
    return (type);
}
/* EOF */
