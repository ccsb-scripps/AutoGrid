/* get_atom_type.cc */



    #include <stdio.h>
    #include "get_atom_type.h"


extern char *programname;
extern FILE *logFile;


int get_atom_type( char aname[4], 
		   char chtype[MAX_TYPES],
           int * type_ct_ptr)

{
    char ch, ch1;
    register int i; 
    int type = -1;
    int type_ct;

    /* DEBUG (void) fprintf (stderr, "DEBUG: aname='%s', chtype='%s'\n", aname, chtype); */

    type_ct = *type_ct_ptr;
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
    for ( i=0; i<type_ct; i++) {
        /* DEBUG (void) fprintf (stderr, "DEBUG: get_atom_type: i= %d\n", i); */
        if (ch == chtype[i]) {
            type = i;
            /* DEBUG (void) fprintf (stderr, "DEBUG: get_atom_type: type= %d\n", type); */
            break; /* break out of this for i loop */
        }
    }    
    if (type == -1) {
        /*append it to chtype and increment type_ct*/
        type = type_ct;
        chtype[type_ct++] = ch;
        for ( i=0; i<type_ct; i++) {
        };
        *type_ct_ptr = type_ct;
    }

    /**type_ct_ptr = type_ct;*/
    if (type == -1) {
        (void) fprintf(stderr, "\n%s: WARNING!  Atom type error, can't find type for \"%s\" in typelist \"%s\".\n", programname, aname, chtype );
        (void) fprintf(logFile, "\n%s: WARNING!  Atom type error, can't find type for \"%s\" in typelist \"%s\".\n", programname, aname, chtype );
        type = UNKNOWN;
    }
    return (type);
}
/* EOF */
