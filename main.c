/*

 $Id: main.c,v 1.2 2003/02/12 19:32:29 lindy Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* main.c */

#include <sys/types.h> 
#include <sys/times.h>
#include <sys/param.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h> /* tolower */
#include <unistd.h> /* long sysconf(int name) */

#include "autogrid.h"
#include "autoglobal.h"
#include "autocomm.h"

extern float idct;

int main( int argc,  char **argv )

/******************************************************************************/
/*      Name: main (executable's name is "autogrid").                         */
/*  Function: Calculation of interaction energy grids for Autodock.           */
/*            Directional H_bonds implemented, after Goodford:                */
/*              Oxygen accepter to probe Hd                       t0/ti model */
/*              Hydroxyl hydrogen donor to probe Oa                  cos**4   */
/*              Hydrogen on nitrogen donor to probe Oa               cos**2   */
/*              Directional nitrogen acceptor                         n/a     */
/*              Disordered hydroxyls                                 cos**4   */
/*            Distance dependent dielectric after Mehler and Solmajer.        */
/*            Solvation term of Stouten et al 1993                            */
/* Copyright: (C) 2000, TSRI                                                  */
/*                                                                            */
/*   Authors: Garrett Matthew Morris, David S. Goodsell                       */
/*                                                                            */
/*            The Scripps Research Institute                                  */
/*            Department of Molecular Biology, MB5                            */
/*            10666 North Torrey Pines Road                                   */
/*            La Jolla, CA 92037.                                             */
/*                                                                            */
/*            e-mail: garrett@scripps.edu                                     */
/*                    goodsell@scripps.edu                                    */
/*                                                                            */
/*            Helpful suggestions and advice:                                 */
/*            Arthur J. Olson                                                 */
/*            Bruce Duncan, Yng Chen, Michael Pique, Victoria Roberts         */
/*                                                                            */
/*      Date: 05/02/00                                                        */
/*                                                                            */
/*    Inputs: Control file, macromolecular PDB file                           */
/*            Note: atomic charge, volume, and solvation parameter must be    */
/*                  included in each atomic record                            */
/*   Returns: Atomic affinity and electrostatic grid maps.                    */
/*   Globals: MAX_DIST, MAX_MAPS, MAX_TYPES                                   */
/*                                                                            */
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 07/06/89 DSG     FORTRAN implementation                                    */
/* 07/05/92 GMM     C translation                                             */
/* 14/05/92 GMM     Time-dependent seed in random-number generation           */
/* 10/07/92 GMM     ASCII grid display of extrema and midpoint                */
/* 13/07/92 GMM     Barycentre calculation of centre of mass                  */
/* 14/07/92 GMM     AVS field format for displaying grids                     */
/* 27/07/92 GMM     Introduced 'atom_maps' -> variable number of maps         */
/* 16/10/92 GMM     H-bonding allowed in any map                              */
/* 17/11/92 GMM     Header information included in grid maps, for checking.   */
/* 06/11/92 GMM     Command line parsing, using Bruce S. Duncan's "setflags". */
/* 04/01/93 GMM     Created for use in makefile.                              */
/* 03/11/94 GMM     Four-fold acceleration due to non-bond cutoff checking    */
/*                  Plus correction of normalization for directional H-bonds  */
/*                  Plus floating grid calculation                            */
/*                  Plus system times                                         */
/* 02/06/95 GMM     New keyword-based interface; added r-eqm & epsilon option */
/* 29/08/95 DSG     Directional H-Bonds on oxygen atoms                       */
/*                  Disordered hydroxyls (automatically found 29/04/96)       */
/* 20/09/95 DSG     Solvation Term of Stouten et al 1993                      */
/******************************************************************************/

{
char atomname[6];
char atmtypstr[MAX_MAPS];
char AVSmaps_fn[MAX_CHARS];
char chtype[MAX_TYPES];
char error_message[LINE_LEN];
/*char extension[5];*/
char filename[MAX_MAPS][MAX_CHARS];
char fname[MAX_CHARS];
char hostnm[MAX_CHARS];
char line[LINE_LEN];
char gpfline[LINE_LEN];
char macromol_fn[MAX_CHARS];
/* char q_str[7]; */
char record[6];
char tempchar = ' ';
char token[5];
char warned = 'F';
char xyz[5];
char xyzfilename[MAX_CHARS];

FILE *macromolFile, 
     *AVSmapsFLD, 
     *xyzFile, 
     *mapFile[MAX_MAPS], 
     *disFile;

double A, epsilon0, rk, lambda, B, lambda_B;
double coord[MAX_ATOMS][XYZ], charge[MAX_ATOMS], q_tot = 0.0;
double vol[MAX_ATOMS], solpar[MAX_ATOMS];
double c[XYZ];
double cext[XYZ];
double cgridmax[XYZ];
double cgridmin[XYZ];
double cmax[XYZ], cmin[XYZ];
double csum[XYZ], cmean[XYZ];
double center[XYZ];
double constant[MAX_MAPS];
/* Cartesian-coordinate of covalent affinity well. */
double covpos[XYZ];
double d[XYZ];
double dc[XYZ];
double diel, invdielcal;
double dxA;
double dxB;
double emax[MAX_MAPS], emin[MAX_MAPS];
double enrg[MAX_MAPS], cA, cB, tmpconst;
double epsilon[MAX_DIST];
double minus_inv_two_sigma_sqd;
double percentdone;
double PI_halved;
double q_max = -BIG,  q_min = BIG;
double rA;
double rB; /* double e; */
/* Distance from current grid point to the covalent attachment point */
double rcov = 0.0;
double ri, inv_rd, rd2, r; /* re, r2, rd, */
double r_min, inv_r, racc, rdon, rsph, cos_theta, theta, tmp;
double rvector[MAX_ATOMS][XYZ], rvector2[MAX_ATOMS][XYZ], cross[XYZ];
double r_smooth = 0.;
double rdot;
double Rij, epsij;
double spacing = 0.375; /* One quarter of a C-C bond length. */
double t0, ti;
double vol_probe[MAX_MAPS], solpar_probe[MAX_MAPS], sigma;
double ln_half = 0.0;
double covhalfwidth = 1.0;
double covbarrier = 1000.0;

double version_num = 3.06;

float timeRemaining = 0.;

int atom_maps_1 = 1;
int atom_maps = 0;
int Ais2B = FALSE;
int Bis2A = FALSE;
int atomtype[MAX_ATOMS], disGrid = FALSE, dddiel = FALSE, disorder_h = FALSE;
int disorder[MAX_ATOMS];
int donecovalent = FALSE;
int elecPE;
/* int covmap; */
int from, to;
int fprintf_retval = 0;
int gpfkeyword = -1;
int iscovalent[ATOM_MAPS];
int hbond[ATOM_MAPS];
int icoord[XYZ]; /* int icenter; */
int indcom = 0;
int infld;
int length = LINE_LEN;
int MD_1 = MAX_DIST - 1;
int nbond;
int ne[XYZ], n1[XYZ], nelements[XYZ];
int nDone = 0;
int problem_wrt = FALSE;
int rexp[MAX_ATOMS], veclen = 0;
int tycount[MAX_TYPES];
int xA, xB;

register int i = 0, ii = 0, j = 0, k = 0, indx_r = 0, i_smooth = 0;
/* register int i = 0, ii = 0, j = 0, jj = 0, k = 0, indx_r = 0, i_smooth = 0;
 */
register int ia = 0, ib = 0, ic = 0, probe = -1, iat = 0, i1 = 0, i2 = 0;

static int natom;
static double energy[MAX_TYPES][MAX_DIST][MAX_MAPS];
static double sol_fn[MAX_DIST], energy_smooth[MAX_DIST];
static long clktck = 0;

Clock     job_start;
Clock     job_end;
struct tms  tms_job_start;
struct tms  tms_job_end;

Clock     grd_start;
Clock     grd_end;
struct tms  tms_grd_start;
struct tms  tms_grd_end;

/*
** Fetch clock ticks per second.
*/

if (clktck == 0) {
    if ( (clktck = sysconf(_SC_CLK_TCK)) < 0) {
        (void) fprintf( stderr, "\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n" );
        (void) fprintf( logFile, "\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n" );
        exit(911);
    } else {
        idct = 1. / (float)clktck;
    }
}

ln_half = (double) log(0.5);

/*
** Get the time at the start of the run...
*/

job_start = times( &tms_job_start );

/*
** Parse the command line...
*/

(void) setflags( argc, argv );

/*
** Initialize the atom type array, 
*/

(void) strcpy( chtype, ATOMTYPE );

(void) strcpy( xyz, "xyz" );

for (i = 0;  i < XYZ;  i++) {
   icoord[i] = 0;
}
 /* Initialize max and min coodinate bins */
for (i = 0;  i < XYZ;  i++) {
        cmax[i] = -BIG;
        cmin[i] = BIG;
        csum[i] = 0.;
        covpos[i] = 0.0;
}

PI_halved = PI/2.;

/*
 * Initialize int tycount[MAX_TYPES] array to 0
 */
for (i=0; i<MAX_TYPES; i++) {
    tycount[i] = 0;
}

/*
** Output the "AutoGrid" banner...
*/

banner( version_num );

/*
** Print the time and date when the file was created...
*/

(void) fprintf( logFile, "This file was created at:\t\t\t" );
printdate( logFile, 1 );

if (gethostname( hostnm, MAX_CHARS ) == 0) {
    (void) fprintf( logFile, "                   using:\t\t\t\"%s\"\n", hostnm );
}
/*
** Read in the GPF...
*/

while( fgets( gpfline, LINE_LEN, GPF) != NULL ) {

    gpfkeyword = gpfparser( gpfline );

    /* First switch figures out how to echo the currently inputted GPF line. */

    switch( gpfkeyword ) {
        case -1:
            (void) fprintf( logFile, "GPF> %s", gpfline );
            prStr( error_message, "%s: WARNING: Unrecognized keyword in grid parameter file.\n", programname );
            (void) fprintf( logFile, "%s", error_message );
            (void) fprintf( stderr,  "%s", error_message );
            continue;  /* while fgets gpfline... */

        case GPF_NULL:
        case GPF_COMMENT:
            (void) fprintf( logFile, "GPF> %s", gpfline );
            break;

        default:
            (void) fprintf( logFile, "GPF> %s", gpfline );
            indcom = strindex( gpfline, "#" );
            if (indcom != -1) {
                gpfline[ indcom ] = '\0'; /* Truncate str. at the comment */
            }
            break;

    } /* switch */

/******************************************************************************/

    /* Second switch interprets the current GPF line. */

    switch( gpfkeyword ) {

/******************************************************************************/
 
        case GPF_NULL:
        case GPF_COMMENT:
            break;

/******************************************************************************/
 
        case GPF_RECEPTOR:
            (void) sscanf( gpfline, "%*s %s", macromol_fn);
            (void) fprintf( logFile, "Macromolecular Input File :\t%s\n\nAtom Type Assignments:\n\n", macromol_fn);
            if ( (macromolFile = fopen(macromol_fn, "r")) == NULL ) {
                (void) fprintf( stderr, "\n%s: can't find or open macromolecule-file %s\n", programname, macromol_fn);
                (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
                (void) fprintf( logFile, "\n%s: can't find or open macromolecule-file %s\n", programname, macromol_fn);
                (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
                exit(911);
            }
            ia = 0;
            while ( (fgets(line, length, macromolFile)) != NULL ) {
                (void) sscanf(line, "%6s", record);
                if (equal(record, "ATOM", 4) || /* Protein atoms */
                    equal(record, "HETA", 4) || /* Non-standard heteroatoms */
                    equal(record, "CHAR", 4)) { /* Partial Atomic Charge - 
                                                 not a PDB record */

                    /* (void) sscanf(&line[12], "%4s", atomname); */

                    (void) strncpy( atomname, &line[12], 4);
                    /* atomname is declared as an array of 6 characters, 
                     * the PDB atom name is 4 characters (C elements 0, 1, 2 and 3)
                     * but let's ensure that the fifth character (C element 4)
                     * is a null character, which terminates the string. */
                    atomname[4] = '\0';

                    /* Output the number of this atom... */
                    (void) fprintf( logFile, "Atom no. %d,  \"%s\" ", ia+1, atomname);
                    (void) fflush( logFile );
                    
                    /* Read in this macromolecule atom's coordinates, 
                     * partial charges, and
                     * solvation parameters in PDBQS format... */
                    
                    (void) sscanf(&line[30], "%lf", &coord[ia][X]);
                    (void) sscanf(&line[38], "%lf", &coord[ia][Y]);
                    (void) sscanf(&line[46], "%lf", &coord[ia][Z]);
                    
                    /* Output the coordinates of this atom... */
                    (void) fprintf( logFile, " at (%.3lf, %.3lf, %.3lf), ", 
                                    coord[ia][X], coord[ia][Y], coord[ia][Z]);
                    (void) fflush( logFile );

                    if (oldpdbq) {
                        (void) sscanf(&line[54], "%lf", &charge[ia] );
                    } else {
                        (void) sscanf(&line[70], "%lf", &charge[ia] );
                        (void) sscanf(&line[76], "%lf", &vol[ia] );
                        (void) sscanf(&line[84], "%lf", &solpar[ia] );

                        /* convert cal/molA**3 to kcal/molA**3 */
                        solpar[ia] *= 0.001;

                    }
                    q_max = max( q_max, charge[ia] );
                    q_min = min( q_min, charge[ia] );

                    if (atomname[0] == ' ') {
                        /* truncate the first character... */
                        atomname[0] = atomname[1];
                        atomname[1] = atomname[2];
                        atomname[2] = atomname[3];
                        atomname[3] = '\0';
                    } else if (atomname[0] == '0' ||
                        atomname[0] == '1' ||
                        atomname[0] == '2' ||
                        atomname[0] == '3' ||
                        atomname[0] == '4' ||
                        atomname[0] == '5' ||
                        atomname[0] == '6' ||
                        atomname[0] == '7' ||
                        atomname[0] == '8' ||
                        atomname[0] == '9') {
                        if (atomname[1] == 'H') {
                            /* Assume this is the 'mangled' name of a hydrogen atom, 
                             * after the atom name has been changed from 'HD21' to '1HD2'
                             * for example.
                             *
                             * [0-9]H\(.\)\(.\)
                             *   0  1  2    3   
                             *   :  :  :    :
                             *   V  V  V    V
                             *  tmp 0  1    2
                             *                 tmp
                             *                  :
                             *                  V
                             *      0  1    2   3
                             *      :  :    :   :
                             *      V  V    V   V
                             *      H\(.\)\(.\)[0-9]
                             */
                            tempchar    = atomname[0];
                            atomname[0] = atomname[1];
                            atomname[1] = atomname[2];
                            atomname[2] = atomname[3];
                            atomname[3] = tempchar;
                        }
                    }
                    if (atomname[0] == 'L') {
                        (void) fprintf( logFile, "\tWARNING: lone pair atom type found; will be treated as atom type %d, \"%c\"\n", METAL+1, chtype[METAL]);
                        /* Treat Lone Pair centres as if they were a Metal...  */
                        atomtype[ia] = METAL;
                    } else {
                        /*
                         * The default atom type is 'X', 
                         * atomtype[ia] = UNKNOWN;
                         */
                        atomtype[ia] = get_atom_type(atomname, chtype);
                    }
                    /* Tell the user what you thought this atom was... */
                    (void) fprintf( logFile, "  was assigned atom type %d, \"%c\".\n", atomtype[ia]+1, chtype[atomtype[ia]]);
                    (void) fflush( logFile );
                    /* Count the number of each atom type */
                    tycount[atomtype[ia]]++;
                    /* Keep track of the extents of the macromolecule */
                    for (i = 0;  i < XYZ;  i++) {
                        cmax[i] = max( cmax[i], coord[ia][i] );
                        cmin[i] = min( cmin[i], coord[ia][i] );
                        csum[i] += coord[ia][i];
                    }
                    /* Total up the partial charges as we go... */
                    q_tot += charge[ia];
                    /* 
                     * Increment the atom counter 
                     */
                    ia++;
                    /* Check that there aren't too many atoms... */
                    if (ia > MAX_ATOMS) {
                        (void) fprintf( logFile, "Error      : Sorry, AutoGrid cannot continue.\n");
                        (void) fprintf( logFile, "           : Too many atoms in macromolecule input file %s;\n", macromol_fn);
                        (void) fprintf( logFile, "           : -- the maximum number of atoms, MAX_ATOMS, allowed is %d.\n", MAX_ATOMS);
                        (void) fprintf( logFile, "Suggestion : Increase the value in the \"#define MAX_ATOMS %d\" line", MAX_ATOMS);
                        (void) fprintf( logFile, "           : in the source file \"autogrid.h\", and re-compile AutoGrid.\n" );
                        (void) fflush( logFile );
                        exit(911);
                    } /* endif */
                } /* endif */
            } /* endwhile */
            (void) fclose( macromolFile );
            /* Update the total number of atoms in the macromolecule */
            natom = ia;
            (void) fprintf( logFile, "\nMaximum partial atomic charge found = %+.3lf e\n", q_max );
            (void) fprintf( logFile, "Minimum partial atomic charge found = %+.3lf e\n\n", q_min );
            (void) fflush( logFile );
            /* Check there are partial charges... */
            if (q_max == 0. && q_min == 0.) {
                (void) fprintf( logFile, "\nWARNING!  Partial atomic charges not found!\n\n" );
                oldpdbq = TRUE;
                (void) fprintf( logFile, "WARNING! It appears that the old PDBQ format is present. Switching \"-o\" flag on...\n");
                (void) fprintf( stderr, "\nWARNING! It appears that the old PDBQ format is present. Switching \"-o\" flag on...\n\n");
                (void) fflush( logFile );
                ia = 0;
                q_tot = 0.0;
                if ( (macromolFile = fopen(macromol_fn, "r")) == NULL ) {
                    (void) fprintf( stderr, "\n%s: For some reason, I can't re-open the macromolecule file %s\n", programname, macromol_fn);
                    (void) fprintf( stderr, "\n%s: For some reason, I can't re-open the macromolecule file %s\n", programname, macromol_fn);
                    (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
                    (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
                    exit(911);
                }
                while ( (fgets(line, length, macromolFile)) != NULL ) {
                    (void) sscanf(line, "%6s", record);
                    if (equal(record, "ATOM", 4) ||         /* Atom - protein atoms */
                        equal(record, "HETA", 4) ||        /* Heteroatom - non-protein atoms */
                        equal(record, "CHAR", 4)) {        /* Charge - non-PDB-standard record */
                        (void) sscanf(&line[54], "%lf", &charge[ia] );
                        q_max = max( q_max, charge[ia] );
                        q_min = min( q_min, charge[ia] );
                        q_tot += charge[ia];
                        ia++;
                    } /* fi */
                } /* while */
                (void) fclose( macromolFile );
                (void) fprintf( logFile, "Maximum partial atomic charge found = %+.3lf e\n", q_max );
                (void) fprintf( logFile, "Minimum partial atomic charge found = %+.3lf e\n\n", q_min );
                (void) fflush( logFile );
            } /* if */
            for (ia = 0;  ia < natom;  ia++) {
                rexp[ia] = 0;
            }
            (void) fprintf( logFile, "Atom\tAtom\tNumber of this Type\n");
            (void) fprintf( logFile, "Type\t ID \t in Macromolecule\n");
            (void) fprintf( logFile, "____\t____\t___________________\n");
            (void) fflush( logFile );
            for (i = 0;  i < NATOMTYPES;  i++) {
                (void) fprintf( logFile, " %d\t %c\t\t%6d\n", (i+1), chtype[i], tycount[i]);
            }
            (void) fprintf( logFile, "\nTotal number of atoms :\t\t%d atoms \n", natom);
            (void) fflush( logFile );
            (void) fprintf( logFile, "Total charge :\t\t\t%.2lf e\n", q_tot);
            (void) fflush( logFile );
            (void) fprintf( logFile, "\n\nMacromolecule coordinates fit within the following volume:\n\n");
            (void) fflush( logFile );
            (void) fprintf( logFile, "                   _______(%.1lf, %.1lf, %.1lf)\n", cmax[X], cmax[Y], cmax[Z]);
            (void) fprintf( logFile, "                  /|     /|\n");
            (void) fprintf( logFile, "                 / |    / |\n");
            (void) fprintf( logFile, "                /______/  |\n");
            (void) fprintf( logFile, "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n", (cmax[X]+cmin[X])/2., (cmax[Y]+cmin[Y])/2., (cmax[Z]+cmin[Z])/2.);
            (void) fprintf( logFile, "                |  /   |  /\n");
            (void) fprintf( logFile, "                | /    | /\n");
            (void) fprintf( logFile, "                |/_____|/\n");
            (void) fprintf( logFile, "(%.1lf, %.1lf, %.1lf)      \n", cmin[X], cmin[Y], cmin[Z]);
            (void) fprintf( logFile, "\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cmax[X], cmax[Y], cmax[Z]);
            (void) fprintf( logFile, "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cmin[X], cmin[Y], cmin[Z]);
            (void) fprintf( logFile, "\n");
            for (i = 0;  i < XYZ;  i++) {
                cmean[i] = csum[i] / (double)natom;
            }
            (void) fflush( logFile );
            break;

/******************************************************************************/
 
        case GPF_GRIDFLD:
            (void) sscanf( gpfline, "%*s %s", AVSmaps_fn);
            infld = strindex( AVSmaps_fn, ".fld" );
            if (infld == -1) {
                (void) fprintf( stderr, "\n\n%s: WARNING:  Grid data file needs the extension \".fld\" for AVS input\n\n", programname);
                (void) fprintf( logFile, "\n\n%s: WARNING:  Grid data file needs the extension \".fld\" for AVS input\n\n", programname);
                exit(911);
            } else {
                infld = strindex( AVSmaps_fn, "fld" );
                (void) strcpy(xyzfilename, AVSmaps_fn);
                xyzfilename[infld] = 'x';
                xyzfilename[infld + 1] = 'y';
                xyzfilename[infld + 2] = 'z';
            }
            if ( (AVSmapsFLD = fopen(AVSmaps_fn, "w")) == NULL ) {
                (void) fprintf( stderr, "\n%s: can't create grid dimensions data file %s\n", programname, AVSmaps_fn);
                (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
                (void) fprintf( logFile, "\n%s: can't create grid dimensions data file %s\n", programname, AVSmaps_fn);
                (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
                exit(911);
            } else {
                (void) fprintf( logFile, "\nCreating (AVS-readable) grid maps file : %s\n", AVSmaps_fn);
            }
            if ( (xyzFile = fopen(xyzfilename, "w")) == NULL ) {
                (void) fprintf( stderr, "\n%s: can't create grid extrema data file %s\n", programname, xyzfilename);
                (void) fprintf( stderr, "%s: SORRY!    unable to create the \".xyz\" file.\n\n", programname);
                (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
                (void) fprintf( logFile, "\n%s: can't create grid extrema data file %s\n", programname, xyzfilename);
                (void) fprintf( logFile, "%s: SORRY!    unable to create the \".xyz\" file.\n\n", programname);
                (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
                exit(911);
            } else {
                (void) fprintf( logFile, "\nCreating (AVS-readable) grid-coordinates extrema file : %s\n\n", xyzfilename);
            }
            break;

/******************************************************************************/

        case GPF_NPTS:
            (void) sscanf( gpfline, "%*s %d %d %d", &nelements[X], &nelements[Y], &nelements[Z] );
            for (i = 0;  i < XYZ;  i++) {
                nelements[i] = check_size(nelements[i], xyz[i]);
                ne[i] = nelements[i] / 2;
                n1[i] = nelements[i] + 1;
            }
            (void) fprintf( logFile, "\n");
            (void) fprintf( logFile, "Number of grid points in x-direction:\t%d\n", n1[X]);
            (void) fprintf( logFile, "Number of grid points in y-direction:\t%d\n", n1[Y]);
            (void) fprintf( logFile, "Number of grid points in z-direction:\t%d\n", n1[Z]);
            (void) fprintf( logFile, "\n");
            percentdone = 100. / (double) n1[Z];
            break;

/******************************************************************************/

        case GPF_SPACING:
            (void) sscanf( gpfline, "%*s %lf", &spacing );
            (void) fprintf( logFile, "Grid Spacing :\t\t\t%.3lf Angstrom\n", spacing);
            (void) fprintf( logFile, "\n");
            break;

/******************************************************************************/

        case GPF_GRIDCENTER:
            (void) sscanf( gpfline, "%*s %s", token);
            if (equal( token, "auto", 4)) {
                for (i = 0;  i < XYZ;  i++) {
                    center[i] = cmean[i];
                }
                (void) fprintf( logFile, "Grid maps will be centered on the center of mass.\n");
                (void) fprintf( logFile, "Coordinates of center of mass : (%.3lf, %.3lf, %.3lf)\n", center[X], center[Y], center[Z]);
            } else {
                (void) sscanf( gpfline, "%*s %lf %lf %lf", &center[X], &center[Y], &center[Z]);
                (void) fprintf( logFile, "\nGrid maps will be centered on user-defined coordinates:\n\n\t\t(%.3lf, %.3lf, %.3lf)\n",  center[X], center[Y], center[Z]);
            }
            /* centering stuff... */
            for (ia = 0;  ia < natom;  ia++) {
                for (i = 0;  i < XYZ;  i++) {
                    coord[ia][i] -= center[i];        /* transform to center of gridmaps */
                }
            }
            for (i = 0;  i < XYZ;  i++) {
                cext[i]     = spacing * (double)ne[i];
                cgridmax[i] = center[i] + cext[i];
                cgridmin[i] = center[i] - cext[i];
            }
            (void) fprintf( logFile, "\nGrid maps will cover the following volume:\n\n");
            (void) fprintf( logFile, "                   _______(%.1lf, %.1lf, %.1lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z] );
            (void) fprintf( logFile, "                  /|     /|\n");
            (void) fprintf( logFile, "                 / |    / |\n");
            (void) fprintf( logFile, "                /______/  |\n");
            (void) fprintf( logFile, "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n", center[X], center[Y], center[Z]);
            (void) fprintf( logFile, "                |  /   |  /\n");
            (void) fprintf( logFile, "                | /    | /\n");
            (void) fprintf( logFile, "                |/_____|/\n");
            (void) fprintf( logFile, "(%.1lf, %.1lf, %.1lf)      \n\n", cgridmin[X], cgridmin[Y], cgridmin[Z] );
            for (i = 0;  i < XYZ;  i++) {
                (void) fprintf( logFile, "Grid map %c-dimension :\t\t%.1lf Angstroms\n", xyz[i], 2.*cext[i] );
            }
            (void) fprintf( logFile, "\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z] );
            (void) fprintf( logFile, "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cgridmin[X], cgridmin[Y], cgridmin[Z] );
            for (i = 0;  i < XYZ;  i++) {
                (void) fprintf(xyzFile, "%.3lf %.3lf\n", cgridmin[i], cgridmax[i]);
            }
            (void) fclose(xyzFile);
            break;

/******************************************************************************/

        case GPF_TYPES:
            (void) sscanf( gpfline, "%*s %s", atmtypstr);
            elecPE = atom_maps = strlen(atmtypstr);
            veclen = atom_maps_1 = atom_maps + 1;
            for (i = 0;  i < atom_maps;  i++) {
                iscovalent[i] = FALSE;
                hbond[i] = FALSE;
            }
            (void) fprintf( logFile, "\nAtom names for types 1-%d and for substrate-atom affinity grid maps :\t", atom_maps);
            for (i = 0;  i < atom_maps;  i++) {
                (void) fprintf( logFile, "%c ", atmtypstr[i]);
                if (atmtypstr[i] == COVALENTTYPE) {
                  iscovalent[i] = TRUE;
                  (void) fprintf( logFile, "\nAtom type number %d will be used to calculate a covalent affinity grid map\n\n", i+1);
                }
            }
            (void) fprintf( logFile, " respectively.\n");
            break;

/******************************************************************************/
 
        case GPF_SOL_PAR:
            /*
            ** read volume and solvation parameter for probe *****************
            */
            if (probe == -1) {
                (void) fprintf( logFile, "ERROR!  You must specify the \"map filename\" and all the \"nbp_r_eps\" parameters before using this command.\n\n" );
                exit(-1);
            }
            (void) sscanf( gpfline, "%*s %lf %lf", &vol_probe[probe], &solpar_probe[probe]);
            (void) fprintf( logFile, "\nProbe solvation parameters: \n\n\tatomic fragmental volume: %.2f A^3\n\tatomic solvation parameter: %.4f cal/mol A^3\n\n", vol_probe[probe], solpar_probe[probe]);
            (void) fflush( logFile ); 
            solpar_probe[probe] *= 0.001; /*convert cal/molA^3 to kcal/molA^3 */
            break; /* end solvation parameter */

/******************************************************************************/

        case GPF_CONSTANT:
            /*
            ** read constant term for probe **********************************
            */
            if (probe == -1) {
                (void) fprintf( logFile, "ERROR!  You must specify the \"map filename\" and all the \"nbp_r_eps\" parameters before using this command.\n\n" );
                exit(-1);
            }
            (void) sscanf( gpfline, "%*s %lf", &constant[probe]);
            (void) fprintf( logFile, "\n\nConstant added to map: %.3f\n\n", constant[probe]);
            (void) fflush( logFile ); 
            break;

/******************************************************************************/

        case GPF_MAP:
            /* probe is the index of the atom type we are calculating a map
               for. */
            ++probe;
            (void) sscanf( gpfline, "%*s %s", filename[ probe ] );
            if ( (mapFile[ probe ] = fopen( filename[ probe ], "w")) == NULL ) {
                (void) fprintf( stderr, "\n%s: can't open grid map \"%s\" for writing.\n", programname, filename[probe] );
                (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
                (void) fprintf( logFile, "\n%s: can't open grid map \"%s\" for writing.\n", programname, filename[probe] );
                (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
                exit(911);
            }
            (void) fprintf( logFile, "\n\nOutput File %d   %s\n\n", (probe+1), filename[probe] );
            (void) fflush( logFile ); 

            /* i is the index of the protein atom type, that probe will
               interact with. */
            for (i = 0;  i < NATOMTYPES;  i++) {
                if ( fgets( line, LINE_LEN, GPF) != NULL ) {
                (void) fprintf( logFile, "GPF> %s", line );
                    switch( gpfparser( line ) ) {
                        case -1:
                            prStr( error_message, "%s: WARNING: Unrecognized keyword -- use either \"nbp_r_eps\" or \"nbp_coeffs\".\n", programname );
                            (void) fprintf( logFile, "%s", error_message );
                            (void) fprintf( stderr,  "%s", error_message );
                            break;
                        case GPF_NBP_COEFFS:
                            (void) sscanf( line, "%*s %lf %lf %d %d", &cA, &cB, &xA, &xB );
                            break;
                        case GPF_NBP_R_EPS:
                            (void) sscanf( line, "%*s %lf %lf %d %d", &Rij, &epsij, &xA, &xB);
                            cA = (tmpconst = epsij / (double)(xA - xB)) * pow( Rij, (double)xA ) * (double)xB;
                            cB = tmpconst * pow( Rij, (double)xB ) * (double)xA;
                            break;
                    } /* switch */
                    dxA = (double) xA;
                    dxB = (double) xB;
                    if (xB == (2*xA)) {
                        Bis2A = TRUE;
                        Ais2B = FALSE;
                    } else if (xA == (2*xB)) {
                        Ais2B = TRUE;
                        Bis2A = FALSE;
                    } else {
                        Ais2B = FALSE;
                        Bis2A = FALSE;
                    }
                    if ( (xA != 12) || (xB != 6) ) {
                        hbond[probe] = TRUE;
                    }
                    (void) fprintf( logFile, "\n             %9.1lf       %9.1lf \n", cA, cB);
                    (void) fprintf( logFile, "    E    =  -----------  -  -----------\n");
                    (void) fprintf( logFile, "     %c, %c         %2d              %2d\n", atmtypstr[probe], chtype[i], xA, xB);
                    (void) fprintf( logFile, "                r               r \n\n");
                    /* loop over distance index, indx_r, from 0 to MAX_DIST */
                    for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
                        r  = angstrom(indx_r);
                        if (Bis2A) {
                            rA = pow( r, dxA );
                            rB = rA * rA;
                        } else if (Ais2B) {
                            rB = pow( r, dxB );
                            rA = rB * rB;
                        } else {
                            rA = pow( r, dxA );
                            rB = pow( r, dxB );
                        } /* if */
                        /* calculate the dispersion-repulsion or hydrogen bonding energy
                           between protein atom type i and ligand atom type
                           probe, and store it for this distance, indx_r. */
                        energy[i][indx_r][probe] = min( EINTCLAMP, (cA/rA - cB/rB) );
                    } /* indx_r */
                    energy[i][0][probe]    = EINTCLAMP;
                    energy[i][MD_1][probe] = 0.;

                    /* smooth with min function */
                    if (i_smooth > 0) {
                        /* loop over distance index, indx_r, from 0 to MAX_DIST */
                        for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
                            /* set energy_smooth at this distance, indx_r, to */
                            /* 100, 000 -- effectively clamping the maximum */
                            /* energy. */
                            energy_smooth[indx_r] = 100000.;
                            /* loop over j, from either 0 or indx_r - i_smooth, whichever is larger,  */
                            /*              to either MAX_DIST or indx_r + i_smooth, whichever is smaller. */
                            /* This effectively goes over all distances from */
                            /* indx_r, plus or minus i_smooth. */
                            for (j = max(0, indx_r-i_smooth);  j < min(MAX_DIST, indx_r+i_smooth);  j++) {
                              /* store the minimum energy within the this */
                              /* range, of this distance, indx_r, plus or minus the smoothing radius, i_smooth */
                              energy_smooth[indx_r] = min(energy_smooth[indx_r], energy[i][j][probe]);
                            }
                        }
                        /* loop over distance index, indx_r, from 0 to MAX_DIST */
                        for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
                            /* update the interatomic energy table with the */
                            /* smoothed value at this distance, indx_r. */
                            energy[i][indx_r][probe] = energy_smooth[indx_r];
                        }
                    } /* end smoothing */
                } /* if */
            } /*  i  */

           /* exponential function for protein and ligand desolvation */
           /* note: the solvation term will not be smoothed */
           sigma = 3.6;
           minus_inv_two_sigma_sqd = -1. / (2. * sigma * sigma);
           for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
                r  = angstrom(indx_r);
                /* sol_fn[indx_r] = exp(-sq(r)/(2.*sigma*sigma)); */
                sol_fn[indx_r] = exp( sq(r) * minus_inv_two_sigma_sqd );
           }
/*
** Print out a table, of distance versus energy...
*/
            (void) fprintf( logFile, "\n  r ");
            for (iat = 0;  iat < NATOMTYPES;  iat++) {
                (void) fprintf( logFile, "    %c    ", chtype[iat]);
            } /* iat */
            (void) fprintf( logFile, "\n ___");
            for (iat = 0;  iat < NATOMTYPES;  iat++) {
                (void) fprintf( logFile, " ________");
            } /* iat */
            (void) fprintf( logFile, "\n");
            for (i = 0;  i <= 500;  i += 10) {
                (void) fprintf( logFile, "%4.1lf", angstrom(i) );
                for (iat = 0;  iat < NATOMTYPES;  iat++) {
                    (void) fprintf( logFile, (energy[iat][i][probe]<100000.)?"%9.2lf":"%9.2lg", energy[iat][i][probe]);
                } /* iat */
                (void) fprintf( logFile, "\n");
            } /* i */
            (void) fprintf( logFile, "\n");
            break;

/******************************************************************************/
 
        case GPF_ELECMAP:
            (void) sscanf( gpfline, "%*s %s", filename[ elecPE ] );
            if ( (mapFile[ elecPE ] = fopen( filename[ elecPE ], "w" )) == NULL){
                (void) fprintf( stderr, "\n%s: can't open grid map \"%s\" for writing.\n", programname, filename[elecPE] );
                (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
                (void) fprintf( logFile, "\n%s: can't open grid map \"%s\" for writing.\n", programname, filename[elecPE] );
                (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
                exit(911);
            }
            (void) fprintf( logFile, "\nElectrostatics File: %s\n\n", filename[ elecPE ]);
            break;

/******************************************************************************/

        case GPF_COVALENTMAP:
            (void) sscanf( gpfline, "%*s %lf %lf %lf %lf %lf", &covhalfwidth, &covbarrier, &(covpos[X]), &(covpos[Y]), &(covpos[Z]) );
            (void) fprintf( logFile, "\ncovalentmap <half-width in Angstroms> <barrier> <x> <y> <z>\n");
            (void) fprintf( logFile, "\nCovalent well's half-width in Angstroms:         %8.3f\n", covhalfwidth);
            (void) fprintf( logFile, "\nCovalent barrier energy in kcal/mol:             %8.3f\n", covbarrier);
            (void) fprintf( logFile, "\nCovalent attachment point will be positioned at: (%8.3f, %8.3f, %8.3f)\n\n", covpos[X], covpos[Y], covpos[Z]);
            for (i = 0;  i < XYZ;  i++) {
                /* center covpos in the grid maps frame of reference, */
                covpos[i] -= center[i];
            }
            break;

/******************************************************************************/
 
        case GPF_DISORDER:
            disorder_h = TRUE;
            (void) fprintf( logFile, "\nHydroxyls will be disordered \n\n");
            break;

/******************************************************************************/
 
        case GPF_SMOOTH:
            (void) sscanf( gpfline, "%*s %lf", &r_smooth );
            (void) fprintf( logFile, "\nPotentials will be smoothed by: %lf\n\n", r_smooth);
            /* Angstrom is divided by A_DIVISOR in look-up table. */
            /* Typical value of r_smooth is 0.5 Angstroms */
            /* so i_smooth = 0.5 * 100. / 2 = 25 */
            i_smooth = r_smooth*A_DIVISOR/2;
            break;

/******************************************************************************/
 
        case GPF_DIEL:
            (void) sscanf( gpfline, "%*s %lf", &diel );
            if (diel < 0.) {
                /* negative... */
                (void) fprintf( logFile, "\nUsing *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.\n\n");
                dddiel = TRUE;
/*____________________________________________________________________________
** Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
**____________________________________________________________________________*/
                lambda   = 0.003627;
                epsilon0 = 78.4;
                A        = (-8.5525);
                B        = epsilon0 - A;
                rk       = 7.7839;
                lambda_B = (-lambda) * B;
                for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
                    epsilon[indx_r] = A + B / (1. + rk*exp(lambda_B*angstrom(indx_r)) );
                }
                if (epsilon[0] < APPROX_ZERO) {
                    epsilon[0] = 1.;
                }
                (void) fprintf( logFile, "  d   Dielectric\n ___  __________\n");
                for (i = 0;  i <= 500;  i += 10) {
                    ri = angstrom(i);
                    (void) fprintf( logFile, "%4.1lf%9.2lf\n", ri, epsilon[i]);
                }
                (void) fprintf( logFile, "\n");
                for (i = 1;  i < MAX_DIST;  i++) {
                    epsilon[i] = 332.0 / epsilon[i]; 
                   /* -diel: weight from free energy survey */
                    epsilon[i] *= -diel;
    /* Really, LHS should be "inv_epsilon[i]", but this way saves memory... */
                }
            } else {
                /* positive or zero... */
                dddiel = FALSE;
                if (diel <= APPROX_ZERO) {
                    diel = 40.;
                }
                (void) fprintf( logFile, "Using a *constant* dielectric of:  %.2f\n", diel );
                invdielcal = 332. / diel;
            }
            break;

/******************************************************************************/
 
        case GPF_FMAP:
            (void) sscanf( gpfline, "%*s %s", fname );
            if ( (disFile = fopen( fname, "w" )) == NULL) {
                (void) fprintf( stderr, "\n%s: can't open grid map \"%s\" for writing.\n", programname, fname );
                (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
                (void) fprintf( logFile, "\n%s: can't open grid map \"%s\" for writing.\n", programname, fname );
                (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
                exit(911);
            }
            (void) fprintf( logFile, "\nFloating Grid file name = %s\n", fname );
            ++veclen;
            disGrid = TRUE;
            break;

/******************************************************************************/
 
        default:
            break;

/******************************************************************************/

    } /* switch */ 
} /* while */

(void) fprintf( logFile, "\n>>> Closing the grid parameter file (GPF)... <<<\n\n" );
(void) fprintf( logFile, UnderLine );
(void) fclose( GPF );

if ( ! disGrid ) {
    (void) fprintf( logFile, "\n\nNo Floating Grid was requested.\n" );
}

(void) fprintf( AVSmapsFLD, "# AVS field file\n#\n");
(void) fprintf( AVSmapsFLD, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
(void) fprintf( AVSmapsFLD, "# Created by %s.\n#\n", programname);
(void) fprintf( AVSmapsFLD, "#SPACING %.3f\n", (float) spacing);
(void) fprintf( AVSmapsFLD, "#NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
(void) fprintf( AVSmapsFLD, "#CENTER %.3f %.3f %.3f\n", (float)center[X], (float)center[Y], (float)center[Z]);
(void) fprintf( AVSmapsFLD, "#MACROMOLECULE %s\n", macromol_fn);
(void) fprintf( AVSmapsFLD, "#GRID_PARAMETER_FILE %s\n#\n", grid_param_fn);
(void) fprintf( AVSmapsFLD, "ndim=3\t\t\t# number of dimensions in the field\n");
(void) fprintf( AVSmapsFLD, "dim1=%d\t\t\t# number of x-elements\n", n1[X]);
(void) fprintf( AVSmapsFLD, "dim2=%d\t\t\t# number of y-elements\n", n1[Y]);
(void) fprintf( AVSmapsFLD, "dim3=%d\t\t\t# number of z-elements\n", n1[Z]);
(void) fprintf( AVSmapsFLD, "nspace=3\t\t# number of physical coordinates per point\n");
(void) fprintf( AVSmapsFLD, "veclen=%d\t\t# number of affinity values at each point\n", veclen );
(void) fprintf( AVSmapsFLD, "data=float\t\t# data type (byte, integer, float, double)\n");
(void) fprintf( AVSmapsFLD, "field=uniform\t\t# field type (uniform, rectilinear, irregular)\n");
for (i = 0;  i < XYZ;  i++) {
    (void) fprintf( AVSmapsFLD, "coord %d file=%s filetype=ascii offset=%d\n", (i+1), xyzfilename, (i*2) );
}
for (i = 0;  i < atom_maps;  i++) {
    (void) fprintf( AVSmapsFLD, "label=%c-affinity\t# component label for variable %d\n", atmtypstr[i], (i+1));
} /* i */
(void) fprintf( AVSmapsFLD, "label=Electrostatics\t# component label for variable %d\n", veclen-1 );
if (disGrid) {
    (void) fprintf( AVSmapsFLD, "label=Floating_Grid\t# component label for variable %d\n", veclen );
}
(void) fprintf( AVSmapsFLD, "#\n# location of affinity grid files and how to read them\n#\n");
for (i = 0;  i < atom_maps;  i++) {
    (void) fprintf( AVSmapsFLD, "variable %d file=%s filetype=ascii skip=6\n", (i+1), filename[i] );
}
(void) fprintf( AVSmapsFLD, "variable %d file=%s filetype=ascii skip=6\n", atom_maps_1, filename[ elecPE ]); 
if (disGrid) {
    (void) fprintf( AVSmapsFLD, "variable %d file=%s filetype=ascii skip=6\n", veclen, fname);
}
(void) fclose( AVSmapsFLD );


/**************************************************
** Loop over all MACROMOLECULE atoms to
** calculate bond vectors for directional H-bonds
***************************************************/

for (ia = 0;  ia < natom;  ia++) {

disorder[ia] = FALSE;  /*initialize disorder flag*/
warned = 'F';

/*
** Set scan limits looking for bonded atoms
*/
    from = max(ia-20, 0);
    to   = min(ia+20, natom-1);

/*
** If 'ia' is a hydrogen atom, it could be a
** MACROMOLECULE HYDROGEN-BOND DONOR, 
*/
    if (atomtype[ia] == HYDROGEN) {

        for ( ib = from; ib <= to; ib++) {
/*
** =>  NH-> or OH->
*/
            if ((atomtype[ib] == NITROGEN) || (atomtype[ib] == OXYGEN)) {
/*
** Calculate the square of the N-H or O-H bond distance, rd2, 
**                            ib-ia  ib-ia
*/
                for (i = 0;  i < XYZ;  i++) {
                    d[i] = coord[ia][i] - coord[ib][i];
                }
                rd2 = sq( d[X] ) + sq( d[Y] ) + sq( d[Z] );
/*
** If ia & ib are less than 1.3 A apart -- they are covalently bonded, 
*/
                if (rd2 < 1.69) {
                    if (rd2 < APPROX_ZERO) {
                        if (rd2 == 0.) {
                            (void) fprintf (stderr, "WARNING! While calculating an H-O or H-N bond vector...\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                            (void) fprintf (logFile, "WARNING! While calculating an H-O or H-N bond vector...\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1./sqrt(rd2);
/*
** N-H: Set exponent rexp to 2 for m/m H-atom, 
*/
                    if (atomtype[ib] == NITROGEN) rexp[ia] = 2;
/*
** O-H: Set exponent rexp to 4 for m/m H-atom, 
** and flag disordered hydroxyls
*/
                    if (atomtype[ib] == OXYGEN) {
                        rexp[ia] = 4;
                        if (disorder_h == TRUE) disorder[ia] = TRUE;
                    }
/*
** Normalize the vector from ib to ia, N->H or O->H...
*/
                    for (i = 0;  i < XYZ;  i++) {
                        rvector[ia][i] = d[i] * inv_rd;
                    }
/*
** First O-H/N-H H-bond-donor found; Go on to next atom, 
*/
                    break;
                } /* Found covalent bond.                                 */
            } /* Found NH or OH in macromolecule.                         */
        } /* Finished scanning for the NH or OH in macromolecule.         */

/*
** If 'ia' is an Oxygen atom, it could be a
** MACROMOLECULE H_BOND ACCEPTOR, 
*/

    } else if (atomtype[ia] == OXYGEN) {
/*
** Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
**        to (ia+5)th m/m-atom
** determine number of atoms bonded to the oxygen
*/
        nbond = 0;
        for ( ib = from; ib <= to; ib++) {
            if ( ib != ia ) {
                rd2 = 0.;

                for (i = 0;  i < XYZ;  i++) {
                    dc[i] = coord[ia][i] - coord[ib][i];
                    rd2 += sq( dc[i] );
                }

                /*
                    for (i = 0;  i < XYZ;  i++) {
                        rd2 += sq(coord[ia][i] - coord[ib][i]);
                    }
                */
                if (((rd2 < 2.89) && (atomtype[ib] == CARBON)) || 
                    ((rd2 < 1.69) && (atomtype[ib] == HYDROGEN))) {
                    if (nbond == 2) {
                        (void) fprintf( logFile, "WARNING! oxygen with three bonded atoms, atom serial number %d\n", ia+1);
                    }
                    if (nbond == 1) {
                        nbond = 2;
                        i2 = ib;
                    }
                    if (nbond == 0) {
                        nbond = 1;
                        i1 = ib;
                    }
                }
            } /* ( ib != ia ) */
        } /*ib-loop*/

        /* if no bonds, something is wrong */

        if (nbond == 0) {
            (void) fprintf( logFile, "WARNING! oxygen with no bonded atoms, atom serial number %d\n", ia+1);
        }

        /* one bond: Carbonyl Oxygen O=C-X */

        if (nbond == 1) {

            /* calculate normalized carbonyl bond vector rvector[ia][] */

            rd2 = 0.;
            for (i = 0;  i < XYZ;  i++) {
                rvector[ia][i] = coord[ia][i]-coord[i1][i];
                rd2 += sq(rvector[ia][i]);
            }
            if (rd2 < APPROX_ZERO) {
                if ((rd2 == 0.) && (warned == 'F')) {
                    (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                    (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                    warned = 'T';
                }
                rd2 = APPROX_ZERO;
            }
            inv_rd = 1./sqrt(rd2);
            for (i = 0;  i < XYZ;  i++) {
                rvector[ia][i] *= inv_rd;
            }

            /* find a second atom (i2) bonded to carbonyl carbon (i1) */
            for ( i2 = from; i2 <= to; i2++) {
                if (( i2 != i1 ) && ( i2 != ia )) {
                    rd2 = 0.;
                    for (i = 0;  i < XYZ;  i++) {
                        dc[i] = coord[i1][i] - coord[i2][i]; /*NEW*/
                        rd2 += sq( dc[i] );
                    }
                    if (((rd2 < 2.89) && (atomtype[i2] != HYDROGEN)) || 
                        ((rd2 < 1.69) && (atomtype[i2] == HYDROGEN))) {

                        /* found one */
                        /* d[i] vector from carbon to second atom */
                        rd2 = 0.;
                        for (i = 0;  i < XYZ;  i++) {
                            d[i] = coord[i2][i]-coord[i1][i]; 
                            rd2 += sq( d[i] );
                        }
                        if (rd2 < APPROX_ZERO) {
                            if ((rd2 == 0.) && (warned == 'F')) {
                                (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                                (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                                warned = 'T';
                            }
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1./sqrt(rd2);
                        for (i = 0;  i < XYZ;  i++) {
                            d[i] *= inv_rd;
                        }

                        /* C=O cross C-X gives the lone pair plane normal */
                        rvector2[ia][0] = rvector[ia][1]*d[2]-rvector[ia][2]*d[1];
                        rvector2[ia][1] = rvector[ia][2]*d[0]-rvector[ia][0]*d[2];
                        rvector2[ia][2] = rvector[ia][0]*d[1]-rvector[ia][1]*d[0];
                        rd2 = 0.;
                        for (i = 0;  i < XYZ;  i++) {
                            rd2 += sq(rvector2[ia][i]);
                        }
                        if (rd2 < APPROX_ZERO) {
                            if ((rd2 == 0.) && (warned == 'F')) {
                                (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                                (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                                warned = 'T';
                            }
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1./sqrt(rd2);
                        for (i = 0;  i < XYZ;  i++) {
                            rvector2[ia][i] *= inv_rd;
                        }
                    }
                }
            }/*i2-loop*/
        } /* endif nbond==1 */

        /* two bonds: Hydroxyl or Ether Oxygen X1-O-X2 */
        if (nbond == 2) {

            /* disordered hydroxyl */

            if ( ((atomtype[i1] == HYDROGEN) || (atomtype[i2] == HYDROGEN)) 
                && (atomtype[i1] != atomtype[i2]) && (disorder_h == TRUE) )  {

                if (atomtype[i1] == CARBON) ib = i1;
                if (atomtype[i2] == CARBON) ib = i2;
                disorder[ia] = TRUE;
                rd2 = 0.;
                for (i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i]-coord[ib][i];
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                        (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1./sqrt(rd2);
                for (i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] *= inv_rd;
                }

            } else {

                /* not a disordered hydroxyl */
                /* normalized X1 to X2 vector, defines lone pair plane */

                rd2 = 0.;
                for (i = 0;  i < XYZ;  i++) {
                    rvector2[ia][i] = coord[i2][i]-coord[i1][i];
                    rd2 += sq(rvector2[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                        (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1./sqrt(rd2);
                for (i = 0;  i < XYZ;  i++) {
                    rvector2[ia][i] *= inv_rd;
                }

                /* vector pointing between the lone pairs:
                ** front of the vector is the oxygen atom, 
                ** X1->O vector dotted with normalized X1->X2 vector plus 
                ** coords of X1 gives the point on the X1-X2 line for the 
                ** back of the vector.
                */
                for (i = 0;  i < XYZ;  i++) {
                    rdot= (coord[ia][i]-coord[i1][i]) * rvector2[ia][i] ;
                }
                rd2 = 0.;
                for (i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i] - ( (rdot*rvector2[ia][i]) + coord[i1][i] ) ;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                        (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1./sqrt(rd2);
                for (i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] *= inv_rd;
                }

            } /* end disordered hydroxyl */

        }  /* end two bonds to Oxygen */

    } /* end test for atom type */

} /* Do Next macromolecular atom... */

/********************************************
** End bond vector loop
*********************************************/

for (k = 0;  k < atom_maps_1;  k++) {
    emax[k] = (double)-BIG;
    emin[k] = (double)BIG; 
}

(void) fprintf( logFile, "Beginning grid calculations.\n");
(void) fprintf( logFile, "\nCalculating %d grids over %d elements, around %d macromolecule atoms.\n\n", veclen, (n1[X]*n1[Y]*n1[Z]), natom );
(void) fflush( logFile);

/*____________________________________________________________________________
** Write out the  correct grid_data '.fld' file_name at the  head of each map
** file, to avoid centering errors in subsequent dockings...
** AutoDock can then  check to see  if the  center of each  map  matches that
** specified in its parameter file...
**____________________________________________________________________________*/

for (k = 0;  k < atom_maps_1;  k++) {
    (void) fprintf( mapFile[k], "GRID_PARAMETER_FILE %s\n", grid_param_fn);
    (void) fprintf( mapFile[k], "GRID_DATA_FILE %s\n", AVSmaps_fn);
    (void) fprintf( mapFile[k], "MACROMOLECULE %s\n", macromol_fn);
    (void) fprintf( mapFile[k], "SPACING %.3lf\n", spacing);
    (void) fprintf( mapFile[k], "NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
    (void) fprintf( mapFile[k], "CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
}
if (disGrid) {
    (void) fprintf( disFile, "GRID_PARAMETER_FILE %s\n", grid_param_fn);
    (void) fprintf( disFile, "GRID_DATA_FILE %s\n", AVSmaps_fn);
    (void) fprintf( disFile, "MACROMOLECULE %s\n", macromol_fn);
    (void) fprintf( disFile, "SPACING %.3lf\n", spacing);
    (void) fprintf( disFile, "NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
    (void) fprintf( disFile, "CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
}

(void) fprintf( logFile, "                    Percent   Estimated Time  Time/this plane\n");
(void) fprintf( logFile, "XY-plane  Z-coord   Done      Remaining       Real, User, System\n");
(void) fprintf( logFile, "            /Ang              /sec            /sec\n");
(void) fprintf( logFile, "________  ________  ________  ______________  __________________________\n\n");

/*
** Iterate over all grid points
*/

ic = 0;

for (icoord[Z] = -ne[Z]; icoord[Z] <= ne[Z]; icoord[Z]++) {
    c[Z] = ((double)icoord[Z]) * spacing;
    grd_start = times( &tms_grd_start );

    for (icoord[Y] = -ne[Y]; icoord[Y] <= ne[Y]; icoord[Y]++) {
        c[Y] = ((double)icoord[Y]) * spacing;

        for (icoord[X] = -ne[X]; icoord[X] <= ne[X]; icoord[X]++) {
            c[X] = ((double)icoord[X]) * spacing;

/*
**  c[] contains the current grid point.
*/
            for (j = 0;  j < atom_maps_1;  j++) {
                enrg[j] = constant[j];
            }
            if (disGrid) {
                r_min = BIG;
            }

            donecovalent = FALSE;
/*
**  Do all Macromolecule (protein, DNA, etc.) atoms...
*/
            for (ia = 0;  ia < natom;  ia++) {
/*
**  Get distance, r, from current grid point, c, to this macromolecule atom, coord, 
*/
                for (i = 0;  i < XYZ;  i++) { 
                    d[i]  = coord[ia][i] - c[i]; 
                }
                r = hypotenuse( d[X], d[Y], d[Z] );
                if (r < APPROX_ZERO) {
                    r = APPROX_ZERO;
                }
                inv_r  = 1./r;
                for (i = 0;  i < XYZ;  i++) {
                    d[i] *= inv_r;
                }
                indx_r = min( lookup(r), MD_1 );

                if (disGrid) {
/* Calculate the so-called "Floating Grid"... */
                    r_min = min( r, r_min );
                }

                if (dddiel) {
/* Distance-dependent dielectric... */
                    enrg[elecPE] += charge[ia] * inv_r * epsilon[indx_r];
/* elecPE is the last grid map, i.e. electrostatics */
                } else {
/* Constant dielectric... */
                    enrg[elecPE] += charge[ia] * inv_r * invdielcal;
                }

/*
** If distance from grid point to atom ia is too large, 
** or if atom is a disordered hydrogen, 
**   add nothing to the grid-point's non-bond energy;
**   just continue to next atom...
*/
                if ( r > NBCUTOFF ) {
                    continue; /* onto the next atom... */
                }
                if ((atomtype[ia] == HYDROGEN) && (disorder[ia] == TRUE)) {
                    continue; /* onto the next atom... */
                }

                /*** racc = rdon = 1.; ***/
                racc = 1.;
                rdon = 1.;
                
                if (atomtype[ia] == HYDROGEN) {
/*
**  ia-th macromolecule atom = Hydrogen ( 4 = H )
**  => Protein H-bond donor, OH or NH.
**  calculate racc for H-bond ACCEPTOR PROBES at this grid pt.
**            ====     ======================
*/
                    cos_theta = 0.;
/*
**  d[] = Unit vector from current grid pt to ia_th m/m atom.
**  cos_theta = d dot rvector == cos(angle) subtended.
*/
                    for (i = 0;  i < XYZ;  i++) {
                        cos_theta -= d[i] * rvector[ia][i];
                    }

                    if (cos_theta <= 0.) {
/*
**  H->current-grid-pt vector >= 90 degrees from
**  N->H or O->H vector, 
*/
                        racc = 0.;
                    } else {
/*
**  racc = [cos(theta)]^2.0 for N-H
**  racc = [cos(theta)]^4.0 for O-H, 
*/
                        switch( rexp[ia] ) {
                            case 1:
                            default:
                                racc = cos_theta;
                                break;
                            case 2:
                                racc = cos_theta*cos_theta;
                                break;
                            case 4:
                                tmp = cos_theta*cos_theta;
                                racc = tmp*tmp;
                                break;
                        }
                        /* racc = pow( cos_theta, (double)rexp[ia] ); */
                    }
                    /* endif (atomtype[ia] == HYDROGEN) */
                } else if ((atomtype[ia] == OXYGEN) && (disorder[ia] == FALSE)) {
                    /*
                    **  ia-th macromolecule atom = Oxygen 
                    **  => Protein H-bond acceptor, oxygen.
                    */

                    rdon = 0.;

                    /* check to see that probe is in front of oxygen, not behind */
                    cos_theta = 0.;
                    for (i = 0;  i < XYZ;  i++) {
                        cos_theta -= d[i] * rvector[ia][i];
                    }
                    /*
                    ** t0 is the angle out of the lone pair plane, calculated 
                    ** as 90 deg - acos (vector to grid point DOT lone pair 
                    ** plane normal)
                    */
                    t0 = 0.;
                    for (i = 0;  i < XYZ;  i++) {
                        t0 += d[i] * rvector2[ia][i];
                    }
                    if (t0 > 1.) {
                        t0 = 1.;
                        (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %.3f, a value greater than 1.\n", t0);
                    } else if (t0 < -1.) {
                        t0 = -1.;
                        (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %.3f, a value less than 1.\n", t0);
                    }
                    t0 = PI_halved - acos(t0);

                    /*
                    ** ti is the angle in the lone pair plane, away from the 
                    ** vector between the lone pairs, 
                    ** calculated as (grid vector CROSS lone pair plane normal)
                    ** DOT C=O vector - 90 deg
                    */
                    cross[0] = d[1] * rvector2[ia][2] - d[2] * rvector2[ia][1];
                    cross[1] = d[2] * rvector2[ia][0] - d[0] * rvector2[ia][2];
                    cross[2] = d[0] * rvector2[ia][1] - d[1] * rvector2[ia][0];
                    rd2 = sq(cross[0])+sq(cross[1])+sq(cross[2]);
                    if (rd2 < APPROX_ZERO) {
                        if ((rd2 == 0.) && (warned == 'F')) {
                            (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                            (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia+1, ib+1);
                            warned = 'T';
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1./sqrt(rd2);
                    ti = 0.;
                    for (i = 0;  i < XYZ;  i++) {
                        ti += cross[i] * inv_rd * rvector[ia][i];
                    }

                    /* rdon expressions from Goodford */
                    rdon = 0.;
                    if (cos_theta >= 0.) {
                        if (ti > 1.) {
                            ti = 1.;
                            (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %.3f, a value greater than 1.\n", ti);
                        } else if (ti < -1.) {
                            ti = -1.;
                            (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %.3f, a value less than 1.\n", ti);
                        }
                        ti = acos(ti) - PI_halved;
                        if (ti < 0.) {
                            ti = -ti;
                        }
                        /* the 2.0*ti can be replaced by (ti+ti) in: rdon = (0.9+0.1*sin(2.0*ti))*cos(t0);*/
                        rdon = (0.9+0.1*sin(ti+ti))*cos(t0);
                    /* 0.34202 = cos (100 deg) */
                    } else if (cos_theta >= -0.34202) {
                        rdon = 562.25*pow(0.116978-sq(cos_theta), 3.)*cos(t0); 
                    }

                    /* endif atomtype == OXYGEN, not disordered */
                } else if ((atomtype[ia] == OXYGEN) && (disorder[ia] == TRUE)) {

                    /* cylindrically disordered hydroxyl */
                    cos_theta = 0.;
                    for (i = 0;  i < XYZ;  i++) {
                        cos_theta -= d[i] * rvector[ia][i];
                    }        
                    if (cos_theta > 1.) {
                        cos_theta = 1.;
                        (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %.3f, a value greater than 1.\n", cos_theta);
                    } else if (cos_theta < -1.) {
                        cos_theta = -1.;
                        (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %.3f, a value less than 1.\n", cos_theta);
                    }
                    theta = acos(cos_theta);
                    racc = 0.;
                    rdon = 0.;
                    if (theta <= 1.24791 + PI_halved) {
                        /* 1.24791 rad = 180 deg minus C-O-H bond angle, 
                        ** 108.5 deg */
                        rdon = pow(cos(theta-1.24791), 4.);
                        racc = rdon;
                    }
                } /* atomtype test */

/*
** For each probe atom-type, 
** Sum pairwise interactions between each probe
** and the current receptor atom...
*/
                for (probe = 0;  probe < atom_maps;  probe++) {

                    if (iscovalent[probe] == TRUE) {
                      
                        if (donecovalent == FALSE) {
                            /* calculate the distance from the current
                             * grid point to the covalent attachment point */
                            for (ii = 0;  ii < XYZ;  ii++) { 
                                d[ii]  = covpos[ii] - c[ii]; 
                            }
                            rcov = hypotenuse( d[X], d[Y], d[Z] );
                            rcov = rcov / covhalfwidth;
                            if (rcov < APPROX_ZERO) {
                                rcov = APPROX_ZERO;
                            }
                            enrg[probe] = covbarrier * (1. - exp(ln_half * rcov * rcov)) + constant[probe];
                            /* prevent unnecessary repetitions of this
                             * calculation as we loop over all the
                             * macromolecule's atoms... */
                            donecovalent = TRUE;
                        }

                    } else {
                        if (hbond[probe] == TRUE) {

                            /* rsph ramps in angular dependence for distances with negative energy */
                            rsph = energy[atomtype[ia]][indx_r][probe]/100.;
                            rsph = max(rsph, 0.);
                            rsph = min(rsph, 1.);
                            
                            if ((atmtypstr[probe] == 'N') ||
                                (atmtypstr[probe] == 'O') ||
                                (atmtypstr[probe] == 'S')) {

                                  /* PROBE can be an H-BOND ACCEPTOR, */
                                if (disorder[ia] == FALSE ) {
                                    enrg[probe] += energy[atomtype[ia]][indx_r][probe] * (racc+(1.-racc)*rsph);
                                } else {
                                    enrg[probe] += energy[HYDROGEN][max(0, indx_r-110)][probe] * (racc+(1.-racc)*rsph);
                                }

                            } else if (atmtypstr[probe] == 'H') {

                                /*  PROBE is H-BOND DONOR, */
                                enrg[probe] += energy[atomtype[ia]][indx_r][probe] *
                                                             (rdon+(1.-rdon)*rsph);
                            }

                        } else {

                                /*  PROBE does not form H-bonds..., */
                                enrg[probe] += energy[atomtype[ia]][indx_r][probe];

                        }/* hbond test */

                        /* add desolvation energy (notice sign) */
                        if ((atomtype[ia] != HYDROGEN) && (atmtypstr[probe] != 'H')) {
                            /* Full ligand and protein Stouten solvation energy */
                            /* enrg[probe] -= (solpar_probe[probe]*vol[ia] + 
                            ** solpar[ia]*vol_probe[probe]) * sol_fn[indx_r]; */

                            /* Ligand only Stouten energy */
                            enrg[probe] -= solpar_probe[probe] * vol[ia] * sol_fn[indx_r];
                        }
                    } /* is not covalent */
                }/* probe */
            }/* ia loop, over all macromolecule atoms... */

/*
** O U T P U T . . .
**
** Now output this grid point's energies to the maps:
**
*/
            for (k = 0;  k < atom_maps_1;  k++) {
                if (!problem_wrt) {
                    if (fabs(enrg[k]) < PRECISION) {
                        fprintf_retval = fprintf(mapFile[k], "0.\n");
                    } else {
                        fprintf_retval = fprintf(mapFile[k], "%.3f\n", (float)enrg[k]);
                    }
                    if (fprintf_retval < 0) {
                        problem_wrt = TRUE;
                    }
                }

                emax[k] = max( emax[k], enrg[k] );
                emin[k] = min( emin[k], enrg[k] );
            }
            if (disGrid) {
                if ((!problem_wrt)&&(fprintf(disFile, "%.3f\n", (float)r_min) < 0)) {
                    problem_wrt = TRUE;
                }
            }

        } /* icoord[X] loop */

    } /* icoord[Y] loop */

    if (problem_wrt) {
        (void) fprintf( stderr, "WARNING! Problems writing grid maps - there may not be enough space.\n");
        (void) fprintf( logFile, "WARNING! Problems writing grid maps - there may not be enough space.\n");
    }
    grd_end = times( &tms_grd_end );
    ++nDone;
    timeRemaining = (float)(grd_end - grd_start) * idct * (float)(n1[Z] - nDone);
    (void) fprintf( logFile, " %6d   %8.3lf   %5.1lf%%   ", icoord[Z], cgridmin[Z]+c[Z], percentdone*(double)++ic);
    prHMSfixed( timeRemaining );
    (void) fprintf( logFile, "  ");
    timesys( grd_end - grd_start, &tms_grd_start, &tms_grd_end );
    (void) fflush( logFile );

} /* icoord[Z] loop */

/*____________________________________________________________________________
** Print a summary of extrema-values from the atomic-affinity and
** electrostatics grid-maps, 
**____________________________________________________________________________*/

(void) fprintf(logFile, "\nGrid\tAtom\tMinimum\t\tMaximum\n");
(void) fprintf(logFile, "Map \tType\tEnergy \t\tEnergy \n");
(void) fprintf(logFile, "\t\t(kcal/mol)\t\t(kcal/mol)\n");
(void) fprintf(logFile, "___\t____\t_____________\t_____________\n");

for (i = 0;  i < atom_maps;  i++) {
    (void) fprintf( logFile, " %d\t %c\t  %6.2lf\t%6.2le\n", i+1, atmtypstr[i], emin[i], emax[i]);
}

(void) fprintf( logFile, " %d\t %c\t  %6.2lf\t%6.2le\tElectrostatic Potential\n", atom_maps+1, 'e', emin[elecPE], emax[i]);

(void) fprintf( logFile, "\n\n * Note:  Every pairwise-atomic interaction was clamped at %.2f\n\n", EINTCLAMP );

/*
** Close all files, ************************************************************
*/

for (i = 0;  i < atom_maps_1;  i++) {
    (void) fclose( mapFile[i] );
}
if (disGrid) {
    (void) fclose(disFile);
}
(void) fprintf( stderr, "\n%s: Successful Completion.\n", programname);
(void) fprintf( logFile, "\n%s: Successful Completion.\n", programname);

job_end = times( &tms_job_end );
timesyshms( job_end - job_start, &tms_job_start, &tms_job_end );

(void) fclose( logFile);

return 0;
}
/*
** End of main function.
*/
