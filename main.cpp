/* main.c */
/*
  $Id: main.cpp,v 1.32 2005/05/24 00:08:46 gillet Exp $
*/


#include <sys/types.h>
#ifndef _WIN32
#include <sys/times.h>
#include <unistd.h> /* long sysconf(int name) */
#include <sys/param.h>
#else
#include "times.h"
#include <Winsock2.h>
#include "util.h"
#endif

#include <math.h>
#include <stdio.h>
#include <search.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h> /* tolower */

#include <stddef.h> 
#include <ctype.h> 


/* the BOINC API header file */
#ifdef BOINC

#include "diagnostics.h"
#include "boinc_api.h" 
#include "filesys.h" 		// boinc_fopen(), etc... */

#endif

#include "autogrid.h"
#include "autoglobal.h"
#include "autocomm.h"

extern float idct;


/* fopen rewrite to either use BOINC api or normal system call */
FILE *ag_fopen(const char *path, const char *mode)
{
  FILE *filep;
#ifdef BOINC
  int rc;
  char resolved_name[512];
  rc = boinc_resolve_filename(path, resolved_name, sizeof(resolved_name));
  if (rc){
      fprintf(stderr, "BOINC_ERROR: cannot open filename.%s\n",path);
      boinc_finish(rc);    /* back to BOINC core */
    }
    // Then open the file with boinc_fopen() not just fopen()
    filep = boinc_fopen(resolved_name, mode);
#else
    filep = fopen(path,mode);
#endif
    return filep;
}

static int get_rec_index(const char key[]);


int main( int argc,  char **argv )

/******************************************************************************/
/*      Name: main (executable's name is "autogrid").                         */
/*  Function: Calculation of interaction energy grids for Autodock.           */
/*            Directional H_bonds from Goodford:                              */
/*            Distance dependent dielectric after Mehler and Solmajer.        */
/*            Charge-based desolvation                                        */
/* Copyright: (C) 2004, TSRI                                                  */
/*                                                                            */
/*   Authors: Garrett Matthew Morris, Ruth Huey, David S. Goodsell            */
/*                                                                            */
/*            The Scripps Research Institute                                  */
/*            Department of Molecular Biology, MB5                            */
/*            10550 North Torrey Pines Road                                   */
/*            La Jolla, CA 92037-1000.                                        */
/*                                                                            */
/*            e-mail: garrett@scripps.edu                                     */
/*                    rhuey@scripps.edu                                       */
/*                    goodsell@scripps.edu                                    */
/*                                                                            */
/*            Helpful suggestions and advice:                                 */
/*            Arthur J. Olson                                                 */
/*            Bruce Duncan, Yng Chen, Michael Pique, Victoria Roberts         */
/*            Lindy Lindstrom                                                 */
/*                                                                            */
/*      Date: 07/07/04                                                        */
/*                                                                            */
/*    Inputs: Control file, receptor PDBQT file, parameter file               */
/*   Returns: Atomic affinity, desolvation and electrostatic grid maps.       */
/*   Globals: MAX_DIST, MAX_MAPS                                              */
/*             increased from 8 to 16  6/4/2004                               */
/*                                                                            */
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 07/06/89 DSG     FORTRAN implementation                                    */
/* 07/05/92 GMM     C translation                                             */
/* 20/09/95 GMM/DSG AutoGrid3                                                 */
/* 07/07/04 DSG/RH  AutoGrid4                                                 */
/******************************************************************************/

/* Note: 21/03/03 GMM note: ATOM_MAPS is no longer used here; was used for
 * is_covalent and is_hbonder, but these are now folded into the MapObject 
 * and arrayed up to MAX_MAPS (currently). MAX_MAPS is always larger than 
 * ATOM_MAPS, so this is safe. */

{

/*  for associative dictionary storing parameters by autogrid 'type'  */
FILE * dataFile;
char dataline[100];
//ENTRY item; 
/*see  atom_parameter_manager.c */
static ParameterEntry thisparm;
ParameterEntry * newparm;
ParameterEntry * found_parm;
char param_filename[MAX_CHARS];



/* LIGAND: 
 *  maximum is MAX_MAPS */
/*each type is now at most two characters plus '\0'*/
/* currently ligand_atom_types is sparse... 
 * some types are not set*/
/*FIX THIS: initialize to "" */
char ligand_types[MAX_MAPS][3];

/*array of ptrs used to parse input line*/
char * ligand_atom_types[MAX_MAPS];

#define NUM_RECEPTOR_TYPES  NUM_ALL_TYPES
/*malloc this after the number of receptor types is parsed*/
/*MAX_DIST is really NBCUTOFF times 100 */
static double energy_lookup[NUM_RECEPTOR_TYPES][MAX_DIST][MAX_MAPS];


typedef struct mapObject {
    int    atom_type; /*corresponds to receptor numbers????*/
    int    map_index;
    int    is_covalent;
    int    is_hbonder;
    FILE   *map_fileptr;
    char   map_filename[MAX_CHARS];
    char   type[3]; /*eg HD or OA or NA or N*/
    double constant; /*this will become obsolete*/
    double energy_max;
    double energy_min;
    double energy;
    double vol_probe;
    double solpar_probe;
    /*new 6/28*/
    double Rij;
    double epsij;
    hbond_type hbond; /*hbonding character: */
    double Rij_hb;
    double epsij_hb;
    /*per receptor type parameters, ordered as in receptor_types*/
    double nbp_r[NUM_RECEPTOR_TYPES]; /*radius of energy-well minimum*/
    double nbp_eps[NUM_RECEPTOR_TYPES];/*depth of energy-well minimum*/
    int xA[NUM_RECEPTOR_TYPES]; /*generally 12*/
    int xB[NUM_RECEPTOR_TYPES]; /*6 for non-hbonders 10 for h-bonders*/
    int hbonder[NUM_RECEPTOR_TYPES];
} MapObject;
/*constant will go away*/

char * maptypeptr; /*ptr for current map->type*/
MapObject *gridmap; /* was statically assigned  MapObject gridmap[MAX_MAPS]; */

/* needed to make regression tests work between platforms*/
FloatOrDouble *dummy_map;

/*variables for RECEPTOR:*/
/*each type is now at most two characters, eg 'NA\0'*/
/*NB: these are sparse arrays, some entries are not set */
/*FIX THIS: initialize to "" */
char receptor_types[NUM_RECEPTOR_TYPES][3];
/* number of different receptor atom types*/
int receptor_types_ct = 0;
/* array of numbers of each type */
/*NB: this is a sparse int array, some entries are 0*/
/*FIX THIS: initialize to 0*/
int receptor_atom_type_count[NUM_RECEPTOR_TYPES];

/*array of ptrs used to parse input line*/
char * receptor_atom_types[NUM_RECEPTOR_TYPES];


/* MAX_ATOMS */
double charge[MAX_ATOMS];
double vol[MAX_ATOMS];
double solpar[MAX_ATOMS];
/*integers are simpler!*/
int atom_type[MAX_ATOMS];
hbond_type hbond[MAX_ATOMS];
int disorder[MAX_ATOMS];
int rexp[MAX_ATOMS];
double coord[MAX_ATOMS][XYZ];
double rvector[MAX_ATOMS][XYZ];
double rvector2[MAX_ATOMS][XYZ];

/*canned atom type number*/
int hydrogen, carbon, arom_carbon, oxygen, nitrogen; 
int nonHB_hydrogen, nonHB_nitrogen, sulphur, nonHB_sulphur;

/* XYZ */
double cross[XYZ];
double c[XYZ];
double cext[XYZ];
double cgridmax[XYZ];
double cgridmin[XYZ];
double cmax[XYZ];
double cmin[XYZ];
double csum[XYZ];
double cmean[XYZ];
double center[XYZ];
double covpos[XYZ]; /* Cartesian-coordinate of covalent affinity well. */
double d[XYZ];
double dc[XYZ];
int icoord[XYZ]; /* int icenter; */
int ne[XYZ];
int n1[XYZ];
int nelements[XYZ];

/* MAX_CHARS */
char AVS_fld_filename[MAX_CHARS];
char floating_grid_filename[MAX_CHARS];
char host_name[MAX_CHARS];
char receptor_filename[MAX_CHARS];
char xyz_filename[MAX_CHARS];

/* LINE_LEN */
char message[LINE_LEN];
char line[LINE_LEN];
char GPF_line[LINE_LEN];
int length = LINE_LEN;

/* MAX_DIST */
double epsilon[MAX_DIST];
int MD_1 = MAX_DIST - 1;
static double sol_fn[MAX_DIST];
double energy_smooth[MAX_DIST];

/*JUNE 29*/
int ctr;
/*JUNE 29*/

char atom_name[6];
/*char extension[5];*/
/* char q_str[7]; */
char record[6];
char temp_char = ' ';
char token[5];
char warned = 'F';
char xyz[5];

FILE *receptor_fileptr,
     *AVS_fld_fileptr,
     *xyz_fileptr,
     *floating_grid_fileptr;

/*initialize ff values*/
double FE_coeff_vdW    = 0.14880;
double FE_coeff_hbond  = 0.12667;
double FE_coeff_estat  = 0.14084;
double FE_coeff_desolv = 0.122;
double FE_coeff_tors   = 0.090277;
/*for NEW3 desolvation terms*/
double solpar_q = .01097;  /*unweighted value restored 3:9:05 */
/*double solpar_q = 0.0013383; =.01097 * 0.122*/


/*for NEW3 desolvation terms*/
/*double solpar_q = 0.0013383; =.01097 * 0.122*/
/*double solpar_q = 0.001118; */

double A, epsilon0, rk, lambda, B, lambda_B;
double q_tot = 0.0;
double diel, invdielcal;
double dxA;
double dxB;
double minus_inv_two_sigma_sqd;
double percentdone;
double PI_halved;
double q_max = -BIG,  q_min = BIG;
double rA;
double rB; /* double e; */
double rcov = 0.0; /* Distance from current grid point to the covalent attachment point */
double ri, inv_rd, rd2, r; /* re, r2, rd, */
double r_min, inv_r, inv_rmax, racc, rdon, rsph, cos_theta, theta, tmp;
double r_smooth = 0.;
double rdot;
double Rij, epsij;
double spacing = 0.375; /* One quarter of a C-C bond length. */
double t0, ti;
double ln_half = 0.0;
double covhalfwidth = 1.0;
double covbarrier = 1000.0;
double cA, cB, tmpconst;
double sigma;
double version_num = 4.00;
/*are these necessary??*/
double temp_vol, temp_solpar;
double temp_hbond_enrg, hbondmin[MAX_MAPS], hbondmax[MAX_MAPS];
double rmin, Hramp;

/*int num_rec_types = 0;*/

float timeRemaining = 0.;

int num_maps = 0;
int num_atom_maps = -1;
int floating_grid = FALSE, dddiel = FALSE, disorder_h = FALSE;
int elecPE = 0;
int dsolvPE = 0;
/* int covmap; */
int from, to;
int fprintf_retval = 0;
int GPF_keyword = -1;
int indcom = 0;
int infld;
int nbond;
int nDone = 0;
int problem_wrt = FALSE;
int xA, xB;
int hbondflag[MAX_MAPS];


#define INIT_NUM_GRID_PTS -1
int num_grid_points_per_map = INIT_NUM_GRID_PTS;
register int i = 0, ii = 0, j = 0, k = 0, indx_r = 0, i_smooth = 0;
/* register int i = 0, ii = 0, j = 0, jj = 0, k = 0, indx_r = 0, i_smooth = 0;
 */
register int ia = 0, ib = 0, ic = 0, map_index = -1, iat = 0, i1 = 0, i2 = 0, i3 = 0;
register int closestH = 0;

static int num_receptor_atoms;
static long clktck = 0;

Clock      job_start;
Clock      job_end;
struct tms tms_job_start;
struct tms tms_job_end;

Clock      grd_start;
Clock      grd_end;
struct tms tms_grd_start;
struct tms tms_grd_end;

#ifdef BOINC
    int flags = 0;
    int rc;
    flags =
      BOINC_DIAG_DUMPCALLSTACKENABLED |
      BOINC_DIAG_HEAPCHECKENABLED |
      BOINC_DIAG_REDIRECTSTDERR |
      BOINC_DIAG_REDIRECTSTDOUT ;
    boinc_init_diagnostics(flags);

#ifdef BOINCCOMPOUND
    BOINC_OPTIONS options;
    options.main_program = false;
    options.check_heartbeat = false;// monitor does check heartbeat
    options.handle_trickle_ups = false;
    options.handle_trickle_downs = false;
    options.handle_process_control = false;
    options.send_status_msgs = true;// only the worker programs (i.e. model) sends status msgs
    options.direct_process_action = true;// monitor handles suspend/quit, but app/model doesn't
    // Initialisation of Boinc 
    rc =  boinc_init_options(options); //return 0 for success
    if( rc ){
      fprintf(stderr,"BOINC_ERROR: boinc_init_options() failed \n");
      exit(rc);
    }

#else
    // All BOINC applications must initialize the BOINC interface:
    rc = boinc_init();
    if (rc){
      fprintf(stderr, "BOINC_ERROR: boinc_init() failed.\n");
      exit(rc);
    }
#endif

#endif

/*
 * Fetch clock ticks per second.
 */
#ifdef _WIN32
if (clktck == 0) {
    if ( (clktck = CLOCKS_PER_SEC) < 0) {
        (void) fprintf( stderr, "\"CLOCKS_PER_SEC\" command failed in \"main.c\"\n");
        (void) fprintf( logFile, "\"CLOCKS_PER_SEC\" command failed in \"main.c\"\n");
        exit(-1);
    } else {
        idct = (float)1. / (float)clktck;
    }
}

#else
if (clktck == 0) {
    if ( (clktck = sysconf(_SC_CLK_TCK)) < 0) {
        (void) fprintf( stderr, "\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n");
        (void) fprintf( logFile, "\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n");
        exit(-1);
    } else {
        idct = 1. / (float)clktck;
    }
}
#endif

ln_half = (double) log(0.5);

/*
 * Get the time at the start of the run...
 */
job_start = times( &tms_job_start);

/*
 * Parse the command line...
 */
(void) setflags( argc, argv);

(void) strcpy( xyz, "xyz");

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
 * Initialize int receptor_atom_type_count[] array to 0
 */
for (i=0; i<NUM_RECEPTOR_TYPES; i++) {
    receptor_atom_type_count[i] = 0;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * Output the "AutoGrid" banner...
 */
banner( version_num);

/*
 * Print out MAX_MAPS - maximum number of maps allowed
 */
(void) fprintf(logFile, "Maximum number of maps that can be computed = %d (defined by MAX_MAPS in \"autocomm.h\").\n\n", MAX_MAPS);


/*
 * Print the time and date when the file was created...
 */
(void) fprintf( logFile, "This file was created at:\t\t\t");
printdate( logFile, 1);

#ifndef _WIN32
if (gethostname( host_name, MAX_CHARS ) == 0) {
    (void) fprintf( logFile, "                   using:\t\t\t\"%s\"\n", host_name);
}
#endif
/******************************************************************************/

/* Read in the grid parameter file...  */

while( fgets( GPF_line, LINE_LEN, GPF_fileptr) != NULL ) {

/******************************************************************************/

    GPF_keyword = gpfparser( GPF_line);

    /* This first "switch" figures out how to echo the current GPF line. */

    switch( GPF_keyword ) {

    case -1:
        (void) fprintf( logFile, "GPF> %s", GPF_line);
        (void) sprintf( message, "%s: WARNING: Unrecognized keyword in grid parameter file.\n", programname);
        (void) fprintf( logFile, "%s", message);
        (void) fprintf( stderr,  "%s", message);
        continue;  /* while fgets GPF_line... */

    case GPF_NULL:
    case GPF_COMMENT:
        (void) fprintf( logFile, "GPF> %s", GPF_line);
        (void) fflush( logFile);
        break;

    default:
        (void) fprintf( logFile, "GPF> %s", GPF_line);
        indcom = strindex( GPF_line, "#");
        if (indcom != -1) {
            GPF_line[ indcom ] = '\0'; /* Truncate str. at the comment */
        }
        (void) fflush( logFile);
        break;

    } /* first switch */

/******************************************************************************/

    /* This second switch interprets the current GPF line. */

    switch( GPF_keyword ) {

/******************************************************************************/

    case GPF_NULL:
    case GPF_COMMENT:
        break;

/******************************************************************************/

    case GPF_RECEPTOR:
        /* read in the receptor filename */

        (void) sscanf( GPF_line, "%*s %s", receptor_filename);
        (void) fprintf( logFile, "Receptor Input File :\t%s\n\nAtom Type Assignments:\n\n", receptor_filename);

        /* try to open receptor file */
        if ( (receptor_fileptr = ag_fopen(receptor_filename, "r")) == NULL ) {
            (void) fprintf( stderr, "\n%s: can't find or open receptor-file %s\n", programname, receptor_filename);
            (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
            (void) fprintf( logFile, "\n%s: can't find or open receptor-file %s\n", programname, receptor_filename);
            (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
            exit(-1);
        }

        /* start to read in the lines of the receptor file */
        ia = 0;
        while ( (fgets(line, length, receptor_fileptr)) != NULL ) {
            (void) sscanf(line, "%6s", record);
            if (equal(record, "ATOM", 4) || /* Amino Acid or DNA/RNA atoms */
                equal(record, "HETA", 4) || /* Non-standard heteroatoms */
                equal(record, "CHAR", 4)) { /* Partial Atomic Charge - not a PDB record */

                (void) strncpy( atom_name, &line[12], 4);
                /* atom_name is declared as an array of 6 characters,
                 * the PDB atom name is 4 characters (C indices 0, 1, 2 and 3)
                 * but let's ensure that the fifth character (C index 4)
                 * is a null character, which terminates the string. */
                atom_name[4] = '\0';

                /* Output the serial number of this atom... */
                (void) fprintf( logFile, "Atom no. %d,  \"%s\" ", ia + 1, atom_name);
                (void) fflush( logFile);

                /* Read in this receptor atom's coordinates,partial charges, and
                 * solvation parameters in PDBQS format... */

                (void) sscanf(&line[30], "%lf", &coord[ia][X]);
                (void) sscanf(&line[38], "%lf", &coord[ia][Y]);
                (void) sscanf(&line[46], "%lf", &coord[ia][Z]);

                /* Output the coordinates of this atom... */
                (void) fprintf( logFile, " at (%.3lf, %.3lf, %.3lf), ",
                                coord[ia][X], coord[ia][Y], coord[ia][Z]);
                (void) fflush( logFile);

                /*1:CHANGE HERE: need to set up vol and solpar*/
                (void) sscanf(&line[70], "%lf", &charge[ia]);
                //printf("new type is: %s\n", &line[77]);
                (void) sscanf(&line[77], "%s", thisparm.autogrid_type);
                found_parm = apm_find(thisparm.autogrid_type);
                if (found_parm != NULL){
                    if (found_parm->rec_index<0){
                        found_parm->rec_index = receptor_types_ct++;
                    }
                    atom_type[ia] = found_parm->rec_index;
                    solpar[ia] = found_parm->solpar;
                    vol[ia] = found_parm->vol;
                    hbond[ia] = found_parm->hbond; /*NON=0, DS,D1, AS, A1, A2*/
#ifdef DEBUG
                    printf("%d:key=%s, type=%d,solpar=%f,vol=%f\n",ia,thisparm.autogrid_type, atom_type[ia],solpar[ia],vol[ia]);
#endif
                    ++receptor_atom_type_count[found_parm->rec_index];
                } else {
                    fprintf(logFile, "\n\nreceptor file contains unknown type: '%s'\nadd parameters for it to the parameter library first\n", thisparm.autogrid_type);
                    exit(-1);
                }
                    
                /* if from pdbqs: convert cal/molA**3 to kcal/molA**3 */
                /*solpar[ia] *= 0.001;*/

                q_max = max(q_max, charge[ia]);
                q_min = min(q_min, charge[ia]);

                if (atom_name[0] == ' ') {
                    /* truncate the first character... */
                    atom_name[0] = atom_name[1];
                    atom_name[1] = atom_name[2];
                    atom_name[2] = atom_name[3];
                    atom_name[3] = '\0';
                } else if (atom_name[0] == '0' ||
                    atom_name[0] == '1' ||
                    atom_name[0] == '2' ||
                    atom_name[0] == '3' ||
                    atom_name[0] == '4' ||
                    atom_name[0] == '5' ||
                    atom_name[0] == '6' ||
                    atom_name[0] == '7' ||
                    atom_name[0] == '8' ||
                    atom_name[0] == '9') {
                    if (atom_name[1] == 'H') {
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
                        temp_char    = atom_name[0];
                        atom_name[0] = atom_name[1];
                        atom_name[1] = atom_name[2];
                        atom_name[2] = atom_name[3];
                        atom_name[3] = temp_char;
                    }
                }

                /*FIX THIS***
                 * ?????????*/

                /*if (atom_name[0] == 'L') {
                    (void) fprintf( logFile, "\tWARNING: lone pair atom type found; will be treated as atom type %d, \"%c\"\n", METAL + 1, receptor_atom_type_string[METAL]);*/
                    /* Treat Lone Pair centres as if they were a Metal...  */
                /*    atom_type[ia] = METAL;
                } else {*/
                    /*
                     * The default atom type is 'X',
                     * atom_type[ia] = UNKNOWN;
                     */
                  /*  atom_type[ia] = get_atom_type(atom_name, receptor_atom_type_string);
                }*/
                /* Tell the user what you thought this atom was... */
                (void) fprintf( logFile, "  was atom type \"%s\" assigned %d. atom_type->%d\n", found_parm->autogrid_type, found_parm->rec_index, atom_type[ia]);
                /*(void) fprintf( logFile, "  was assigned atom type %d, \"%c\".\n", atom_type[ia] + 1, receptor_atom_type_string[atom_type[ia]]);*/
                (void) fflush( logFile);

                /* Count the number of each atom type */
                /*++receptor_atom_type_count[ atom_type[ia] ];*/

                /* Keep track of the extents of the receptor */
                for (i = 0;  i < XYZ;  i++) {
                    cmax[i] = max(cmax[i], coord[ia][i]);
                    cmin[i] = min(cmin[i], coord[ia][i]);
                    csum[i] += coord[ia][i];
                }
                /* Total up the partial charges as we go... */
                q_tot += charge[ia];

                /* Increment the atom counter */
                ia++;

                /* Check that there aren't too many atoms... */
                if (ia > MAX_ATOMS) {
                    (void) fprintf( logFile, "Error      : Sorry, AutoGrid cannot continue.\n");
                    (void) fprintf( logFile, "           : Too many atoms in receptor input file %s;\n", receptor_filename);
                    (void) fprintf( logFile, "           : -- the maximum number of atoms, MAX_ATOMS, allowed is %d.\n", MAX_ATOMS);
                    (void) fprintf( logFile, "Suggestion : Increase the value in the \"#define MAX_ATOMS %d\" line", MAX_ATOMS);
                    (void) fprintf( logFile, "           : in the source file \"autogrid.h\", and re-compile AutoGrid.\n");
                    (void) fflush( logFile);
                    exit(-1);
                } /* endif */
            } /* endif */
        } /* endwhile */
        (void) fclose( receptor_fileptr);
        /* Update the total number of atoms in the receptor */
        num_receptor_atoms = ia;
        (void) fprintf( logFile, "\nMaximum partial atomic charge found = %+.3lf e\n", q_max);
        (void) fprintf( logFile, "Minimum partial atomic charge found = %+.3lf e\n\n", q_min);
        (void) fflush( logFile);
        /* Check there are partial charges... */
        if (q_max == 0. && q_min == 0.) {
            (void) fprintf( logFile, "\nERROR!  Partial atomic charges not found!\n\n");
            exit(-1);
        } /* if  there are no charges EXIT*/

        for (ia = 0;  ia < num_receptor_atoms;  ia++) {
            rexp[ia] = 0;
        }
        (void) fprintf( logFile, "Atom\tAtom\tNumber of this Type\n");
        (void) fprintf( logFile, "Type\t ID \t in Receptor\n");
        (void) fprintf( logFile, "____\t____\t___________________\n");
        (void) fflush( logFile);
        /*2. CHANGE HERE: need to count number of each receptor_type*/
        for (ia = 0;  ia < NUM_RECEPTOR_TYPES;  ia++) {
            i = 0;
            if(receptor_atom_type_count[ia]!=0){
                (void) fprintf( logFile, " %d\t %s\t\t%6d\n", (ia), receptor_types[ia], receptor_atom_type_count[ia]);
                i++;
            };
        }
        (void) fprintf( logFile, "\nTotal number of atoms :\t\t%d atoms \n", num_receptor_atoms);
        (void) fflush( logFile);
        (void) fprintf( logFile, "Total charge :\t\t\t%.2lf e\n", q_tot);
        (void) fflush( logFile);
        (void) fprintf( logFile, "\n\nReceptor coordinates fit within the following volume:\n\n");
        (void) fflush( logFile);
        (void) fprintf( logFile, "                   _______(%.1lf, %.1lf, %.1lf)\n", cmax[X], cmax[Y], cmax[Z]);
        (void) fprintf( logFile, "                  /|     /|\n");
        (void) fprintf( logFile, "                 / |    / |\n");
        (void) fprintf( logFile, "                /______/  |\n");
        (void) fprintf( logFile, "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n", (cmax[X] + cmin[X])/2., (cmax[Y] + cmin[Y])/2., (cmax[Z] + cmin[Z])/2.);
        (void) fprintf( logFile, "                |  /   |  /\n");
        (void) fprintf( logFile, "                | /    | /\n");
        (void) fprintf( logFile, "                |/_____|/\n");
        (void) fprintf( logFile, "(%.1lf, %.1lf, %.1lf)      \n", cmin[X], cmin[Y], cmin[Z]);
        (void) fprintf( logFile, "\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cmax[X], cmax[Y], cmax[Z]);
        (void) fprintf( logFile, "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cmin[X], cmin[Y], cmin[Z]);
        (void) fprintf( logFile, "\n");
        cmean[0] = csum[0] / (double)num_receptor_atoms;
        cmean[1] = csum[1] / (double)num_receptor_atoms;
        cmean[2] = csum[2] / (double)num_receptor_atoms;
        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_GRIDFLD:
        (void) sscanf( GPF_line, "%*s %s", AVS_fld_filename);
        infld = strindex( AVS_fld_filename, ".fld");
        if (infld == -1) {
            (void) fprintf( stderr, "\n\n%s: WARNING:  Grid data file needs the extension \".fld\" for AVS input\n\n", programname);
            (void) fprintf( logFile, "\n\n%s: WARNING:  Grid data file needs the extension \".fld\" for AVS input\n\n", programname);
            exit(-1);
        } else {
            infld = strindex( AVS_fld_filename, "fld");
            (void) strcpy(xyz_filename, AVS_fld_filename);
            xyz_filename[infld] = 'x';
            xyz_filename[infld + 1] = 'y';
            xyz_filename[infld + 2] = 'z';
        }
        if ( (AVS_fld_fileptr = ag_fopen(AVS_fld_filename, "w")) == NULL ) {
            (void) fprintf( stderr, "\n%s: can't create grid dimensions data file %s\n", programname, AVS_fld_filename);
            (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
            (void) fprintf( logFile, "\n%s: can't create grid dimensions data file %s\n", programname, AVS_fld_filename);
            (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
            exit(-1);
        } else {
            (void) fprintf( logFile, "\nCreating (AVS-readable) grid maps file : %s\n", AVS_fld_filename);
        }
        if ( (xyz_fileptr = ag_fopen(xyz_filename, "w")) == NULL ) {
            (void) fprintf( stderr, "\n%s: can't create grid extrema data file %s\n", programname, xyz_filename);
            (void) fprintf( stderr, "%s: SORRY!    unable to create the \".xyz\" file.\n\n", programname);
            (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
            (void) fprintf( logFile, "\n%s: can't create grid extrema data file %s\n", programname, xyz_filename);
            (void) fprintf( logFile, "%s: SORRY!    unable to create the \".xyz\" file.\n\n", programname);
            (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
            exit(-1);
        } else {
            (void) fprintf( logFile, "\nCreating (AVS-readable) grid-coordinates extrema file : %s\n\n", xyz_filename);
        }
        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_NPTS:
        (void) sscanf( GPF_line, "%*s %d %d %d", &nelements[X], &nelements[Y], &nelements[Z]);
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
        num_grid_points_per_map = n1[X] * n1[Y] * n1[Z];
        percentdone = 100. / (double) n1[Z];
        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_SPACING:
        (void) sscanf( GPF_line, "%*s %lf", &spacing);
        (void) fprintf( logFile, "Grid Spacing :\t\t\t%.3lf Angstrom\n", spacing);
        (void) fprintf( logFile, "\n");
        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_GRIDCENTER:
        (void) sscanf( GPF_line, "%*s %s", token);
        if (equal( token, "auto", 4)) {
            for (i = 0;  i < XYZ;  i++) {
                center[i] = cmean[i];
            }
            (void) fprintf( logFile, "Grid maps will be centered on the center of mass.\n");
            (void) fprintf( logFile, "Coordinates of center of mass : (%.3lf, %.3lf, %.3lf)\n", center[X], center[Y], center[Z]);
        } else {
            (void) sscanf( GPF_line, "%*s %lf %lf %lf", &center[X], &center[Y], &center[Z]);
            (void) fprintf( logFile, "\nGrid maps will be centered on user-defined coordinates:\n\n\t\t(%.3lf, %.3lf, %.3lf)\n",  center[X], center[Y], center[Z]);
        }
        /* centering stuff... */
        for (ia = 0;  ia < num_receptor_atoms;  ia++) {
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
        (void) fprintf( logFile, "                   _______(%.1lf, %.1lf, %.1lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z]);
        (void) fprintf( logFile, "                  /|     /|\n");
        (void) fprintf( logFile, "                 / |    / |\n");
        (void) fprintf( logFile, "                /______/  |\n");
        (void) fprintf( logFile, "                |  |___|__| Midpoint = (%.1lf, %.1lf, %.1lf)\n", center[X], center[Y], center[Z]);
        (void) fprintf( logFile, "                |  /   |  /\n");
        (void) fprintf( logFile, "                | /    | /\n");
        (void) fprintf( logFile, "                |/_____|/\n");
        (void) fprintf( logFile, "(%.1lf, %.1lf, %.1lf)      \n\n", cgridmin[X], cgridmin[Y], cgridmin[Z]);
        for (i = 0;  i < XYZ;  i++) {
            (void) fprintf( logFile, "Grid map %c-dimension :\t\t%.1lf Angstroms\n", xyz[i], 2.*cext[i]);
        }
        (void) fprintf( logFile, "\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z]);
        (void) fprintf( logFile, "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cgridmin[X], cgridmin[Y], cgridmin[Z]);
        for (i = 0;  i < XYZ;  i++) {
            (void) fprintf(xyz_fileptr, "%.3lf %.3lf\n", cgridmin[i], cgridmax[i]);
        }
        (void) fclose(xyz_fileptr);
        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_LIGAND_TYPES:
        // Read in the list of atom types in the ligand.
        // GPF_line e.g.: "ligand_types N O A C HH NH"
        num_atom_maps = parsetypes(GPF_line, ligand_atom_types, MAX_ATOM_TYPES);
        for (i=0; i<num_atom_maps; i++) {
            strcpy(ligand_types[i], ligand_atom_types[i]);
#ifdef DEBUG
            (void) fprintf(logFile, "%d %s ->%s\n",i, ligand_atom_types[i], ligand_types[i]);
#endif
        }
        for(i=0; i<num_atom_maps; i++){
            found_parm = apm_find(ligand_types[i]);
            if (found_parm != NULL) {
                found_parm->map_index = i;
#ifdef DEBUG
                (void) fprintf(logFile, "found ligand type: %-6s%2d\n",
                        found_parm->autogrid_type,
                        found_parm->map_index );
#endif
            } 
            else {
                 // return error here
                (void) fprintf( logFile, "unknown ligand atom type %s\nadd parameters for it to the parameter library first!\n", ligand_atom_types[i]);
                 exit(-1);
             }
        };


        elecPE = num_atom_maps;
        dsolvPE = elecPE + 1;

        /* num_maps is the number of maps to be created:
         * the number of ligand atom types, plus 1 for the electrostatic map.
         * AutoDock can only read in MAX_MAPS maps, which must include
         * the ligand atom maps and electrostatic map */
        num_maps = num_atom_maps + 2;

        /* Check to see if there is enough memory to store these map objects */
        gridmap = (MapObject *)malloc(sizeof(MapObject) * num_maps);

        /* Check to see if the number of grid points requested will be
         * feasible; give warning if not enough memory. */
        if (num_grid_points_per_map != INIT_NUM_GRID_PTS) {
            dummy_map = (FloatOrDouble *)malloc(sizeof(FloatOrDouble) * (num_maps * num_grid_points_per_map));
            if (!dummy_map) {
                /* Too many maps requested */
                (void) sprintf(message, "\n%s: WARNING:  There will not be enough memory to store these grid maps in AutoDock; \ntry reducing the number of ligand atom types (you have %d including electrostatics) \nor reducing the size of the grid maps (you asked for %d x %d x %d grid points); \n or try running AutoDock on a machine with more RAM than this one.\n", programname, num_maps, n1[X], n1[Y], n1[Z]);
                (void) fprintf(stderr, "%s\n", message);
                (void) fprintf(logFile, "%s\n", message);
            } else {
                /* free up this memory right away; we were just testing to
                 * see if we had enough when we try to run AutoDock */
                free(dummy_map);
            }
        } else {
            (void) sprintf(message, "\n%s: ERROR:  You need to set the number of grid points using \"npts\" before setting the ligand atom types, using \"types\".\n", programname);
            (void) fprintf(stderr, "%s\n", message);
            (void) fprintf(logFile, "%s\n", message);
            (void) sprintf(message, "\n%s:  Unsuccessful completion.\n\n", programname);
            (void) fprintf(stderr, "%s\n", message);
            (void) fprintf(logFile, "%s\n", message);
            exit(-1);
        } /* ZZZZZZZZZZZZZZZZZ*/
        if (!gridmap) {
            (void) sprintf(message, "\n%s: ERROR:  Too many ligand atom types; there is not enough memory to create these maps.  Try using fewer atom types than %d.\n", programname, num_atom_maps);
            (void) fprintf(stderr, "%s\n", message);
            (void) fprintf(logFile, "%s\n", message);
            (void) sprintf(message, "\n%s:  Unsuccessful completion.\n\n", programname);
            (void) fprintf(stderr, "%s\n", message);
            (void) fprintf(logFile, "%s\n", message);
            exit(-1);
        }
        /*FIX THIS:
         *
         * !!!*/

        /* Check number of ligand atom maps requested * /
        if (num_maps > MAX_MAPS) {
            / * Too many maps requested * /
            (void) sprintf(message, "\n%s: ERROR:  The number of ligand atom maps requested (%d) plus 1 for the electrostatics map, is too many;  this program can only calculate %d maps.  Change the definition of \"MAX_MAPS\" in \"autocomm.h\".\n", programname, num_atom_maps, MAX_MAPS);
            (void) fprintf(stderr, "%s\n", message);
            (void) fprintf(logFile, "%s\n", message);
            (void) sprintf(message, "\n%s:  Unsuccessful completion.\n\n", programname);
            (void) fprintf(stderr, "%s\n", message);
            (void) fprintf(logFile, "%s\n", message);
            exit(-1);
        }
        */
        for (i = 0;  i < num_atom_maps;  i++) {
            gridmap[i].is_covalent = FALSE;
            gridmap[i].is_hbonder = FALSE;
            gridmap[i].map_index = i;
            found_parm = apm_find(ligand_types[i]);
            gridmap[i].atom_type = found_parm->map_index;
            gridmap[i].solpar_probe = found_parm->solpar;
            gridmap[i].vol_probe = found_parm->vol;
            gridmap[i].Rij = found_parm->Rij;
            gridmap[i].epsij = found_parm->epsij;
            gridmap[i].hbond = found_parm->hbond;
            gridmap[i].Rij_hb = found_parm->Rij_hb;
            gridmap[i].epsij_hb = found_parm->epsij_hb;
            if (gridmap[i].hbond>0){
                gridmap[i].is_hbonder=TRUE;}

#ifdef DEBUG
            (void) fprintf(logFile, " setting ij parms for map %d \n",i);
            (void) fprintf(logFile, "for gridmap[%d], type->%s,Rij->%6.4f, epsij->%6.4f, hbond->%d\n",i,found_parm->autogrid_type, gridmap[i].Rij, gridmap[i].epsij,gridmap[i].hbond);
#endif
            //for (j=0; j<NUM_RECEPTOR_TYPES;j++){
            for (j=0; j<receptor_types_ct;j++){
                /*SET THIS UP HERE!!!!*/
                found_parm = apm_find(receptor_types[j]);
                //strcpy(thisparm.autogrid_type, receptor_types[j]);
                //item.key = thisparm.autogrid_type;

                gridmap[i].nbp_r[j] = (gridmap[i].Rij + found_parm->Rij)/2.;
                gridmap[i].nbp_eps[j] = sqrt(gridmap[i].epsij * found_parm->epsij);
                /*apply the vdW forcefield parameter/weight here */
                gridmap[i].nbp_eps[j] *= FE_coeff_vdW;
                gridmap[i].xA[j] = 12;
                /*setup hbond dependent stuff*/
                gridmap[i].xB[j] = 6;
                gridmap[i].hbonder[j] = 0;
                if ((int)(gridmap[i].hbond)>2 &&
                        ((int)found_parm->hbond==1||(int)found_parm->hbond==2)){ /*AS,A1,A2 map vs DS,D1 probe*/
                    gridmap[i].xB[j] = 10;
                    gridmap[i].hbonder[j] = 1;
                    gridmap[i].is_hbonder = TRUE;
                    /*Rij and epsij for this hb interaction in
                     * parm_data.dat file as Rii and epsii for heavy atom
                     * hb factors*/
                    gridmap[i].nbp_r[j] = gridmap[i].Rij_hb;
                    gridmap[i].nbp_eps[j] = gridmap[i].epsij_hb;

                    /*apply the hbond forcefield parameter/weight here */
                    gridmap[i].nbp_eps[j] *= FE_coeff_hbond;
#ifdef DEBUG
                    (void) fprintf(logFile, "set %d-%d hb eps to %6.4f*%6.4f=%6.4f\n",i,j,gridmap[i].epsij_hb,found_parm->epsij_hb, gridmap[i].nbp_eps[j]);
#endif
                } else if (((int)gridmap[i].hbond==1||(int)gridmap[i].hbond==2) &&
                            ((int)found_parm->hbond>2)) { /*DS,D1 map vs AS,A1,A2 probe*/
                    gridmap[i].xB[j] = 10;
                    gridmap[i].hbonder[j] = 1;
                    gridmap[i].is_hbonder = TRUE;
                    /*Rij and epsij for this hb interaction in
                     * parm_data.dat file as Rii and epsii for heavy atom
                     * hb factors*/
                    gridmap[i].nbp_r[j] = found_parm->Rij_hb;
                    gridmap[i].nbp_eps[j] = found_parm->epsij_hb;

                    /*apply the hbond forcefield parameter/weight here */
                    gridmap[i].nbp_eps[j] *= FE_coeff_hbond;
#ifdef DEBUG
                    (void) fprintf(logFile, "2: set %d-%d hb eps to %6.4f*%6.4f=%6.4f\n",i,j,gridmap[i].epsij_hb,found_parm->epsij_hb, gridmap[i].nbp_eps[j]);
#endif
                } 
#ifdef DEBUG
                (void) fprintf(logFile, "vs receptor_type[%d]:type->%s, hbond->%d ",j,found_parm->autogrid_type, (int)found_parm->hbond);
                (void) fprintf(logFile, "nbp_r->%6.4f, nbp_eps->%6.4f,xB=%d,hbonder=%d\n",gridmap[i].nbp_r[j], gridmap[i].nbp_eps[j],gridmap[i].xB[j], gridmap[i].hbonder[j]);
#endif
              //};
            }; /*initialize energy parms for each possible receptor type*/
        } /*for each map*/
        (void) fprintf( logFile, "\nAtom names for ligand_types 1-%d and for ligand-atom affinity grid maps :\n", num_atom_maps);
        for (i = 0;  i < num_atom_maps;  i++) {
            (void) fprintf( logFile, "\t\t\t%d->%s\n", gridmap[i].map_index, ligand_types[i]);
            /*FIX THIS!!!*/
            /*if (gridmap[i].atom_type == COVALENTTYPE) {
              gridmap[i].is_covalent = TRUE;
              (void) fprintf( logFile, "\nAtom type number %d will be used to calculate a covalent affinity grid map\n\n", i + 1);
            }*/
        }
        (void) fprintf( logFile, " respectively.\n");
        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_RECEPTOR_TYPES:
        // Read in the list of atom types in the receptor.
        // GPF_line e.g.: "receptor_types N O A C HH NH"
        receptor_types_ct = parsetypes(GPF_line, receptor_atom_types,  MAX_ATOM_TYPES);

#ifdef DEBUG
        printf("receptor_types_ct=%d\n",receptor_types_ct);
#endif
        for(i=0; i<receptor_types_ct; i++){
            strcpy(receptor_types[i], receptor_atom_types[i]);
#ifdef DEBUG
            printf("%d %s  ->%s\n",i, receptor_atom_types[i],  receptor_types[i]);
#endif
        }
        for (i=0; i<receptor_types_ct; i++) {
            found_parm = apm_find(receptor_atom_types[i]);
            if (found_parm != NULL){
                found_parm->rec_index = i;
            } else {
                fprintf(stderr, "unknown receptor type:%s\n add parameters for it to the parameter library first!\n", receptor_atom_types[i]);
                exit(-1);
            }
        };
        /*at this point set up hydrogen, carbon, oxygen and nitrogen*/
        hydrogen = get_rec_index("HD");
        nonHB_hydrogen = get_rec_index("H");
        carbon = get_rec_index("C");
        arom_carbon = get_rec_index("A");
        oxygen = get_rec_index("OA");
        nitrogen = get_rec_index("NA");
        nonHB_nitrogen = get_rec_index("N");
        sulphur = get_rec_index("SA");
        nonHB_sulphur = get_rec_index("S");

#ifdef DEBUG
        printf("assigned receptor types:arom_carbon->%d, hydrogen->%d,nonHB_hydrogen->%d, carbon->%d, oxygen->%d, nitrogen->%d\n, nonHB_nitrogen->%d, sulphur->%d, nonHB_sulphur->%d\n",arom_carbon,hydrogen, nonHB_hydrogen, carbon,oxygen, nitrogen, nonHB_nitrogen, sulphur, nonHB_sulphur);
#endif

        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_SOL_PAR:  //THIS IS OBSOLETE!!!
        /*
        ** Read volume and solvation parameter for probe:
        */
        (void) sscanf( GPF_line, "%*s %s %lf %lf", thisparm.autogrid_type, &temp_vol, &temp_solpar );
        found_parm = apm_find(thisparm.autogrid_type);
        if (found_parm != NULL) {
            found_parm->vol = temp_vol;
            found_parm->solpar = temp_solpar;
            i = found_parm->map_index;
            if (i>=0){
                /*DON'T!!!*/
                /*convert cal/molA^3 to kcal/molA^3 */
                /*gridmap[i].solpar_probe = temp_solpar * 0.001;*/
                gridmap[i].solpar_probe = temp_solpar ;
                (void) fprintf( logFile, "\nProbe %s solvation parameters: \n\n\tatomic fragmental volume: %.2f A^3\n\tatomic solvation parameter: %.4f cal/mol A^3\n\n", found_parm->autogrid_type, found_parm->vol,found_parm->solpar);
                (void) fflush( logFile);
            }
        } else {
            (void) fprintf( logFile, "%s key not found\n", thisparm.autogrid_type);
        };
        (void) fflush( logFile);
        break; /* end solvation parameter */

/******************************************************************************/

    /*case GPF_CONSTANT:*/
        /*break;*/

/******************************************************************************/

    case GPF_MAP: 
        /*  */
        /* The variable "map_index" is the 0-based index of the ligand atom type
         * we are calculating a map for. 
         * If the "types" line was CNOSH, there would be 5 ligand atom maps to calculate,
         * and since "map_index" is initialized to -1, map_index will increment
         * each time there is a "map" keyword in the GPF.  The value of
         * map_index should therefore go from 0 to 4 for each "map" keyword.  
         * In this example, num_atom_maps would be 5, and num_atom_maps-1 would be 
         * 4, so if map_index is > 4, there is something wrong in the number of
         * "map" keywords. */
        ++map_index;
        if (map_index > num_atom_maps - 1) {
             (void) sprintf(message, "\n%s: ERROR:  Too many \"map\" keywords (%d);  the \"types\" command declares only %d maps.\nRemove a \"map\" keyword from the GPF.\n", programname, map_index + 1, num_atom_maps);
            (void) fprintf(stderr, "%s\n", message);
            (void) fprintf(logFile, "%s\n", message);
            (void) sprintf(message, "\n%s:  Unsuccessful completion.\n\n", programname);
            (void) fprintf(stderr, "%s\n", message);
            (void) fprintf(logFile, "%s\n", message);
            exit(-1);
        }
        /* Read in the filename for this grid map */ /* GPF_MAP */
        (void) sscanf( GPF_line, "%*s %s", gridmap[map_index].map_filename);
        if ( (gridmap[map_index].map_fileptr = ag_fopen( gridmap[map_index].map_filename, "w")) == NULL ) {
            (void) fprintf( stderr, "\n%s: can't open grid map \"%s\" for writing.\n", programname, gridmap[map_index].map_filename);
            (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
            (void) fprintf( logFile, "\n%s: can't open grid map \"%s\" for writing.\n", programname, gridmap[map_index].map_filename);
            (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
            exit(-1);
        }
        (void) fprintf( logFile, "\n\nOutput File %d   %s\n\n", (map_index + 1), gridmap[map_index].map_filename);
        (void) fflush( logFile);

        break;

/******************************************************************************/
    case GPF_ELECMAP:
        (void) sscanf( GPF_line, "%*s %s", gridmap[elecPE].map_filename);
        if ( (gridmap[elecPE].map_fileptr = ag_fopen( gridmap[elecPE].map_filename, "w" )) == NULL){
            (void) fprintf( stderr, "\n%s: can't open grid map \"%s\" for writing.\n", programname, gridmap[elecPE].map_filename);
            (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
            (void) fprintf( logFile, "\n%s: can't open grid map \"%s\" for writing.\n", programname, gridmap[elecPE].map_filename);
            (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
            exit(-1);
        }
        (void) fprintf( logFile, "\nElectrostatics File: %s\n\n", gridmap[elecPE].map_filename);
        break;

/******************************************************************************/
    case GPF_DSOLVMAP:
        (void) sscanf( GPF_line, "%*s %s", gridmap[dsolvPE].map_filename);
        if ( (gridmap[dsolvPE].map_fileptr = ag_fopen( gridmap[dsolvPE].map_filename, "w" )) == NULL){
            (void) fprintf( stderr, "\n%s: can't open grid map \"%s\" for writing.\n", programname, gridmap[dsolvPE].map_filename);
            (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
            (void) fprintf( logFile, "\n%s: can't open grid map \"%s\" for writing.\n", programname, gridmap[dsolvPE].map_filename);
            (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
            exit(-1);
        }
        (void) fprintf( logFile, "\nDesolvation File: %s\n\n", gridmap[dsolvPE].map_filename);
        break;
/******************************************************************************/

    case GPF_COVALENTMAP:
        (void) sscanf( GPF_line, "%*s %lf %lf %lf %lf %lf", &covhalfwidth, &covbarrier, &(covpos[X]), &(covpos[Y]), &(covpos[Z]));
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
        (void) sscanf( GPF_line, "%*s %lf", &r_smooth);
        (void) fprintf( logFile, "\nPotentials will be smoothed by: %.3lf\n\n", r_smooth);
        /* Angstrom is divided by A_DIVISOR in look-up table. */
        /* Typical value of r_smooth is 0.5 Angstroms */
        /* so i_smooth = 0.5 * 100. / 2 = 25 */
        i_smooth = (int) (r_smooth*A_DIVISOR/2.);
        break;

/******************************************************************************/

    case GPF_QASP:
        (void) sscanf( GPF_line, "%*s %lf", &solpar_q);
        (void) fprintf( logFile, "\nCharge component of the atomic solvation parameter: %.3lf\n\n", solpar_q);
        /* Typical value of solpar_q is 0.001118 */
        break;

/******************************************************************************/
    case GPF_DIEL:
        (void) sscanf( GPF_line, "%*s %lf", &diel);
        if (diel < 0.) {
            /* negative... */
            (void) fprintf( logFile, "\nUsing *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.\n\n");
            dddiel = TRUE;
            /*____________________________________________________________________________
             * Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
             *____________________________________________________________________________*/
            lambda   = 0.003627;
            epsilon0 = 78.4;
            A        = (-8.5525);
            B        = epsilon0 - A;
            rk       = 7.7839;
            lambda_B = (-lambda) * B;
            for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
                epsilon[indx_r] = A + B / (1. + rk*exp(lambda_B*angstrom(indx_r)));
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
                /*3/9/2005: do NOT apply it here*/
                /*epsilon[i] *= -diel;*/
                /* Really, LHS should be "inv_epsilon[i]", but this way saves memory... */
            }
        } else {
            /* positive or zero... */
            dddiel = FALSE;
            if (diel <= APPROX_ZERO) {
                diel = 40.;
            }
            (void) fprintf( logFile, "Using a *constant* dielectric of:  %.2f\n", diel);
            invdielcal = 332. / diel;
        }
        break;

/******************************************************************************/

    case GPF_FMAP:
        (void) sscanf( GPF_line, "%*s %s", floating_grid_filename);
        if ( (floating_grid_fileptr = ag_fopen( floating_grid_filename, "w" )) == NULL) {
            (void) fprintf( stderr, "\n%s: can't open grid map \"%s\" for writing.\n", programname, floating_grid_filename);
            (void) fprintf( stderr, "\n%s: Unsuccessful completion.\n\n", programname);
            (void) fprintf( logFile, "\n%s: can't open grid map \"%s\" for writing.\n", programname, floating_grid_filename);
            (void) fprintf( logFile, "\n%s: Unsuccessful completion.\n\n", programname);
            exit(-1);
        }
        (void) fprintf( logFile, "\nFloating Grid file name = %s\n", floating_grid_filename);
        ++num_maps;
        floating_grid = TRUE;
        break;

/******************************************************************************/

    case GPF_PARAM_FILE:
        /*open and read the parm_data.dat file*/

        (void) sscanf( GPF_line, "%*s %s ", param_filename);
        /*flag that data from a parameter file are available*/

        if ((dataFile = ag_fopen(param_filename, "r"))==NULL){
             fprintf(stderr,"Sorry, I can't find or open %s\n",param_filename);
             exit(911);
        }

        while( fgets(dataline, sizeof(dataline), dataFile)!=NULL){
            char * commentbegin;
            int nfields;
            /*terminate line at commentbegin*/
            
            commentbegin = strchr(dataline, '#');
            if (commentbegin!=NULL) *commentbegin = '\0';
            /*NB: scanf doesn't try to write missing fields*/
            /*WISH LIST:
             * read autogrid_type into a temporary string and test 
             * that it's not too long*/
            /*check for FF_coeff lines*/
            if (strncmp(dataline, "FE_coeff_vdW", 12)==0) {
                sscanf(dataline, "%*s %lf", &FE_coeff_vdW);
            } else if (strncmp(dataline, "FE_coeff_hbond", 14)==0) {
                sscanf(dataline, "%*s %lf", &FE_coeff_hbond);
            } else if (strncmp(dataline, "FE_coeff_estat", 14)==0) {
                sscanf(dataline, "%*s %lf", &FE_coeff_estat);
            } else if (strncmp(dataline, "FE_coeff_desolv", 15)==0) {
                sscanf(dataline, "%*s %lf", &FE_coeff_desolv);
            } else if (strncmp(dataline, "FE_coeff_tors", 13)==0) {
                sscanf(dataline, "%*s %lf", &FE_coeff_tors);
            } else { /*skip 'atom_parm' at start of line*/
            nfields = sscanf(dataline, "%*s %s %lf %lf %lf %lf %lf %lf %d %d %d %d",
                    thisparm.autogrid_type,
                    &thisparm.Rij,
                    &thisparm.epsij,
                    &thisparm.vol,
                    &thisparm.solpar,
                    &thisparm.Rij_hb,
                    &thisparm.epsij_hb,
                    &thisparm.hbond,
                    &thisparm.rec_index,
                    &thisparm.map_index,
                    &thisparm.bond_index);
                    
            if (nfields<2) continue; /*skip lines without enough info*/
            apm_enter(thisparm.autogrid_type, thisparm);
#ifdef DEBUG
            printf("%-6s%7.2f%8.3f%8.3f%8.3f %8.3f%8.3f %d %d %d\n",thisparm.autogrid_type, thisparm.Rij, thisparm.epsij, thisparm.vol, thisparm.solpar,  thisparm.Rij_hb, thisparm.epsij_hb, thisparm.hbond, thisparm.rec_index, thisparm.map_index);
#endif
           } /*end of scan a new parameter line*/
        }
        break;



/******************************************************************************/

    default:
        break;

/******************************************************************************/

    } /* second switch */
} /* while */

(void) fprintf( logFile, "\n>>> Closing the grid parameter file (GPF)... <<<\n\n");
(void) fprintf( logFile, UnderLine);
(void) fclose( GPF_fileptr);

if ( ! floating_grid ) {
    (void) fprintf( logFile, "\n\nNo Floating Grid was requested.\n");
}

(void) fprintf( AVS_fld_fileptr, "# AVS field file\n#\n");
(void) fprintf( AVS_fld_fileptr, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
(void) fprintf( AVS_fld_fileptr, "# Created by %s.\n#\n", programname);
(void) fprintf( AVS_fld_fileptr, "#SPACING %.3f\n", (float) spacing);
(void) fprintf( AVS_fld_fileptr, "#NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
(void) fprintf( AVS_fld_fileptr, "#CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
(void) fprintf( AVS_fld_fileptr, "#MACROMOLECULE %s\n", receptor_filename);
(void) fprintf( AVS_fld_fileptr, "#GRID_PARAMETER_FILE %s\n#\n", GPF_filename);
(void) fprintf( AVS_fld_fileptr, "ndim=3\t\t\t# number of dimensions in the field\n");
(void) fprintf( AVS_fld_fileptr, "dim1=%d\t\t\t# number of x-elements\n", n1[X]);
(void) fprintf( AVS_fld_fileptr, "dim2=%d\t\t\t# number of y-elements\n", n1[Y]);
(void) fprintf( AVS_fld_fileptr, "dim3=%d\t\t\t# number of z-elements\n", n1[Z]);
(void) fprintf( AVS_fld_fileptr, "nspace=3\t\t# number of physical coordinates per point\n");
(void) fprintf( AVS_fld_fileptr, "veclen=%d\t\t# number of affinity values at each point\n", num_maps);
(void) fprintf( AVS_fld_fileptr, "data=float\t\t# data type (byte, integer, float, double)\n");
(void) fprintf( AVS_fld_fileptr, "field=uniform\t\t# field type (uniform, rectilinear, irregular)\n");
for (i = 0;  i < XYZ;  i++) {
    (void) fprintf( AVS_fld_fileptr, "coord %d file=%s filetype=ascii offset=%d\n", (i + 1), xyz_filename, (i*2));
}
for (i = 0;  i < num_atom_maps;  i++) {
    (void) fprintf( AVS_fld_fileptr, "label=%s-affinity\t# component label for variable %d\n", gridmap[i].type, (i + 1));
} /* i */
(void) fprintf( AVS_fld_fileptr, "label=Electrostatics\t# component label for variable %d\n", num_maps-2);
(void) fprintf( AVS_fld_fileptr, "label=Desolvation\t# component label for variable %d\n", num_maps-1);
if (floating_grid) {
    (void) fprintf( AVS_fld_fileptr, "label=Floating_Grid\t# component label for variable %d\n", num_maps);
}
(void) fprintf( AVS_fld_fileptr, "#\n# location of affinity grid files and how to read them\n#\n");
for (i = 0;  i < num_atom_maps;  i++) {
    (void) fprintf( AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", (i + 1), gridmap[i].map_filename);
}
(void) fprintf( AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_atom_maps + 1, gridmap[elecPE].map_filename);
(void) fprintf( AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_atom_maps + 2, gridmap[dsolvPE].map_filename);
if (floating_grid) {
    (void) fprintf( AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_maps, floating_grid_filename);
}
(void) fclose( AVS_fld_fileptr);


#ifdef BOINCCOMPOUND
 boinc_fraction_done(0.1);
#endif

/**************************************************
 * do the map stuff here: 
 * set up xA, xB, npb_r, npb_eps and hbonder 
 * before this pt
 **************************************************/
for (ia=0; ia<num_atom_maps; ia++){
    if (gridmap[ia].is_covalent == FALSE) {
        /* i is the index of the receptor atom type, that 
         * the ia type ligand probe will interact with. */ /* GPF_MAP */
#ifdef DEBUG
        printf("receptor_types_ct=%d\n", receptor_types_ct);
#endif
        for (i = 0;  i < receptor_types_ct;  i++) {
            /*for each receptor_type*/
            xA = gridmap[ia].xA[i];
            xB = gridmap[ia].xB[i];
            Rij = gridmap[ia].nbp_r[i];
            epsij = gridmap[ia].nbp_eps[i];
#ifdef DEBUG
            printf("%d-%d-built xA=%d, xB=%d, npb_r=%6.3lf, nbp_eps=%10.8f for %s\n",ia,i,xA, xB, Rij,epsij,gridmap[ia].type);
#endif
            /*for each receptor_type get its parms and fill in tables*/
            cA = (tmpconst = epsij / (double)(xA - xB)) * pow( Rij, (double)xA ) * (double)xB;
            cB = tmpconst * pow( Rij, (double)xB ) * (double)xA;
            /*printf("tmpconst = %6.4f, cA = %6.4f, cB = %6.4f\n",tmpconst, cA, cB);*/
            dxA = (double) xA;
            dxB = (double) xB;

        (void) fprintf( logFile, "\n             %9.1lf       %9.1lf \n", cA, cB);
        (void) fprintf( logFile, "    E    =  -----------  -  -----------\n");
        (void) fprintf( logFile, "     %s, %s         %2d              %2d\n", ligand_types[ia], receptor_types[i], xA, xB);
        (void) fprintf( logFile, "                r               r \n\n");
        /* loop over distance index, indx_r, from 0 to MAX_DIST */ /* GPF_MAP */
        (void) fprintf( logFile, "calculating energy table %d->%s, %d->%s\n",ia, gridmap[ia].type,i, receptor_types[i]);
        for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
            r  = angstrom(indx_r);
            rA = pow( r, dxA);
            rB = pow( r, dxB);
            energy_lookup[i][indx_r][ia] = min(EINTCLAMP, (cA/rA - cB/rB));
            } /*for each distance*/ 
            energy_lookup[i][0][ia]    = EINTCLAMP;
            energy_lookup[i][MD_1][ia] = 0.;
        /*PRINT OUT INITIAL VALUES before smoothing here
        (void) fprintf( logFile, "before smoothing\n  r ");
        for (iat = 0;  iat < receptor_types_ct;  iat++) {
            (void) fprintf( logFile, "    %s    ", receptor_types[iat]);
        } 
        (void) fprintf( logFile, "\n ___");
        for (iat = 0;  iat < receptor_types_ct;  iat++) {
            (void) fprintf( logFile, " ________");
        }
        (void) fprintf( logFile, "\n");
        for (j = 0;  j <= 500;  j += 10) {
            (void) fprintf( logFile, "%4.1lf", angstrom(j));
            for (iat = 0;  iat < receptor_types_ct;  iat++) {
                (void) fprintf( logFile, (energy_lookup[iat][j][ia]<100000.)?"%9.2lf":"%9.2lg", energy_lookup[iat][j][ia]);
            } 
            (void) fprintf( logFile, "\n");
        } 
        (void) fprintf( logFile, "\n");*/




            /* smooth with min function */ /* GPF_MAP */
            if (i_smooth > 0) {
                for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
                    energy_smooth[indx_r] = 100000.;
                    for (j = max(0, indx_r - i_smooth);  j < min(MAX_DIST, indx_r + i_smooth);  j++) {
                      energy_smooth[indx_r] = min(energy_smooth[indx_r], energy_lookup[i][j][ia]);
                    }
                }
                for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
                    energy_lookup[i][indx_r][ia] = energy_smooth[indx_r];
                }
                } /* end smoothing */
            } /* for i in receptor types: build energy table for this map */

       /*
        * Print out a table, of distance versus energy...
        */ /* GPF_MAP */
        (void) fprintf( logFile, "after smoothing\n  r ");
        for (iat = 0;  iat < receptor_types_ct;  iat++) {
            (void) fprintf( logFile, "    %s    ", receptor_types[iat]);
            /*(void) fprintf( logFile, "    %c    ", receptor_atom_type_string[iat]);*/
        } /* iat */
        (void) fprintf( logFile, "\n ___");
        for (iat = 0;  iat < receptor_types_ct;  iat++) {
            (void) fprintf( logFile, " ________");
        } /* iat */
        (void) fprintf( logFile, "\n");
        for (j = 0;  j <= 500;  j += 10) {
            (void) fprintf( logFile, "%4.1lf", angstrom(j));
            for (iat = 0;  iat < receptor_types_ct;  iat++) {
                (void) fprintf( logFile, (energy_lookup[iat][j][ia]<100000.)?"%9.2lf":"%9.2lg", energy_lookup[iat][j][ia]);
            } /* iat */
            (void) fprintf( logFile, "\n");
        } /* j */
        (void) fprintf( logFile, "\n");
    } else {
        /* parsing for intnbp not needed for covalent maps */
        (void) fprintf( logFile, "\nAny internal non-bonded parameters will be ignored for this map, since this is a covalent map.\n");
    } /*end of else parsing intnbp*/
} /*end of loop over all the maps*/
/* exponential function for receptor and ligand desolvation */
/* note: the solvation term will not be smoothed */
sigma = 3.6;
minus_inv_two_sigma_sqd = -1. / (2. * sigma * sigma);
for (indx_r = 1;  indx_r < MAX_DIST;  indx_r++) {
     r  = angstrom(indx_r);
     /* sol_fn[indx_r] = exp(-sq(r)/(2.*sigma*sigma)); */
     sol_fn[indx_r] = exp( sq(r) * minus_inv_two_sigma_sqd);
     sol_fn[indx_r] *= FE_coeff_desolv;
}

/**************************************************
 * Loop over all RECEPTOR atoms to
 * calculate bond vectors for directional H-bonds
 **************************************************/

/*7:CHANGE HERE: scan the 'map_index' from the input*/
for (ia=0; ia<num_receptor_atoms; ia++) {  /*** ia = i_receptor_atom_a ***/

    disorder[ia] = FALSE;  /* initialize disorder flag. */
    warned = 'F';

    /*
     * Set scan limits looking for bonded atoms
     */
    from = max(ia-20, 0);
    to   = min(ia + 20, num_receptor_atoms-1);

    /*
     * If 'ia' is a hydrogen atom, it could be a
     * RECEPTOR hydrogen-BOND DONOR,
     */
    /*8:CHANGE HERE: fix the atom_type vs atom_types problem in following*/
    if ((int)hbond[ia] == 2) { /*D1 hydrogen bond donor*/

        for ( ib = from; ib <= to; ib++) {       /*** ib = i_receptor_atom_b ***/
            /*
             * =>  NH-> or OH->
             */
            /*if ((atom_type[ib] == nitrogen) || (atom_type[ib]==nonHB_nitrogen) ||(atom_type[ib] == oxygen)||(atom_type[ib] == sulphur)||(atom_type[ib]==nonHB_sulphur)) {*/
    
                /*
                 * Calculate the square of the N-H or O-H bond distance, rd2,
                 *                            ib-ia  ib-ia
                 */
                for (i = 0;  i < XYZ;  i++) {
                    d[i] = coord[ia][i] - coord[ib][i];
                }
                rd2 = sq( d[X] ) + sq( d[Y] ) + sq( d[Z]);
                /*
                 * If ia & ib are less than 1.3 A apart -- they are covalently bonded,
                 */
                if (rd2 < 1.90) { /*INCREASED for H-S bonds*/
                    if (rd2 < APPROX_ZERO) {
                        if (rd2 == 0.) {
                            (void) fprintf (stderr, "WARNING! While calculating an H-O or H-N bond vector...\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                            (void) fprintf (logFile, "WARNING! While calculating an H-O or H-N bond vector...\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1./sqrt(rd2);
                    /*
                     * N-H: Set exponent rexp to 2 for m/m H-atom,
                     */
                    /*if (atom_type[ib] == nitrogen) rexp[ia] = 2;*/
                    if ((atom_type[ib] != oxygen)&&(atom_type[ib] != sulphur)) rexp[ia] = 2;

                    /*
                     * O-H: Set exponent rexp to 4 for m/m H-atom,
                     * and flag disordered hydroxyls
                     */
                    if ((atom_type[ib] == oxygen)||(atom_type[ib] == sulphur)) {
                        rexp[ia] = 4;
                        if (disorder_h == TRUE) disorder[ia] = TRUE;
                    }
                    /*
                     * Normalize the vector from ib to ia, N->H or O->H...
                     */
                    for (i = 0;  i < XYZ;  i++) {
                        rvector[ia][i] = d[i] * inv_rd;
                    }
                    /*
                     * First O-H/N-H H-bond-donor found; Go on to next atom,
                     */
                    break;
                } /* Found covalent bond. */
            /*}  Found NH or OH in receptor. */
        } /* Finished scanning for the NH or OH in receptor. */

    /*
     * If 'ia' is an Oxygen atom, it could be a
     * RECEPTOR H_BOND ACCEPTOR,
     */

    } else if (hbond[ia] == 5) { /*A2*/
        /*
         * Scan from at most, (ia-20)th m/m atom, or ia-th (if ia<20)
         *        to (ia + 5)th m/m-atom
         * determine number of atoms bonded to the oxygen
         */
        nbond = 0;
        for ( ib = from; ib <= to; ib++) {
            if ( ib != ia ) {
                rd2 = 0.;

                for (i = 0;  i < XYZ;  i++) {
                    dc[i] = coord[ia][i] - coord[ib][i];
                    rd2 += sq( dc[i]);
                }

                /*
                    for (i = 0;  i < XYZ;  i++) {
                        rd2 += sq(coord[ia][i] - coord[ib][i]);
                    }
                */
                if (((rd2 < 3.61) && ((atom_type[ib] != hydrogen)&&(atom_type[ib]!=nonHB_hydrogen))) ||
                    ((rd2 < 1.69) && ((atom_type[ib] == hydrogen)||(atom_type[ib]==nonHB_hydrogen)))) {
                    if (nbond == 2) {
                        (void) fprintf( logFile, "WARNING! h-bonding atom with three bonded atoms, atom serial number %d\n", ia + 1);
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
            (void) fprintf( logFile, "WARNING! oxygen with no bonded atoms, atom serial number %d,atom_type->%d\n", ia + 1, atom_type[ia]);
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
                    (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                    (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
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
                        rd2 += sq( dc[i]);
                    }
                    if (((rd2 < 2.89) && (atom_type[i2] != hydrogen)) ||
                        ((rd2 < 1.69) && (atom_type[i2] == hydrogen))) {

                        /* found one */
                        /* d[i] vector from carbon to second atom */
                        rd2 = 0.;
                        for (i = 0;  i < XYZ;  i++) {
                            d[i] = coord[i2][i]-coord[i1][i];
                            rd2 += sq( d[i]);
                        }
                        if (rd2 < APPROX_ZERO) {
                            if ((rd2 == 0.) && (warned == 'F')) {
                                (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                                (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                                warned = 'T';
                            }
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1./sqrt(rd2);
                        for (i = 0;  i < XYZ;  i++) {
                            d[i] *= inv_rd;
                        }

                        /* C=O cross C-X gives the lone pair plane normal */
                        rvector2[ia][0] = rvector[ia][1]*d[2] - rvector[ia][2]*d[1];
                        rvector2[ia][1] = rvector[ia][2]*d[0] - rvector[ia][0]*d[2];
                        rvector2[ia][2] = rvector[ia][0]*d[1] - rvector[ia][1]*d[0];
                        rd2 = 0.;
                        for (i = 0;  i < XYZ;  i++) {
                            rd2 += sq(rvector2[ia][i]);
                        }
                        if (rd2 < APPROX_ZERO) {
                            if ((rd2 == 0.) && (warned == 'F')) {
                                (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                                (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
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

            if ( ((atom_type[i1] == hydrogen) || (atom_type[i2] == hydrogen))
                && (atom_type[i1] != atom_type[i2]) && (disorder_h == TRUE) )  {

                if ((atom_type[i1] == carbon)||(atom_type[i1] == arom_carbon)) ib = i1;
                if ((atom_type[i2] == carbon)||(atom_type[i1] == arom_carbon)) ib = i2;
                disorder[ia] = TRUE;
                rd2 = 0.;
                for (i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i] - coord[ib][i];
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                        (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
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
                    rvector2[ia][i] = coord[i2][i] - coord[i1][i];
                    rd2 += sq(rvector2[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                        (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
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
                rdot = 0.;
                for (i = 0;  i < XYZ;  i++) {
                    rdot += (coord[ia][i] - coord[i1][i]) * rvector2[ia][i] ;
                }
                rd2 = 0.;
                for (i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i] - ( (rdot*rvector2[ia][i]) + coord[i1][i] ) ;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                        (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
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
        /* NEW Directional N Acceptor */
    } else if (hbond[ia] == 4) {/*A1*/
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
		if (((rd2 < 2.89) && ((atom_type[ib] != hydrogen)&&(atom_type[ib]!=nonHB_hydrogen))) || 
		    ((rd2 < 1.69) && ((atom_type[ib] == hydrogen)||(atom_type[ib]==nonHB_hydrogen)))) {
	        if (nbond == 2) {
			    nbond = 3;
			    i3 = ib;
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
	    (void) fprintf( logFile, "WARNING! nitrogen with no bonded atoms, atom serial number %d\n",ia);
	}

	/* one bond: Azide Nitrogen :N=C-X */

	if (nbond == 1) {

	    /* calculate normalized N=C bond vector rvector[ia][] */

	    rd2 = 0.;
	    for (i = 0;  i < XYZ;  i++) {
		    rvector[ia][i] = coord[ia][i]-coord[i1][i];
		    rd2 += sq(rvector[ia][i]);
	    }
	    if (rd2 < APPROX_ZERO) {
		    if ((rd2 == 0.) && (warned == 'F')) {
		        (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
		        (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
		        warned = 'T';
		    }
		    rd2 = APPROX_ZERO;
	    }
	    inv_rd = 1./sqrt(rd2);
	    for (i = 0;  i < XYZ;  i++) {
		    rvector[ia][i] *= inv_rd;
	    }
	} /* endif nbond==1 */

	/* two bonds: X1-N=X2 */
	if (nbond == 2) {
		/* normalized vector from Nitrogen to midpoint between X1 and X2 */

		rd2 = 0.;
		for (i = 0;  i < XYZ;  i++) {
		    rvector[ia][i] = coord[ia][i]-(coord[i2][i]+coord[i1][i])/2.;
		    rd2 += sq(rvector[ia][i]);
		}
		if (rd2 < APPROX_ZERO) {
		    if ((rd2 == 0.) && (warned == 'F')) {
			(void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
			(void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
			warned = 'T';
		    }
		    rd2 = APPROX_ZERO;
		}
		inv_rd = 1./sqrt(rd2);
		for (i = 0;  i < XYZ;  i++) {
		    rvector[ia][i] *= inv_rd;
		}

	}  /* end two bonds for nitrogen*/

	/* three bonds: X1,X2,X3 */
	if (nbond == 3) {
		/* normalized vector from Nitrogen to midpoint between X1, X2, and X3 */

		rd2 = 0.;
		for (i = 0;  i < XYZ;  i++) {
		    rvector[ia][i] = coord[ia][i]-(coord[i1][i]+coord[i2][i]+coord[i3][i])/3.;
		    rd2 += sq(rvector[ia][i]);
		}
		if (rd2 < APPROX_ZERO) {
		    if ((rd2 == 0.) && (warned == 'F')) {
			(void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
			(void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
			warned = 'T';
		    }
		    rd2 = APPROX_ZERO;
		}
		inv_rd = 1./sqrt(rd2);
		for (i = 0;  i < XYZ;  i++) {
		    rvector[ia][i] *= inv_rd;
		}

	}  /* end three bonds for Nitrogen */
    /* endNEW directional N Acceptor */


    } /* end test for atom type */

} /* Do Next receptor atom... */

/********************************************
 * End bond vector loop
 ********************************************/
for (k = 0;  k < num_atom_maps + 1;  k++) {
    gridmap[k].energy_max = (double)-BIG;
    gridmap[k].energy_min = (double)BIG;
}

(void) fprintf( logFile, "Beginning grid calculations.\n");
(void) fprintf( logFile, "\nCalculating %d grids over %d elements, around %d receptor atoms.\n\n", num_maps, num_grid_points_per_map, num_receptor_atoms);
(void) fflush( logFile);

/*____________________________________________________________________________
 * Write out the  correct grid_data '.fld' file_name at the  head of each map
 * file, to avoid centering errors in subsequent dockings...
 * AutoDock can then  check to see  if the  center of each  map  matches that
 * specified in its parameter file...
 *____________________________________________________________________________*/
/*change num_atom_maps +1 to num_atom_maps + 2 for new dsolvPE map*/
for (k = 0;  k < num_atom_maps+2;  k++) {
    (void) fprintf( gridmap[k].map_fileptr, "GRID_PARAMETER_FILE %s\n", GPF_filename);
    (void) fprintf( gridmap[k].map_fileptr, "GRID_DATA_FILE %s\n", AVS_fld_filename);
    (void) fprintf( gridmap[k].map_fileptr, "MACROMOLECULE %s\n", receptor_filename);
    (void) fprintf( gridmap[k].map_fileptr, "SPACING %.3lf\n", spacing);
    (void) fprintf( gridmap[k].map_fileptr, "NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
    (void) fprintf( gridmap[k].map_fileptr, "CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
}
if (floating_grid) {
    (void) fprintf( floating_grid_fileptr, "GRID_PARAMETER_FILE %s\n", GPF_filename);
    (void) fprintf( floating_grid_fileptr, "GRID_DATA_FILE %s\n", AVS_fld_filename);
    (void) fprintf( floating_grid_fileptr, "MACROMOLECULE %s\n", receptor_filename);
    (void) fprintf( floating_grid_fileptr, "SPACING %.3lf\n", spacing);
    (void) fprintf( floating_grid_fileptr, "NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
    (void) fprintf( floating_grid_fileptr, "CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
}

(void) fprintf( logFile, "                    Percent   Estimated Time  Time/this plane\n");
(void) fprintf( logFile, "XY-plane  Z-coord   Done      Remaining       Real, User, System\n");
(void) fprintf( logFile, "            /Ang              /sec            /sec\n");
(void) fprintf( logFile, "________  ________  ________  ______________  __________________________\n\n");

/*
 * Iterate over all grid points, Z( Y ( X ) ) (X is fastest)...
 */

ic = 0;

ctr = 0;
for (icoord[Z] = -ne[Z]; icoord[Z] <= ne[Z]; icoord[Z]++) {
    /*
     *  c[0:2] contains the current grid point.
     */
    c[Z] = ((double)icoord[Z]) * spacing;
    grd_start = times( &tms_grd_start);

    for (icoord[Y] = -ne[Y]; icoord[Y] <= ne[Y]; icoord[Y]++) {
        c[Y] = ((double)icoord[Y]) * spacing;

        for (icoord[X] = -ne[X]; icoord[X] <= ne[X]; icoord[X]++) {
            c[X] = ((double)icoord[X]) * spacing;

            for (j = 0;  j < num_atom_maps + 2;  j++) {
                if (gridmap[j].is_covalent == TRUE) {
                    /* Calculate the distance from the current
                     * grid point, c, to the covalent attachment point, covpos */
                    for (ii = 0;  ii < XYZ;  ii++) {
                        d[ii]  = covpos[ii] - c[ii];
                    }
                    rcov = hypotenuse( d[X], d[Y], d[Z] );
                    rcov = rcov / covhalfwidth;
                    if (rcov < APPROX_ZERO) {
                        rcov = APPROX_ZERO;
                    }
                    gridmap[j].energy = covbarrier * (1. - exp(ln_half * rcov * rcov));
                } else {
                    gridmap[j].energy = 0.; /*used to initialize to 'constant'for this gridmap*/
                } /* is not covalent */
            }
            if (floating_grid) {
                r_min = BIG;
            }
            
            /* Initialize Min Hbond variables  for each new point*/
            for (map_index = 0; map_index < num_atom_maps; map_index++){	
                hbondmin[map_index] = 999999.;
                hbondmax[map_index] = -999999.;
                hbondflag[map_index] = FALSE;
            }
                        
            /* NEW2: Find Closest Hbond */
            rmin=999999.;
            closestH=0;
            for (ia = 0;  ia < num_receptor_atoms;  ia++) {
                if ((hbond[ia]==1)||(hbond[ia]==2))  {/*DS or D1*/
                    for (i = 0;  i < XYZ;  i++) { 
                        d[i]  = coord[ia][i] - c[i]; 
                    }
                    r = hypotenuse( d[X],d[Y],d[Z] );
                    if (r < rmin) {
                        rmin = r;
                        closestH = ia; }
                } /* Hydrogen test */
            } /* ia loop */
            /* END NEW2: Find Min Hbond */

            /*
             *  Do all Receptor (protein, DNA, etc.) atoms...
             */
            for (ia = 0;  ia < num_receptor_atoms;  ia++) {
                /*
                 *  Get distance, r, from current grid point, c, to this receptor atom, coord,
                 */
                for (i = 0;  i < XYZ;  i++) {
                    d[i]  = coord[ia][i] - c[i];
                }
                r = hypotenuse( d[X], d[Y], d[Z]);
                if (r < APPROX_ZERO) {
                    r = APPROX_ZERO;
                }
                inv_r  = 1./r;
                inv_rmax = 1./max(r, 0.5);

                for (i = 0;  i < XYZ;  i++) {
                    d[i] *= inv_r;
                }
                /*make sure lookup index is in the table*/
                indx_r = min(lookup(r), MD_1);


                if (floating_grid) {
                    /* Calculate the so-called "Floating Grid"... */
                    r_min = min(r, r_min);
                }

                /* elecPE is the next-to-last last grid map, i.e. electrostatics */
                if (dddiel) {
                    /* Distance-dependent dielectric... */
                    /*gridmap[elecPE].energy += charge[ia] * inv_r * epsilon[indx_r];*/

                    /*apply the estat forcefield coefficient/weight here */
                    gridmap[elecPE].energy += charge[ia] * inv_rmax * epsilon[indx_r] * FE_coeff_estat;

                    
                } else {
                    /* Constant dielectric... */
                    /*gridmap[elecPE].energy += charge[ia] * inv_r * invdielcal;*/
                    gridmap[elecPE].energy += charge[ia] * inv_rmax * invdielcal * FE_coeff_estat;
                }

                /*
                 * If distance from grid point to atom ia is too large,
                 * or if atom is a disordered hydrogen,
                 *   add nothing to the grid-point's non-bond energy;
                 *   just continue to next atom...
                 */
                if ( r > NBCUTOFF ) {
                    continue; /* onto the next atom... */
                }
                if ((atom_type[ia] == hydrogen) && (disorder[ia] == TRUE)) {
                    continue; /* onto the next atom... */
                }

                /*** racc = rdon = 1.; ***/
                racc = 1.;
                rdon = 1.;
/* NEW2 Hramp ramps in Hbond acceptor probes */
		        Hramp = 1.;
/* END NEW2 Hramp ramps in Hbond acceptor probes */


                if (hbond[ia] == 2) {/*D1*/
                    /*
                     *  ia-th receptor atom = Hydrogen ( 4 = H )
                     *  => receptor H-bond donor, OH or NH.
                     *  calculate racc for H-bond ACCEPTOR PROBES at this grid pt.
                     *            ====     ======================
                     */
                    cos_theta = 0.;
                    /*
                     *  d[] = Unit vector from current grid pt to ia_th m/m atom.
                     *  cos_theta = d dot rvector == cos(angle) subtended.
                     */
                    for (i = 0;  i < XYZ;  i++) {
                        cos_theta -= d[i] * rvector[ia][i];
                    }

                    if (cos_theta <= 0.) {
                        /*
                         *  H->current-grid-pt vector >= 90 degrees from
                         *  N->H or O->H vector,
                         */
                        racc = 0.;
                    } else {
                        /*
                         *  racc = [cos(theta)]^2.0 for N-H
                         *  racc = [cos(theta)]^4.0 for O-H,
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
                        /* racc = pow( cos_theta, (double)rexp[ia]); */

 /* NEW2 calculate dot product of bond vector with bond vector of best hbond */
                if (ia == closestH) {
                    Hramp = 1.;
                } else {
                    cos_theta = 0.;
                    for (i = 0;  i < XYZ;  i++) {
                    cos_theta += rvector[closestH][i] * rvector[ia][i];
                    }
                   cos_theta = min(cos_theta, 1.0); 
                   cos_theta = max(cos_theta, -1.0);
                   theta = acos(cos_theta);
                   Hramp = 0.5-0.5*cos(theta * 120./90.);
                } /* ia test */
 /* END NEW2 calculate dot product of bond vector with bond vector of best hbond */
	
                    }
                    /* endif (atom_type[ia] == hydrogen) */
                    /* NEW Directional N acceptor */
                } else if (hbond[ia] == 4) {/*A1*/
            /*
            **  ia-th macromolecule atom = Nitrogen ( 4 = H )
            **  calculate rdon for H-bond Donor PROBES at this grid pt.
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
            **  X->N vector,
            */
                        rdon = 0.;
                        } else {
            /*
            **  racc = [cos(theta)]^2.0 for H->N
            */
                        rdon = cos_theta*cos_theta;
                        }
                        /* endif (atom_type[ia] == nitrogen) */
            /* end NEW Directional N acceptor */

                } else if ((hbond[ia] == 5) && (disorder[ia] == FALSE)) {/*A2*/
                    /*
                    **  ia-th receptor atom = Oxygen
                    **  => receptor H-bond acceptor, oxygen.
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
                        (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", t0);
                    } else if (t0 < -1.) {
                        t0 = -1.;
                        (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", t0);
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
                    rd2 = sq(cross[0]) + sq(cross[1]) + sq(cross[2]);
                    if (rd2 < APPROX_ZERO) {
                        if ((rd2 == 0.) && (warned == 'F')) {
                            (void) fprintf (stderr, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                            (void) fprintf (logFile, "WARNING! Attempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
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
                            (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", ti);
                        } else if (ti < -1.) {
                            ti = -1.;
                            (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", ti);
                        }
                        ti = acos(ti) - PI_halved;
                        if (ti < 0.) {
                            ti = -ti;
                        }
                        /* the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);*/
                        rdon = (0.9 + 0.1*sin(ti + ti))*cos(t0);
                    /* 0.34202 = cos (100 deg) */
                    } else if (cos_theta >= -0.34202) {
                        rdon = 562.25*pow(0.116978 - sq(cos_theta), 3.)*cos(t0);
                    }

                    /* endif atom_type == oxygen, not disordered */
                } else if ((hbond[ia] == 5) && (disorder[ia] == TRUE)) {/*A2*/

                    /* cylindrically disordered hydroxyl */
                    cos_theta = 0.;
                    for (i = 0;  i < XYZ;  i++) {
                        cos_theta -= d[i] * rvector[ia][i];
                    }
                    if (cos_theta > 1.) {
                        cos_theta = 1.;
                        (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", cos_theta);
                    } else if (cos_theta < -1.) {
                        cos_theta = -1.;
                        (void) fprintf( logFile, "WARNING! I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", cos_theta);
                    }
                    theta = acos(cos_theta);
                    racc = 0.;
                    rdon = 0.;
                    if (theta <= 1.24791 + PI_halved) {
                        /* 1.24791 rad = 180 deg minus C-O-H bond angle,
                        ** 108.5 deg */
                        rdon = pow(cos(theta - 1.24791), 4.);
                        racc = rdon;
                    }
                } /* atom_type test */

                /*
                 * For each probe atom-type,
                 * Sum pairwise interactions between each probe
                 * at this grid point (c[0:2])
                 * and the current receptor atom, ia...
                 */
                for (map_index = 0;  map_index < num_atom_maps;  map_index++) {
                    /* We do not want to change the current enrg value 
                     * for any covalent maps, make sure iscovalent is
                     * false... */                    
                    maptypeptr = gridmap[map_index].type;

                    if (gridmap[map_index].is_covalent == FALSE) {
                        if (gridmap[map_index].is_hbonder == TRUE) {
                            /*  PROBE forms H-bonds... */

                            /* rsph ramps in angular dependence for distances with negative energy */
                            rsph = energy_lookup[atom_type[ia]][indx_r][map_index]/100.;
                            rsph = max(rsph, 0.);
                            rsph = min(rsph, 1.);
                            if ((gridmap[map_index].hbond==3||gridmap[map_index].hbond==5) /*AS or A2*/
                              &&(hbond[ia]==1||hbond[ia]==2)){/*DS or D1*/
                                  /* PROBE can be an H-BOND ACCEPTOR, */
                                if (disorder[ia] == FALSE ) {
                                    gridmap[map_index].energy += energy_lookup[atom_type[ia]][indx_r][map_index] * Hramp * (racc + (1. - racc)*rsph);

                                } else {
                                    gridmap[map_index].energy += energy_lookup[hydrogen][max(0, indx_r - 110)][map_index] * Hramp * (racc + (1. - racc)*rsph);

                                }
                            } else if ((gridmap[map_index].hbond==4) /*A1*/
                              &&(hbond[ia]==1||hbond[ia]==2)) { /*DS,D1*/
	                                hbondmin[map_index] = min( hbondmin[map_index],energy_lookup[atom_type[ia]][indx_r][map_index] * (racc+(1.-racc)*rsph));
				                    hbondmax[map_index] = max( hbondmax[map_index],energy_lookup[atom_type[ia]][indx_r][map_index] * (racc+(1.-racc)*rsph));
				                    hbondflag[map_index] = TRUE;
                            } else if ((gridmap[map_index].hbond==1||gridmap[map_index].hbond==2)                                       && (hbond[ia]>2)){/*DS,D1 vs AS,A1,A2*/

                                /*  PROBE is H-BOND DONOR, */
                                temp_hbond_enrg = energy_lookup[atom_type[ia]][indx_r][map_index] * (rdon + (1. - rdon)*rsph);
                                hbondmin[map_index] = min( hbondmin[map_index], temp_hbond_enrg);
			                    hbondmax[map_index] = max( hbondmax[map_index], temp_hbond_enrg);
				                hbondflag[map_index] = TRUE;
                            } else {
                                /*  hbonder PROBE-ia cannot form a H-bond..., */
                                gridmap[map_index].energy += energy_lookup[atom_type[ia]][indx_r][map_index];
                            }

                        } else { /*end of is_hbonder*/
                            /*  PROBE does not form H-bonds..., */
                            gridmap[map_index].energy += energy_lookup[atom_type[ia]][indx_r][map_index];

                        }/* end hbonder test */
                        /* add desolvation energy  */
                        /* forcefield desolv coefficient/weight in sol_fn*/
                        gridmap[map_index].energy += gridmap[map_index].solpar_probe * vol[ia]*sol_fn[indx_r] + 
                                        (solpar[ia]+solpar_q*fabs(charge[ia]))*gridmap[map_index].vol_probe*sol_fn[indx_r];
                } /* is not covalent */
            }/* map_index */
            gridmap[dsolvPE].energy += solpar_q * vol[ia] * sol_fn[indx_r];
            }/* ia loop, over all receptor atoms... */
            for (map_index = 0; map_index < num_atom_maps; map_index++) {
			    if (hbondflag[map_index]) {
				    gridmap[map_index].energy += hbondmin[map_index]; 
				    gridmap[map_index].energy += hbondmax[map_index];
                };
            }


            /*
             * O U T P U T . . .
             *
             * Now output this grid point's energies to the maps:
             *
             */
            /*2 includes new dsolvPE*/
            for (k = 0;  k < num_atom_maps+2;  k++) {
                if (!problem_wrt) {
                    if (fabs(gridmap[k].energy) < PRECISION) {
                        fprintf_retval = fprintf(gridmap[k].map_fileptr, "0.\n");
                    } else {
                        fprintf_retval = fprintf(gridmap[k].map_fileptr, "%.3f\n", (float)gridmap[k].energy);
                    }
                    if (fprintf_retval < 0) {
                        problem_wrt = TRUE;
                    }
                }

                gridmap[k].energy_max = max(gridmap[k].energy_max, gridmap[k].energy);
                gridmap[k].energy_min = min(gridmap[k].energy_min, gridmap[k].energy);
            }
            if (floating_grid) {
                if ((!problem_wrt)&&(fprintf(floating_grid_fileptr, "%.3f\n", (float)r_min) < 0)) {
                    problem_wrt = TRUE;
                }
            }
        ctr++;
        } /* icoord[X] loop */

    } /* icoord[Y] loop */

    if (problem_wrt) {
        (void) fprintf( stderr, "WARNING! Problems writing grid maps - there may not be enough space.\n");
        (void) fprintf( logFile, "WARNING! Problems writing grid maps - there may not be enough space.\n");
    }
    grd_end = times( &tms_grd_end);
    ++nDone;
    timeRemaining = (float)(grd_end - grd_start) * idct * (float)(n1[Z] - nDone);
    (void) fprintf( logFile, " %6d   %8.3lf   %5.1lf%%   ", icoord[Z], cgridmin[Z] + c[Z], percentdone*(double)++ic);
    prHMSfixed( timeRemaining);
    (void) fprintf( logFile, "  ");
    timesys( grd_end - grd_start, &tms_grd_start, &tms_grd_end);
    (void) fflush( logFile);

} /* icoord[Z] loop */

#ifdef BOINCCOMPOUND
 boinc_fraction_done(0.9);
#endif

/*____________________________________________________________________________
 * Print a summary of extrema-values from the atomic-affinity and
 * electrostatics grid-maps,
 *____________________________________________________________________________*/

(void) fprintf(logFile, "\nGrid\tAtom\tMinimum   \tMaximum\n");
(void) fprintf(logFile, "Map \tType\tEnergy    \tEnergy \n");
(void) fprintf(logFile, "\t\t(kcal/mol)\t(kcal/mol)\n");
(void) fprintf(logFile, "___\t____\t_____________\t_____________\n");

for (i = 0;  i < num_atom_maps;  i++) {
    (void) fprintf( logFile, " %d\t %s\t  %6.2lf\t%6.2le\n", i + 1, ligand_types[i], gridmap[i].energy_min, gridmap[i].energy_max);
    /*(void) fprintf( logFile, " %d\t %c\t  %6.2lf\t%6.2le\n", i + 1, gridmap[i].atom_type, gridmap[i].energy_min, gridmap[i].energy_max);*/
}

(void) fprintf( logFile, " %d\t %c\t  %6.2lf\t%6.2le\tElectrostatic Potential\n", num_atom_maps + 1, 'e', gridmap[elecPE].energy_min, gridmap[i].energy_max);

(void) fprintf( logFile, " %d\t %c\t  %6.2lf\t%6.2le\tDesolvation Potential\n", num_atom_maps + 2, 'd', gridmap[dsolvPE].energy_min, gridmap[i+1].energy_max);
(void) fprintf( logFile, "\n\n * Note:  Every pairwise-atomic interaction was clamped at %.2f\n\n", EINTCLAMP);

/*
 * Close all files, ************************************************************
 */

for (i = 0;  i < num_atom_maps+2;  i++) {
    (void) fclose( gridmap[i].map_fileptr);
}
if (floating_grid) {
    (void) fclose(floating_grid_fileptr);
}
/* Free up the memory allocated to the gridmap objects... */
free(gridmap);

(void) fprintf( stderr, "\n%s: Successful Completion.\n", programname);
(void) fprintf( logFile, "\n%s: Successful Completion.\n", programname);

job_end = times( &tms_job_end);
timesyshms( job_end - job_start, &tms_job_start, &tms_job_end);

(void) fclose( logFile);

#ifdef BOINCCOMPOUND
 boinc_fraction_done(1.);
#endif

#ifdef BOINC	
   
    boinc_finish(0);       /* should not return */
#endif

return 0;
}
/*
 * End of main function.
 */

static int get_rec_index(const char key[]) {
    ParameterEntry * found_parm;
    found_parm = apm_find(key);
    if (found_parm != NULL)
        return found_parm->rec_index;
    return -1;
}


#ifdef BOINC
/*  Dummy graphics API entry points.
 *  This app does not do graphics, but it still must provide these callbacks.
 */

void app_graphics_render(int xs, int ys, double time_of_day) {}
void app_graphics_reread_prefs(){}
void boinc_app_mouse_move(int x, int y, bool left, bool middle, bool right ){}
void boinc_app_mouse_button(int x, int y, int which, bool is_down){}
void boinc_app_key_press(int wParam, int lParam){}
void boinc_app_key_release(int wParam, int lParam){}
#endif


/* Windows entry point WinMain() */

#ifdef _WIN32 

/*******************************************************
 * Windows:  Unix applications begin with main() while Windows applications
 * begin with WinMain, so this just makes WinMain() process the command line
 * and then invoke main()
 */

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst,
                   LPSTR Args, int WinMode)
{
    LPSTR command_line;
    char* argv[100];
    int argc;
    
    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}

#endif

/*
 * EOF
 */ 
