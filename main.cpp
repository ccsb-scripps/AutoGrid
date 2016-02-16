/*

 $Id: main.cpp,v 1.129 2016/02/16 04:56:06 mp Exp $

 AutoGrid 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

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

#include <sys/param.h>
#include <unistd.h> /* long sysconf(int name) */
#include <sys/types.h>
#include <ctype.h> /* tolower,isascii,isdigit */
#ifdef NOTNEEDED
#ifdef _WIN32
#include <Winsock2.h>
#include "util.h"
#endif
#endif

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <search.h>
#include <string.h>
#include <strings.h>  // for bzero() on Solaris
#include <time.h>
#include <stdlib.h>
#ifndef HAVE_SYSCONF
#include "mingw_sysconf.h" // for sysconf(_SC_CLK_TCK) and possibly gethostname
#endif

#include <stddef.h> 
#include <ctype.h> 
#include <cassert> 


/* the BOINC API header file */
#ifdef BOINC

#include "diagnostics.h"
#include "boinc_api.h" 
#include "filesys.h"                 // boinc_fopen(), etc... */

#endif

#include "autogrid.h"
#include "autoglobal.h"
#include "autocomm.h"
#include "constants.h"
#include "distdepdiel.h"
#include "read_parameter_library.h"
#include "timesys.h"
#include "timesyshms.h"

// M Sanner 2015 add bhtree to speed up calculation using spatial hashing
#include "bhtree.h"

extern Real idct;

// round() is a C99 function and not universally available
// Required to round %.3f consistently on different platforms
#ifdef HAVE_ROUND
#define round3dp(x) ((round((x)*1000.0L))/1000.0L)
#else
#define round3dp(x) (( floor((x)*1000.0 + 0.5)) / 1000.0)
#endif


// print_error() is used with error_level where
// error_level is defined in autogrid.h

void print_error( FILE *fileptr, int error_level, char *message) 
    // print an error or informational message to a file-pointer or
    // standard error
{
    char output_message[LINE_LEN];
    char *tag;

    switch ( error_level ) {
        case AG_ERROR:
        case FATAL_ERROR:
            tag =  "ERROR";
            break;
        case WARNING:
            tag =  "WARNING";
            break;
        default:
        case INFORMATION:
            tag = "INFORMATION";
            break;
        case SUGGESTION:
            tag =  "SUGGESTION";
            break;
    }

    (void) sprintf( output_message, "\n%s: %s:  %s\n", programname, tag, message);

    // Records all messages in the logFile.
    (void) fprintf( logFile, "%s\n", output_message);

    // Only send errors, fatal errors and warnings to standard error, stderr.
    switch ( error_level ) {
        case AG_ERROR:
        case FATAL_ERROR:
        case WARNING:
            (void) fprintf( stderr, "%s\n", output_message);
            break;
    }

    // If this is a fatal error, exit now.
    if (error_level == FATAL_ERROR) {
        (void) fprintf( logFile, "\n%s: Unsuccessful Completion.\n", programname); 
        exit( EXIT_FAILURE ); // POSIX, defined in stdlib.h (usually 1)
    }
}

/* fopen rewrite to either use BOINC api or normal system call */
FILE *ad_fopen(const char *path, const char *mode, FILE *logFile)
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
// to support use_vina_potential
static int get_map_index(const char key[]);

#ifdef _OPENMP
/* M Pique */
#include <omp.h>
 // MAXTHREADS is max number of hardware threads to use for computation.
#define MAXTHREADS 16
#else
#define omp_get_thread_num() (0)
#define omp_get_max_threads() (1)
#define MAXTHREADS 1
#endif

int main( int argc,  char **argv )

/******************************************************************************/
/*      Name: main (executable's name is "autogrid").                         */
/*  Function: Calculation of interaction energy grids for Autodock.           */
/*            Directional H_bonds from Goodford:                              */
/*            Distance dependent dielectric after Mehler and Solmajer.        */
/*            Charge-based desolvation                                        */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
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
/*   Globals: NDIEL, MAX_MAPS                                              */
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
/* use vina potential function instead of autodock4 potential function for grid calculations */
/* EXPERIMENTAL FUNCTION: NOT SUPPORTED */
const int use_vina_potential = FALSE;

/*  for associative dictionary storing parameters by autogrid 'type'  */
/*see  atom_parameter_manager.c */
static ParameterEntry thisparm;
ParameterEntry * found_parm;
char FN_parameter_library[MAX_CHARS];  // the AD4 parameters .dat file name
int parameter_library_found = 0;



/* LIGAND: 
 *  maximum is MAX_MAPS (really ATOM_MAPS) */
/*each type is now at most two characters plus '\0'*/
char ligand_types[MAX_MAPS][3];

/*array of ptrs used to parse input line*/
char * ligand_atom_types[MAX_MAPS];

/*malloc this after the number of receptor types is parsed*/
static EnergyTables et;


/* the mapObject gridmap[] holds metadata about each map. 
   AutoGrid map working storage is a subset of the map planes (in "energy") */
typedef struct mapObject {
    int    atom_type; /*corresponds to receptor numbers????*/
    int    map_index;
    int    is_covalent;
    int    is_hbonder;
    FILE   *map_fileptr;
    char   map_filename[MAX_CHARS];
    char   type[3]; /*eg HD or OA or NA or N*/
    double energy_max;
    double energy_min;
    double *(energy[MAXTHREADS]); // each possible thread is allocated n1[X]*n1[Y]*sizeof(map element [double])
    double vol_probe;
    double solpar_probe;
    /*new 6/28*/
    double Rij;
    double epsij;
    hbond_type hbond; /*hbonding character: */
    double Rij_hb;
    double epsij_hb;
    /*per receptor type parameters, ordered as in receptor_types*/
    double cA[NUM_RECEPTOR_TYPES], cB[NUM_RECEPTOR_TYPES];/*coefficients if specified in gpf*/
    double nbp_r[NUM_RECEPTOR_TYPES]; /*radius of energy-well minimum*/
    double nbp_eps[NUM_RECEPTOR_TYPES];/*depth of energy-well minimum*/
    int xA[NUM_RECEPTOR_TYPES]; /*generally 12*/
    int xB[NUM_RECEPTOR_TYPES]; /*6 for non-hbonders 10 for h-bonders*/
    int hbonder[NUM_RECEPTOR_TYPES];
} MapObject;

MapObject *gridmap = NULL; /* was statically assigned  MapObject gridmap[MAX_MAPS]; */
// convenience macro for indexing into energy array (single plane only for now)
// usage:  gridmap[mapnum].energy[plane][(mapindex(ix,iy)] = value;
#define mapindex(ix,iy) ( (ix) + ((iy)*n1[X]) )

/*variables for RECEPTOR:*/
/*each type is now at most two characters, eg 'NA\0'*/
/*NB: these are sparse arrays, some entries are not set */
char receptor_types[NUM_RECEPTOR_TYPES][3];
/* number of different receptor atom types declared on receptor_types line in GPF */
int receptor_types_gpf_ct = 0;
int has_receptor_types_in_gpf = 0;
/* number of different receptor atom types actually found in receptor PDBQT */
int receptor_types_ct = 0;
/* array of numbers of each type */
/*NB: this is a sparse int array, some entries are 0*/
int receptor_atom_type_count[NUM_RECEPTOR_TYPES];
int lc; /* receptor pdbqt file line counter */

/*array of ptrs used to parse input line*/
char * receptor_atom_types[NUM_RECEPTOR_TYPES];

// M Sanner 2015 BHTREE
#ifdef USE_BHTREE
 BHtree *bht;
 BHpoint **BHat;
#define BH_collision_dist 2.0
#define BH_cutoff_dist NBC
//#define BH_collision_dist -1.0
 int *closeAtomsIndices[MAXTHREADS];
 float *closeAtomsDistances[MAXTHREADS];
// M Sanner 2015 BHTREE END
#endif

/* AG_MAX_ATOMS */
/* changed these from "double" to "static" to reduce stack usage - MPique 2012 */
static double charge[AG_MAX_ATOMS];
static double vol[AG_MAX_ATOMS];
static double solpar[AG_MAX_ATOMS];
/*integers are simpler!*/
static int atom_type[AG_MAX_ATOMS];
static hbond_type hbond[AG_MAX_ATOMS];
static int disorder[AG_MAX_ATOMS];
static int rexp[AG_MAX_ATOMS];
static double coord[AG_MAX_ATOMS][XYZ];
static double rvector[AG_MAX_ATOMS][XYZ];
static double rvector2[AG_MAX_ATOMS][XYZ];
static int hbond_12[AG_MAX_ATOMS], nhbond_12;  // indices of hbond[ia]==1 or 2; count


/* ligand atom type number for autodock4 potential*/
static int hydrogen, nonHB_hydrogen, carbon, arom_carbon,
   oxygen, nitrogen, nonHB_nitrogen, sulphur, nonHB_sulphur,
   bromine, chlorine, fluorine, iodine;

/*canned ligand atom type number for vina_potential*/
int map_hydrogen, map_carbon, map_arom_carbon, map_oxygen, map_nitrogen; 
int map_nonHB_hydrogen, map_nonHB_nitrogen, map_sulphur, map_nonHB_sulphur;
int map_chlorine, map_bromine, map_fluorine, map_iodine;
/* XYZ */
double cext[XYZ];
double cgridmax[XYZ];
double cgridmin[XYZ];
double cmax[XYZ];
double cmin[XYZ];
double csum[XYZ];
double cmean[XYZ]={0.,0.,0.}; /* receptor mean atom coords */
double center[XYZ];
double covpos[XYZ]={0.,0.,0.}; /* Cartesian-coordinate of covalent affinity well. */
int nelements[XYZ];
static int n1[XYZ]; /* nelements[i]+1: static to detect 'not set yet' */
int ne[XYZ]; /* floor (nelements[i]/2) */

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

/* NDIEL (old name MAX_DIST) for dielectric and desolvation interactions  */
/* NEINT - for vdW and Hb interactions */
double energy_smooth[NEINT];

int nthreads=1; /* for OpenMP */
int iz; /* plane-counter */

char atom_name[6];
char token[LINE_LEN];
char warned = 'F';
static const char xyz[] = "xyz"; // used to print headings

static FILE *receptor_fileptr,
     *AVS_fld_fileptr,
     *xyz_fileptr,
     *floating_grid_fileptr;

/*for NEW3 desolvation terms*/
double solpar_q = .01097;  /*unweighted value restored 3:9:05 */
/*double solpar_q = 0.0013383; =.01097 * 0.122*/

Linear_FE_Model AD4; // set in setup_parameter_library and read_parameter_library
double q_tot = 0.0;
double diel, invdielcal=0.;//expected never used if uninitialized
double percentdone=0.0;
double PI_halved;
double q_max = -BIG,  q_min = BIG;
double ri, inv_rd, rd2, r; /* re, r2, rd, */
double inv_r, inv_rmax, racc, rdon;
double r_smooth = 0.5; //NEW ON BY DEFAULT Feb2012
//MP double Rij, epsij;
double spacing = 0.375; /* One quarter of a C-C bond length. */
double covhalfwidth = 1.0;
double covbarrier = 1000.0;
//MP double cA, cB;
const double ln_half = log(0.5);

#ifndef PACKAGE_VERSION
static char * version_num = "4.2.2";
#else
static char * version_num = PACKAGE_VERSION;
#endif

/*are these necessary??*/
double temp_vol, temp_solpar;
double Hramp=0; /*Used to cover case where vina potential is used;initialized here to quiet compiler warning*/
double factor=332.0L;  /* Used to convert between calories and SI units */

/*int num_rec_types = 0;*/

float timeRemaining = 0.;

int num_atom_maps = 0; /* number of ligand atom types, from "ligand_types" keyword */
int num_maps = 0; /* number of "map", "elecmap", or "dsolvmap" keywords handled so far */
int elecPE = -1; /* index (num_maps value) for electrostatic map, if requested */
int dsolvPE = -1; /* index (num_maps value) for desolvation map, if requested */
bool floating_grid = FALSE; // Untested for a long time - M Pique 2015
bool dddiel = FALSE, disorder_h = FALSE;
bool map_receptor_interior = FALSE; // default is to fast-skip over grid points within the receptor
 // set the value to be reported in the atomic affinity maps for interior grid points
#define INTERIOR_VALUE  9999.
int fprintf_retval = 0;
int GPF_keyword = -1;
int indcom = 0;
int infld;
int nbond;
int nDone = 0;
int problem_wrt = FALSE;
//MP int xA, xB;


int outlev = LOGFORADT;

#define INIT_NUM_GRID_PTS -1
int num_grid_points_per_map = INIT_NUM_GRID_PTS;
int i_smooth;
int ib;

int num_receptor_atoms=0;
static long clktck = 0;

Clock      job_start;
Clock      job_end;
struct tms tms_job_start;
struct tms tms_job_end;



for (int i=0; i<MAX_MAPS; i++) {
    /* initialize to "" */
    strcpy(ligand_types[i], "");
}
for (int i=0; i<NUM_RECEPTOR_TYPES; i++) {
    /* initialize to "" */
    strcpy(receptor_types[i], "");
    receptor_atom_type_count[i]=0;
}

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
    // Initialization of Boinc 
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
if (clktck == 0) {
    if ( (clktck = sysconf(_SC_CLK_TCK)) < 0) {
        (void) fprintf( stderr, "\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n");
        (void) fprintf( logFile, "\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n");
        exit(EXIT_FAILURE);
    } else {
        idct = 1. / (float)clktck;
    }
}


/*
 * Get the time at the start of the run...
 */
job_start = times( &tms_job_start);

/*
 * Parse the command line...
 */
(void) setflags( argc, argv, version_num,
#ifdef USE_BHTREE
 TRUE
#else
 FALSE
#endif
,
#ifdef _OPENMP
 TRUE
#else
 FALSE
#endif
);

 /* Initialize max and min coodinate bins */
for (int i = 0;  i < XYZ;  i++) {
        cmax[i] = -BIG;
        cmin[i] = BIG;
        csum[i] = 0.;
}

PI_halved = PI/2.;
racc = 1.; /*to quiet compiler warnings*/
rdon = 1.; /*to quiet compiler warnings*/

/*
 * Initialize int receptor_atom_type_count[] array to 0
 */
for (int i=0; i<NUM_RECEPTOR_TYPES; i++) {
    receptor_atom_type_count[i] = 0;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * Output the "AutoGrid" banner...
 */
banner( version_num);

/* report compilation options: this is mostly a duplicate of code in setflags.cpp */
(void) fprintf(logFile, "                           $Revision: 1.129 $\n");
(void) fprintf(logFile, "Compilation parameters:  NUM_RECEPTOR_TYPES=%d NEINT=%d\n",
    NUM_RECEPTOR_TYPES, NEINT);
(void) fprintf(logFile, "  AG_MAX_ATOMS=%d  MAX_MAPS=%d NDIEL=%d MAX_ATOM_TYPES=%d\n",
    AG_MAX_ATOMS, MAX_MAPS, NDIEL, MAX_ATOM_TYPES);

fprintf(logFile,"        e_vdW_Hb table has %8ld entries of size %ld\n",
  (long)(sizeof et.e_vdW_Hb/sizeof ***et.e_vdW_Hb), (long)(sizeof ***et.e_vdW_Hb));
/*
 * Print out MAX_MAPS - maximum number of maps allowed
 */
(void) fprintf(logFile, "Maximum number of maps that can be computed = %d (defined by MAX_MAPS in \"autocomm.h\").\n", MAX_MAPS);
	    fprintf(logFile, "  Non-bond cutoff for internal energy calculation (NBC): %.2f\n", NBC);
            fprintf(logFile, "  Optimize internal energy scoring (USE_8A_NBCUTOFF): ");
#ifdef USE_8A_NBCUTOFF
	    fprintf(logFile, " yes\n");
#else
	    fprintf(logFile, " no\n");
#endif
            fprintf(logFile, "  Faster search for nearby atoms (USE_BHTREE): ");
#ifdef USE_BHTREE
	    fprintf(logFile, " yes\n");
#else
	    fprintf(logFile, " no\n");
#endif
            fprintf(logFile, "  Run calculations in parallel if possible (_OPENMP): ");
#ifdef _OPENMP
	    fprintf(logFile, " yes\n");
#else
	    fprintf(logFile, " no\n");
#endif


/*
 * Print the time and date when the file was created...
 */
(void) fprintf( logFile, "This file was created at:\t\t\t");
printdate( logFile, 1);


(void) strcpy(host_name, "unknown_host");
#ifdef HAVE_GETHOSTNAME
gethostname( host_name, sizeof host_name );
#endif
(void) fprintf( logFile, "                   using:\t\t\t\"%s\"\n", host_name);


(void) fprintf( logFile, "\n\n");

//______________________________________________________________________________
//
// Read in default parameters
//
setup_parameter_library(logFile, outlev, "default Unbound_Same_As_Bound", Unbound_Same_As_Bound, &AD4);


/******************************************************************************/

/* Read in the grid parameter file...  */

while( fgets( GPF_line, LINE_LEN, GPF ) != NULL ) {

/******************************************************************************/

    GPF_keyword = gpfparser( GPF_line);

    /* This first "switch" figures out how to echo the current GPF line. */

    switch( GPF_keyword ) {

    case -1:
        (void) fprintf( logFile, "GPF> %s", GPF_line);
        print_error( logFile, WARNING, "Unrecognized keyword in grid parameter file.\n" );
        continue;  /* while fgets GPF_line... */

    case GPF_NULL:
    case GPF_COMMENT:
        (void) fprintf( logFile, "GPF> %s", GPF_line);
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
        (void) fprintf( logFile, "\nReceptor Input File :\t%s\n\nReceptor Atom Type Assignments:\n\n", receptor_filename);

        /* try to open receptor file */
        if ( (receptor_fileptr = ad_fopen(receptor_filename, "r", logFile)) == NULL ) {
            (void) sprintf( message, "can't find or open receptor PDBQT file \"%s\".\n", receptor_filename);
            print_error( logFile, FATAL_ERROR, message );
        }

        /* start to read in the lines of the receptor file */
        num_receptor_atoms = 0;
	lc = 0; /* receptor pdbqt file line counter */
        while ( (fgets(line, length, receptor_fileptr)) != NULL ) {
	    lc++;
            if (equal(line, "ATOM  ") || /* Amino Acid or DNA/RNA atoms */
                equal(line, "HETATM") || /* Non-standard heteroatoms */
                equal(line, "CHAR")) { /* Partial Atomic Charge - not a PDB record */
                /* Check that there aren't too many atoms... */
                if (num_receptor_atoms >= AG_MAX_ATOMS) {
                    (void) sprintf( message, "Too many atoms in receptor PDBQT file %s line %d;", receptor_filename, lc );
                    print_error( logFile, AG_ERROR, message );
                    (void) sprintf( message, "-- the maximum number of atoms, AG_MAX_ATOMS, allowed is %d.", AG_MAX_ATOMS );
                    print_error( logFile, AG_ERROR, message );
                    (void) sprintf( message, "Increase the value in the \"#define AG_MAX_ATOMS %d\" line", AG_MAX_ATOMS );
                    print_error( logFile, SUGGESTION, message );
                    print_error( logFile, SUGGESTION, "in the source file \"autogrid.h\", and re-compile AutoGrid." );
                    (void) fflush( logFile);
                    // FATAL_ERROR will cause AutoGrid to exit...
                    print_error( logFile, FATAL_ERROR, "Sorry, AutoGrid cannot continue.");
                } /* endif */
                /* Check that line is long enough */
                if (strlen(line) < 78) {
                    (void) sprintf( message, 
			"ATOM/HETATM line is too short in receptor PDBQT file %s line %d;",
			receptor_filename, lc );
                    print_error( logFile, AG_ERROR, message );
		}

                (void) strncpy( atom_name, &line[12], 4);
                /* atom_name is declared as an array of 6 characters,
                 * the PDB atom name is 4 characters (C indices 0, 1, 2 and 3)
                 * but let's ensure that the fifth character (C index 4)
                 * is a null character, which terminates the string. */
                atom_name[4] = '\0';

                /* Output the serial number of this atom... */
                (void) fprintf( logFile, "Atom no. %2d, \"%s\"", num_receptor_atoms + 1, atom_name);
                /* Read in this receptor atom's coordinates,partial charges, and
                 * solvation parameters in PDBQS format... */
		char field[10], field1[10], field2[10], field3[10];

		(void) strncpy(field1, &line[30], 8); field[8] = '\0';
		(void) strncpy(field2, &line[38], 8); field[8] = '\0';
		(void) strncpy(field3, &line[46], 8); field[8] = '\0';
		if ( 3 != 
                 sscanf(field1, "%lf", &coord[num_receptor_atoms][X]) +
                 sscanf(field2, "%lf", &coord[num_receptor_atoms][Y]) +
                 sscanf(field3, "%lf", &coord[num_receptor_atoms][Z]) ) {
                    (void) sprintf( message, "ATOM/HETATM line bad x,y,z in receptor PDBQT file %s line %d;", receptor_filename, lc);
                    print_error( logFile, AG_ERROR, message );
		}

                /* Output the coordinates of this atom... */
                (void) fprintf( logFile, " at (%.3lf, %.3lf, %.3lf), ",
                                coord[num_receptor_atoms][X], coord[num_receptor_atoms][Y], coord[num_receptor_atoms][Z]);

                /*1:CHANGE HERE: need to set up vol and solpar*/
		(void) strncpy(field, &line[70], 6); field[6] = '\0';
                (void) sscanf(field, "%lf", &charge[num_receptor_atoms]);
                //printf("new type is: %s\n", &line[77]);
		(void) strncpy(field, &line[77], 2); field[2] = '\0';
                (void) sscanf(field, "%s", thisparm.autogrid_type);
                found_parm = apm_find(thisparm.autogrid_type);
                if ( found_parm != NULL ) {
                    //(void) fprintf ( logFile, "DEBUG: found_parm->rec_index = %d, ->xs_radius= %f", found_parm->rec_index, found_parm->xs_radius);
                    if ( found_parm->rec_index < 0 ) {
                        strcpy( receptor_types[ receptor_types_ct ], found_parm->autogrid_type );
                        found_parm->rec_index = receptor_types_ct++;
                        //(void) fprintf ( logFile, "DEBUG: found_parm->rec_index => %d", found_parm->rec_index );
                    }
                    atom_type[num_receptor_atoms] = found_parm->rec_index;
                    solpar[num_receptor_atoms] = found_parm->solpar;
                    vol[num_receptor_atoms] = found_parm->vol;
                    hbond[num_receptor_atoms] = found_parm->hbond; /*NON=0, DS,D1, AS, A1, A2, AD */ /* N3P: added AD*/
		    // build list of "hbond" atoms to speed search M Pique Oct 2015
		    if(hbond[num_receptor_atoms]==1 || hbond[num_receptor_atoms]==2) 
			hbond_12[nhbond_12++] = num_receptor_atoms;
#ifdef DEBUG
                    printf("%d:key=%s, type=%d,solpar=%f,vol=%f\n",num_receptor_atoms,thisparm.autogrid_type, atom_type[num_receptor_atoms],solpar[num_receptor_atoms],vol[num_receptor_atoms]);
#endif
                    ++receptor_atom_type_count[found_parm->rec_index];
                } else {
		    char message[1000];
                    sprintf(message, 
		"\n\nreceptor file contains unknown type: '%s line %d'\nadd parameters for it to the parameter library first\n",
		    thisparm.autogrid_type, lc);
                    print_error(logFile, FATAL_ERROR, message);
                }
                    
                /* if from pdbqs: convert cal/molA**3 to kcal/molA**3 */
                /*solpar[num_receptor_atoms] *= 0.001;*/

                q_max = max(q_max, charge[num_receptor_atoms]);
                q_min = min(q_min, charge[num_receptor_atoms]);

                if (atom_name[0] == ' ') {
                    /* truncate the first character... */
                    atom_name[0] = atom_name[1];
                    atom_name[1] = atom_name[2];
                    atom_name[2] = atom_name[3];
                    atom_name[3] = '\0';
                } else if (isascii(atom_name[0])&& isdigit(atom_name[0]) &&
                    atom_name[1] == 'H' ) {
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
                        char temp_char    = atom_name[0];
                        atom_name[0] = atom_name[1];
                        atom_name[1] = atom_name[2];
                        atom_name[2] = atom_name[3];
                        atom_name[3] = temp_char;
                }

                /* Tell the user what you thought this atom was... */
		if(outlev>LOGRUNV)
                (void) fprintf( logFile, " was assigned atom type \"%s\" (rec_index= %d, atom_type= %d).\n", found_parm->autogrid_type, found_parm->rec_index, atom_type[num_receptor_atoms]);

                /* Count the number of each atom type */
                /*++receptor_atom_type_count[ atom_type[num_receptor_atoms] ];*/

                /* Keep track of the extents of the receptor */
                for (int i = 0;  i < XYZ;  i++) {
                    cmax[i] = max(cmax[i], coord[num_receptor_atoms][i]);
                    cmin[i] = min(cmin[i], coord[num_receptor_atoms][i]);
                    csum[i] += coord[num_receptor_atoms][i];
                }
                /* Total up the partial charges as we go... */
                q_tot += charge[num_receptor_atoms];

                /* Increment the atom counter */
                num_receptor_atoms++;

            } /* endif */
        } /* endwhile */
        /* Finished reading in the lines of the receptor file */
        (void) fclose( receptor_fileptr);
        if ( has_receptor_types_in_gpf == 1 ) {
            // Check that the number of atom types found in the receptor PDBQT
            // file match the number parsed in by the "receptor_types" command
            // in the GPF; if they do not match, exit!
            if ( receptor_types_ct != receptor_types_gpf_ct ) {
                (void) sprintf( message, "The number of atom types found in the receptor PDBQT (%d) does not match the number specified by the \"receptor_types\" command (%d) in the GPF!\n\n", receptor_types_ct, receptor_types_gpf_ct );
                print_error( logFile, FATAL_ERROR, message );
                // FATAL_ERROR will cause AutoGrid to exit...
            }
        }
        /* Update the total number of atoms in the receptor */
        (void) fprintf( logFile, "\nMaximum partial atomic charge found = %+.3lf e\n", q_max);
        (void) fprintf( logFile, "Minimum partial atomic charge found = %+.3lf e\n\n", q_min);
        (void) fflush( logFile);
        /* Check there are partial charges... */
        if (q_max == 0. && q_min == 0.) {
            (void) sprintf( message, "No partial atomic charges were found in the receptor PDBQT file %s!\n\n", receptor_filename );
            print_error( logFile, FATAL_ERROR, message );
            // FATAL_ERROR will cause AutoGrid to exit...
        } /* if  there are no charges EXIT*/

        for (int ia = 0;  ia < num_receptor_atoms;  ia++) {
            rexp[ia] = 0;
        }

        (void) fprintf( logFile, "Atom\tAtom\tNumber of this Type\n");
        (void) fprintf( logFile, "Type\t ID \t in Receptor\n");
        (void) fprintf( logFile, "____\t____\t___________________\n");
        /*2. CHANGE HERE: need to count number of each receptor_type*/
        for (int ia = 0;  ia < receptor_types_ct;  ia++) {
            //i = 0;
            if(receptor_atom_type_count[ia]!=0){
                (void) fprintf( logFile, " %d\t %s\t\t%6d\n", (ia), receptor_types[ia], receptor_atom_type_count[ia]);
                //i++;
            };
        }
        (void) fprintf( logFile, "\nTotal number of atoms :\t\t%d atoms \n", num_receptor_atoms);
        (void) fprintf( logFile, "Total charge :\t\t\t%.2lf e\n", q_tot);
        (void) fprintf( logFile, "\n\nReceptor coordinates fit within the following volume:\n\n");
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
            print_error( logFile, FATAL_ERROR, "Grid data file needs the extension \".fld\" for AVS input\n\n" );
        } else {
            infld = strindex( AVS_fld_filename, "fld");
            (void) strcpy(xyz_filename, AVS_fld_filename);
            xyz_filename[infld] = 'x';
            xyz_filename[infld + 1] = 'y';
            xyz_filename[infld + 2] = 'z';
        }
        if ( (AVS_fld_fileptr = ad_fopen(AVS_fld_filename, "w", logFile)) == NULL ) {
            (void) sprintf( message, "can't create grid dimensions data file %s\n", AVS_fld_filename);
            print_error( logFile, FATAL_ERROR, message );
        } else {
            (void) fprintf( logFile, "\nCreating (AVS-readable) grid maps file : %s\n", AVS_fld_filename);
        }
        if ( (xyz_fileptr = ad_fopen(xyz_filename, "w", logFile)) == NULL ) {
            (void) sprintf( message, "can't create grid extrema data file %s\n", xyz_filename);
            print_error( logFile, FATAL_ERROR, message );
        } else {
            (void) fprintf( logFile, "\nCreating (AVS-readable) grid-coordinates extrema file : %s\n\n", xyz_filename);
        }
        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_NPTS:
        (void) sscanf( GPF_line, "%*s %d %d %d", &nelements[X], &nelements[Y], &nelements[Z]);
        for (int i = 0;  i < XYZ;  i++) {
            nelements[i] = check_size(nelements[i], xyz[i]); // will be even: 0,2,4, etc: 
            ne[i] = nelements[i] / 2; // how many non-zero indices in each of positive and negative
            n1[i] = nelements[i] + 1; // true size of array - odd number
	    // example: user specifies 'npts 8'. check_size will return nelements=8.
	    //  ne will be 4
	    //  n1 will be 9
	    //  the first element in the map will be -3 units below the origin.
	    // example: user specifies 'npts 1'. check_size will return nelements=0.
	    //  ne will be 0
	    //  n1 will be 1
	    //  the first element in the map will be 0 units below the origin.
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
        break;

/******************************************************************************/

    case GPF_GRIDCENTER:
       if ( AVS_fld_fileptr == NULL || xyz_fileptr == NULL ) {
            print_error( logFile, FATAL_ERROR, 
            "You need to set the \"gridfld\" file before setting the grid center\".\n" );
        }

        (void) sscanf( GPF_line, "%*s %s", token);
        if (equal( token, "auto")) {
            for (int i = 0;  i < XYZ;  i++) {
                center[i] = cmean[i];
            }
            (void) fprintf( logFile, "Grid maps will be centered on the center of mass.\n");
            (void) fprintf( logFile, "Coordinates of center of mass : (%.3lf, %.3lf, %.3lf)\n", center[X], center[Y], center[Z]);
        } else {
            (void) sscanf( GPF_line, "%*s %lf %lf %lf", &center[X], &center[Y], &center[Z]);
            (void) fprintf( logFile, "\nGrid maps will be centered on user-defined coordinates:\n\n\t\t(%.3lf, %.3lf, %.3lf)\n",  center[X], center[Y], center[Z]);
        }
        /* centering stuff... */
        for (int ia = 0;  ia < num_receptor_atoms;  ia++) {
            for (int i = 0;  i < XYZ;  i++) {
                coord[ia][i] -= center[i];        /* transform to center of gridmaps */
            }
        }
        for (int i = 0;  i < XYZ;  i++) {
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
        for (int i = 0;  i < XYZ;  i++) {
            (void) fprintf( logFile, "Grid map %c-dimension :\t\t%.1lf Angstroms\n", xyz[i], 2.*cext[i]);
        }
        (void) fprintf( logFile, "\nMaximum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n", cgridmax[X], cgridmax[Y], cgridmax[Z]);
        (void) fprintf( logFile, "Minimum coordinates :\t\t(%.3lf, %.3lf, %.3lf)\n\n", cgridmin[X], cgridmin[Y], cgridmin[Z]);
        for (int i = 0;  i < XYZ;  i++) {
            (void) fprintf(xyz_fileptr, "%.3lf %.3lf\n", cgridmin[i], cgridmax[i]);
        }
        (void) fclose(xyz_fileptr);
	xyz_fileptr = NULL;
        (void) fflush( logFile);
        break;

/******************************************************************************/

    case GPF_LIGAND_TYPES:
        // Read in the list of atom types in the ligand.
        // GPF_line e.g.: "ligand_types N O A C HH NH"
       if ( n1[X] == 0 ) {
            print_error( logFile, FATAL_ERROR, 
            "You need to set the \"npts\" before setting the ligand types\".\n" );
        }
        num_atom_maps = parsetypes(GPF_line, ligand_atom_types, MAX_ATOM_TYPES);
        for (int i=0; i<num_atom_maps; i++) {
            strcpy(ligand_types[i], ligand_atom_types[i]);
#ifdef DEBUG
            (void) fprintf(logFile, "%d %s ->%s\n",i, ligand_atom_types[i], ligand_types[i]);
#endif
        }
        for (int i=0; i<num_atom_maps; i++){
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
		char message[1000];
                (void) sprintf( message, "unknown ligand atom type %s\nadd parameters for it to the parameter library first!\n", ligand_atom_types[i]);
                 print_error(logFile, FATAL_ERROR, message);
             }
        };

        /* num_maps is the number of maps to be created:
         * the number of ligand atom types, 
	 * plus (optionally) 1 for the electrostatic map 
	 * plus (optionally) 1 for the desolvation map.
	 * plus (optionally) 1 for the floating grid map.
         * AutoDock can only read in MAX_MAPS maps, which must include
         * the ligand atom maps and electrostatic, desolvation, and floating grid maps */
      
        /* Check to see if there is enough memory to store these map objects */
        gridmap = (MapObject *)calloc(sizeof(MapObject), num_atom_maps+3);

        if ( gridmap == NULL ) {
            (void) sprintf( message, "Too many ligand atom types; there is not enough memory to create these maps.  Try using fewer atom types than %d.\n", num_atom_maps);
            print_error( logFile, FATAL_ERROR, message);
        }


	nthreads = min(MAXTHREADS, min(n1[Z], omp_get_max_threads()));
	fprintf(logFile, "%d CPU thread(s) will be used for calculation\n", nthreads);

        // Initialize the gridmap MapObject: the floating grid does not have one
	// (The "energy" is working storage for each thread - an XY plane)
        for (int i=0; i<num_atom_maps+2; i++) {
            gridmap[i].atom_type = 0; /*corresponds to receptor numbers????*/
            gridmap[i].map_index = 0;
            gridmap[i].is_covalent = 0;
            gridmap[i].is_hbonder = 0;
            gridmap[i].map_fileptr = (FILE *)NULL;
            strcpy(gridmap[i].map_filename, "");
            strcpy(gridmap[i].type,""); /*eg HD or OA or NA or N*/
            gridmap[i].energy_max = 0.0L;
            gridmap[i].energy_min = 0.0L;
	    for (int p=0; p<nthreads; p++) {
               gridmap[i].energy[p] = (double *) malloc( n1[X]*n1[Y]*sizeof(double));
               if(0==gridmap[i].energy[p]) {	
		fprintf(logFile, "MALLOC_ERROR map i=%d plane p=%d *energy=%ld\n", i, p, (long)gridmap[i].energy[p]);
		print_error(logFile, FATAL_ERROR, "unable to allocate plane energy storage");
	        }
	    }
            gridmap[i].vol_probe = 0.0L;
            gridmap[i].solpar_probe = 0.0L;
            gridmap[i].Rij = 0.0L;
            gridmap[i].epsij = 0.0L;
            gridmap[i].hbond = NON; /*hbonding character: */
            gridmap[i].Rij_hb = 0.0L;
            gridmap[i].epsij_hb = 0.0L;
            /*per gridmap[i].receptor type parameters, ordered as in receptor_types*/
            for (int j=0; j<NUM_RECEPTOR_TYPES; j++) {
                gridmap[i].cA[j] =0; /*default is to automatically calculate from r,eps*/
                gridmap[i].cB[j] =0; /*ditto*/
                gridmap[i].nbp_r[j] = 0.0L; /*radius of energy-well minimum*/
                gridmap[i].nbp_eps[j] = 0.0L;/*depth of energy-well minimum*/
                gridmap[i].xA[j] =0; /*generally 12*/
                gridmap[i].xB[j] =0; /*6 for non-hbonders 10 for h-bonders*/
                gridmap[i].hbonder[j] =0;
            } // j
        } // i


        if (num_grid_points_per_map == INIT_NUM_GRID_PTS) {
            print_error( logFile, FATAL_ERROR, "You need to set the number of grid points using \"npts\" before setting the ligand atom types, using \"ligand_types\".\n" );
        }

        for (int i = 0;  i < num_atom_maps;  i++) {
            gridmap[i].is_covalent = FALSE;
            gridmap[i].is_hbonder = FALSE;
            gridmap[i].map_index = i;
            strcpy(gridmap[i].type, ligand_types[i]); /*eg HD or OA or NA or N*/
            found_parm = apm_find(ligand_types[i]);
            if (strcmp(ligand_types[i],"Z")==0){
                fprintf(logFile, "Found covalent map atomtype\n");
                gridmap[i].is_covalent = TRUE;}
            gridmap[i].atom_type = found_parm->map_index;
            gridmap[i].solpar_probe = found_parm->solpar;
            gridmap[i].vol_probe = found_parm->vol;
            gridmap[i].Rij = found_parm->Rij;
            gridmap[i].epsij = found_parm->epsij;
            gridmap[i].hbond = found_parm->hbond;
            gridmap[i].Rij_hb = found_parm->Rij_hb;
            gridmap[i].epsij_hb = found_parm->epsij_hb;
            if (gridmap[i].hbond>0){       //enum: NON,DS,D1,AS,A1,A2,AD /* N3P added AD type*/
                gridmap[i].is_hbonder=TRUE;}

#ifdef DEBUG
            (void) fprintf(logFile, " setting ij parms for map %d \n",i);
            (void) fprintf(logFile, "for gridmap[%d], type->%s,Rij->%6.4f, epsij->%6.4f, hbond->%d\n",i,found_parm->autogrid_type, gridmap[i].Rij, gridmap[i].epsij,gridmap[i].hbond);
#endif
            for (int j=0; j<receptor_types_ct; j++){
                found_parm = apm_find(receptor_types[j]);
                gridmap[i].nbp_r[j] = (gridmap[i].Rij + found_parm->Rij)/2.;
                gridmap[i].nbp_eps[j] = sqrt(gridmap[i].epsij * found_parm->epsij);
                /*apply the vdW forcefield parameter/weight here */
                // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_vdW;
                gridmap[i].xA[j] = 12;
                /*setup hbond dependent stuff*/
                gridmap[i].xB[j] = 6;
                gridmap[i].hbonder[j] = 0;
                if ( ((int)(gridmap[i].hbond)>2 || (int)(gridmap[i].hbond==6) )&&
                        ((int)found_parm->hbond==1||(int)found_parm->hbond==2||(int)found_parm->hbond==6)){ /*AS,A1,A2, AD map vs DS,D1,AD probe N3P modified */
                    gridmap[i].xB[j] = 10;
                    gridmap[i].hbonder[j] = 1;
                    gridmap[i].is_hbonder = TRUE;
                    /*Rij and epsij for this hb interaction in
                     * parm_data.dat file as Rii and epsii for heavy atom
                     * hb factors*/
                    gridmap[i].nbp_r[j] = gridmap[i].Rij_hb;
                    gridmap[i].nbp_eps[j] = gridmap[i].epsij_hb;

                    /*apply the hbond forcefield parameter/weight here */
                    // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_hbond;
#ifdef DEBUG
                    (void) fprintf(logFile, "set %d-%d hb eps to %6.4f*%6.4f=%6.4f\n",i,j,gridmap[i].epsij_hb,found_parm->epsij_hb, gridmap[i].nbp_eps[j]);
#endif
                } else if ( ( (int)gridmap[i].hbond==1||(int)gridmap[i].hbond==2||(int)gridmap[i].hbond==6) &&
                            (((int)found_parm->hbond>2))) { /*DS,D1,AS map vs AS,A1,A2,AS probe N3P: modified*/
                    gridmap[i].xB[j] = 10;
                    gridmap[i].hbonder[j] = 1;
                    gridmap[i].is_hbonder = TRUE;
                    /*Rij and epsij for this hb interaction in
                     * parm_data.dat file as Rii and epsii for heavy atom
                     * hb factors*/
                    gridmap[i].nbp_r[j] = found_parm->Rij_hb;
                    gridmap[i].nbp_eps[j] = found_parm->epsij_hb;

                    /*apply the hbond forcefield parameter/weight here */
                    // This was removed because "setup_p_l" does this for us... gridmap[i].nbp_eps[j] *= FE_coeff_hbond;
#ifdef DEBUG
                    (void) fprintf(logFile, "2: set %d-%d hb eps to %6.4f*%6.4f=%6.4f\n",i,j,gridmap[i].epsij_hb,found_parm->epsij_hb, gridmap[i].nbp_eps[j]);
#endif
                } 
#ifdef DEBUG
                (void) fprintf(logFile, "vs receptor_type[%d]:type->%s, hbond->%d ",j,found_parm->autogrid_type, (int)found_parm->hbond);
                (void) fprintf(logFile, "nbp_r->%6.4f, nbp_eps->%6.4f,xB=%d,hbonder=%d\n",gridmap[i].nbp_r[j], gridmap[i].nbp_eps[j],gridmap[i].xB[j], gridmap[i].hbonder[j]);
#endif
            } /*initialize energy parms for each possible receptor type*/
        } /*for each map*/
        (void) fprintf( logFile, "\nAtom type names for ligand atom types 1-%d used for ligand-atom affinity grid maps:\n\n", num_atom_maps);
        for (int i = 0;  i < num_atom_maps;  i++) {
            (void) fprintf( logFile, "\t\t\tAtom type number %d corresponds to atom type name \"%s\".\n", gridmap[i].map_index, gridmap[i].type);
            if (gridmap[i].is_covalent == TRUE) {
              (void) fprintf( logFile, "\nAtom type number %d will be used to calculate a covalent affinity grid map\n\n", i + 1);
            }
        }
        // at this point set up map_hydrogen, map_carbon, map_oxygen and map_nitrogen etc for vina potential
        map_hydrogen = get_map_index("HD");
        map_nonHB_hydrogen = get_map_index("H");
        map_carbon = get_map_index("C");
        map_arom_carbon = get_map_index("A");
        map_oxygen = get_map_index("OA");
        map_nitrogen = get_map_index("NA");
        map_nonHB_nitrogen = get_map_index("N");
        map_sulphur = get_map_index("SA");
        map_nonHB_sulphur = get_map_index("S");
        map_fluorine = get_map_index("F");
        map_bromine = get_map_index("Br");
        map_chlorine = get_map_index("Cl");
        map_iodine = get_map_index("I");
        (void) fprintf( logFile, "\n\n");
        break;

/******************************************************************************/

    case GPF_RECEPTOR_TYPES:
        // Read in the list of atom types in the receptor.
        // GPF_line e.g.: "receptor_types N O A C HH NH"
        //
        // NOTE:  This line is not guaranteed to match the actual
        //        atom types present in the receptor PDBQT file
        //        specified by the "receptor" command.
        receptor_types_ct = parsetypes(GPF_line, receptor_atom_types,  MAX_ATOM_TYPES);
        receptor_types_gpf_ct = receptor_types_ct;
        has_receptor_types_in_gpf = 1;
#ifdef DEBUG
        printf("receptor_types_gpf_ct=%d\n",receptor_types_gpf_ct);
        printf("receptor_types_ct=%d\n",receptor_types_ct);
#endif
        for (int i=0; i<receptor_types_ct; i++){
            strcpy(receptor_types[i], receptor_atom_types[i]);
#ifdef DEBUG
            printf("%d %s  ->%s\n",i, receptor_atom_types[i],  receptor_types[i]);
#endif
        }
        for (int i=0; i<receptor_types_ct; i++) {
            found_parm = apm_find(receptor_atom_types[i]);
            if (found_parm != NULL){
                found_parm->rec_index = i;
            } else {
                (void) sprintf( message, "Unknown receptor type: \"%s\"\n -- Add parameters for it to the parameter library first!\n", receptor_atom_types[i]);
                print_error( logFile, FATAL_ERROR, message );
            }
        }
        // at this point set up hydrogen, carbon, oxygen and nitrogen
        hydrogen = get_rec_index("HD");
        nonHB_hydrogen = get_rec_index("H");
        carbon = get_rec_index("C");
        arom_carbon = get_rec_index("A");
        oxygen = get_rec_index("OA");
        nitrogen = get_rec_index("NA");
        nonHB_nitrogen = get_rec_index("N");
        sulphur = get_rec_index("SA");
        nonHB_sulphur = get_rec_index("S");
        bromine = get_rec_index("Br");
        chlorine = get_rec_index("Cl");
        fluorine = get_rec_index("Fl");
        iodine = get_rec_index("I");
#ifdef DEBUG
        printf("assigned receptor types:arom_carbon->%d, hydrogen->%d,nonHB_hydrogen->%d, carbon->%d, oxygen->%d, nitrogen->%d\n, nonHB_nitrogen->%d, sulphur->%d, nonHB_sulphur->%d\n",arom_carbon,hydrogen, nonHB_hydrogen, carbon,oxygen, nitrogen, nonHB_nitrogen, sulphur, nonHB_sulphur);
#endif
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
            int mapi = found_parm->map_index;
            if (mapi>=0){
                /*DON'T!!!*/
                /*convert cal/molA^3 to kcal/molA^3 */
                /*gridmap[mapi].solpar_probe = temp_solpar * 0.001;*/
                gridmap[mapi].solpar_probe = temp_solpar ;
                (void) fprintf( logFile, "\nProbe %s solvation parameters: \n\n\tatomic fragmental volume: %.2f A^3\n\tatomic solvation parameter: %.4f cal/mol A^3\n\n", found_parm->autogrid_type, found_parm->vol,found_parm->solpar);
            }
        } else {
            (void) fprintf( logFile, "%s key not found\n", thisparm.autogrid_type);
        };
        break; /* end solvation parameter */

/******************************************************************************/

    /*case GPF_CONSTANT:*/
        /*break;*/

/******************************************************************************/
/******************************************************************************/

    case GPF_USE_VINA_POTENTIAL:
        //use_vina_potential = TRUE;
        (void) fprintf( logFile, "\n Using Vina potential for calculation.\n\n");
        print_error( logFile, FATAL_ERROR, "Vina potential not implemented" );
        break;

/******************************************************************************/


    case GPF_MAP: 
        /*  */
        /* The variable "num_maps" is the 0-based index of the ligand atom type
         * we are calculating a map for. 
         * If the "types" line was C N O S H, there would be 5 ligand atom maps to calculate,
         * num_maps will increment
         * each time there is a "map" keyword in the GPF has been processed.  The value of
         * num_maps should therefore go from 1 to 5 after each "map" keyword.  
         * In this example, num_atom_maps would be 5, 
         * so if num_maps is > 4, there is something wrong in the number of
         * "map" keywords. */
        if (num_maps+1 > num_atom_maps) {
             (void) sprintf(message, "Too many \"map\" keywords (%d);  the \"ligand_types\" command declares only %d atom types.\nRemove a \"map\" keyword from the GPF.\n", num_maps + 1, num_atom_maps);
            print_error( logFile, FATAL_ERROR, message );
        }
        /* Read in the filename for this grid map */ /* GPF_MAP */
        (void) sscanf( GPF_line, "%*s %s", gridmap[num_maps].map_filename);
        if ( (gridmap[num_maps].map_fileptr = ad_fopen( gridmap[num_maps].map_filename, "w", logFile)) == NULL ) {
            (void) sprintf( message, "Cannot open grid map \"%s\" for writing.", gridmap[num_maps].map_filename);
            print_error( logFile, FATAL_ERROR, message );
        }
        (void) fprintf( logFile, "\nOutput Grid Map %d:   %s\n\n", (num_maps + 1), gridmap[num_maps].map_filename);
	num_maps++;

        break;

/******************************************************************************/
    case GPF_ELECMAP:
	if(elecPE>=0) {
            print_error( logFile, FATAL_ERROR, "Duplicate \"elecmap\" request");
		}
	if(num_maps<num_atom_maps) {
            print_error( logFile, FATAL_ERROR, "\"elecmap\" request appeared before all atom map requests had appeared.");
		}
	elecPE = num_maps;
        (void) sscanf( GPF_line, "%*s %s", gridmap[elecPE].map_filename);
        if ( (gridmap[elecPE].map_fileptr = ad_fopen( gridmap[elecPE].map_filename, "w", logFile)) == NULL){
            (void) sprintf( message, "can't open grid map \"%s\" for writing.\n", gridmap[elecPE].map_filename);
            print_error( logFile, FATAL_ERROR, message );
        }
        (void) fprintf( logFile, "\nOutput Electrostatic Potential Energy Grid Map: %s\n\n", gridmap[elecPE].map_filename);
	num_maps++;
        break;

/******************************************************************************/
    case GPF_DSOLVMAP:
	if(dsolvPE>=0) {
            print_error( logFile, FATAL_ERROR, "Duplicate \"dsolvmap\" request");
		}
	if(num_maps<num_atom_maps) {
            print_error( logFile, FATAL_ERROR, "\"desolvmap\" request appeared before all atom map requests had appeared.");
		}
	dsolvPE = num_maps;
        (void) sscanf( GPF_line, "%*s %s", gridmap[dsolvPE].map_filename);
        if ( (gridmap[dsolvPE].map_fileptr = ad_fopen( gridmap[dsolvPE].map_filename, "w", logFile)) == NULL){
            (void) sprintf( message, "can't open grid map \"%s\" for writing.\n", gridmap[dsolvPE].map_filename);
            print_error( logFile, FATAL_ERROR, message );
        }
        (void) fprintf( logFile, "\nOutput Desolvation Free Energy Grid Map: %s\n\n", gridmap[dsolvPE].map_filename);
	num_maps++;
        break;
/******************************************************************************/

    case GPF_COVALENTMAP:
        (void) sscanf( GPF_line, "%*s %lf %lf %lf %lf %lf", &covhalfwidth, &covbarrier, &(covpos[X]), &(covpos[Y]), &(covpos[Z]));
        (void) fprintf( logFile, "\ncovalentmap <half-width in Angstroms> <barrier> <x> <y> <z>\n");
        (void) fprintf( logFile, "\nCovalent well's half-width in Angstroms:         %8.3f\n", covhalfwidth);
        (void) fprintf( logFile, "\nCovalent barrier energy in kcal/mol:             %8.3f\n", covbarrier);
        (void) fprintf( logFile, "\nCovalent attachment point will be positioned at: (%8.3f, %8.3f, %8.3f)\n\n", covpos[X], covpos[Y], covpos[Z]);
        for (int i = 0;  i < XYZ;  i++) {
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
        (void) fprintf( logFile, "\nPotentials will be smoothed by: %.3lf Angstrom\n\n", r_smooth);
        break;

/******************************************************************************/

    case GPF_OUTLEV:
        (void) sscanf( GPF_line, "%*s %d", &outlev);
        (void) fprintf( logFile, "\nOutput level: %d\n\n", outlev);
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
            dddiel = TRUE;
            /* calculate ddd of Mehler & Solmajer */
            (void) fprintf( logFile, "\nUsing *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.\n\n");
            et.epsilon_fn[0] = 1.0;
            for (int indx_r = 1;  indx_r < NDIEL;  indx_r++) {
                et.epsilon_fn[indx_r] = calc_ddd_Mehler_Solmajer( angstrom(indx_r), APPROX_ZERO );
            }
            (void) fprintf( logFile, "  d   Dielectric\n ___  __________\n");
            for (int i = 0;  i <= min(500,NDIEL);  i += 10) {
                ri = angstrom(i);
                (void) fprintf( logFile, "%4.1lf%9.2lf\n", ri, et.epsilon_fn[i]);
            }
            (void) fprintf( logFile, "\n");
            /* convert epsilon to factor / epsilon */
            for (int i = 0;  i < NDIEL;  i++) {
                et.r_epsilon_fn[i] = factor/et.epsilon_fn[i];
            }

        } else {
            /* positive or zero... */
            dddiel = FALSE;
            if (diel <= APPROX_ZERO) {
                diel = 40.;
            }
            (void) fprintf( logFile, "Using a *constant* dielectric of:  %.2f\n", diel);
            invdielcal = factor / diel;
        }
        break;

/******************************************************************************/
    case GPF_MAP_RECEPTOR_INTERIOR:
	map_receptor_interior = TRUE;
        (void) fprintf( logFile, "Mapping receptor interior instead of reporting collision.");
	break;

/******************************************************************************/

    case GPF_FMAP:
	// Caution: not tested! M Pique 2015
        (void) sscanf( GPF_line, "%*s %s", floating_grid_filename);
        if ( (floating_grid_fileptr = ad_fopen( floating_grid_filename, "w", logFile)) == NULL) {
            (void) sprintf( message, "can't open grid map \"%s\" for writing.\n", floating_grid_filename);
            print_error( logFile, FATAL_ERROR, message );
        }
        (void) fprintf( logFile, "\nFloating Grid file name = %s\n", floating_grid_filename);
	/* do NOT increment num_maps because (for unknown reasons) this grid isn't counted - MPique */
        floating_grid = TRUE;
        break;

/******************************************************************************/

    case GPF_PARAM_FILE:
        /* open and read the AD4 parameters .dat file */

        parameter_library_found = sscanf( GPF_line, "%*s %s ", FN_parameter_library);

        read_parameter_library(logFile, outlev, FN_parameter_library, &AD4);

        break;


/******************************************************************************/

    case GPF_NBP_COEFFS:
    case GPF_NBP_R_EPS:
        /* 
        ** nbp_r_eps
        ** override energy parameters:
        ** Lennard-Jones and Hydrogen Bond Potentials,
        ** GPF_NBP_REQM_EPS: Using epsilon and r-equilibrium values...
        ** for the interaction of the specified types 
        */

	{  /* block for local storage allocation name space */
        Real epsij;
        Real Rij;
        Real cA, cB;
        char param[2][LINE_LEN];
        int xA, xB, nfields;

        nfields = sscanf( GPF_line, "%*s " FDFMT2 " %d %d %s %s", &Rij, &epsij, &xA, &xB, param[0], param[1] );
        if(nfields!=6) {
         (void) sprintf( message,  "syntax error, not 6 values in NBP_R_EPS line");
         print_error(logFile, FATAL_ERROR, message);
        }
        /* check values? 
        if ((Rij < RIJ_MIN) || (Rij > RIJ_MAX)) {
            (void) fprintf( logFile,
	    "WARNING: pairwise distance, Rij, %.2f, is not a very reasonable value for the equilibrium separation of two atoms! (%.2f Angstroms <= Rij <= %.2f Angstroms)\n\n", Rij, RIJ_MIN, RIJ_MAX);        
	}
        */

        cA = Rij;
        cB = epsij;
        if ( GPF_keyword == GPF_NBP_R_EPS ) {
           // Calculate the coefficients from Rij and epsij                      
           /* Defend against division by zero... */
           if (xA != xB) { 
    	       double tmpconst = epsij / (Real)(xA - xB);
               cA = tmpconst * pow( (double)Rij, (double)xA ) * (Real)xB;
               cB = tmpconst * pow( (double)Rij, (double)xB ) * (Real)xA;
           } else {           
               (void) sprintf( message, "exponents xA and xB cannot be equal.\n");
               print_error( logFile, FATAL_ERROR, message );
          } 
	}
 
        for (int i=0;i<2;i++) {
	    /* try both orderings of "ligand,receptor" and "receptor,ligand": not error if not found  */
            int ligtype, rectype;
            ligtype = get_map_index(param[i%2]);
            rectype = get_rec_index(param[(i+1)%2]);
            if (ligtype>=0 && rectype>=0){
                pr(logFile, "\n nbp_r_eps or nbp_coeffs: map_index(%s)= %d  rec_index(%s)= %d\n",param[i%2],ligtype, param[(i+1)%2],rectype);
                gridmap[ligtype].cA[rectype] = cA;
                gridmap[ligtype].cB[rectype] = cB;
                gridmap[ligtype].nbp_r[rectype] = Rij;
                gridmap[ligtype].nbp_eps[rectype] = epsij;
                gridmap[ligtype].xA[rectype] = xA;
                gridmap[ligtype].xB[rectype] = xB;
            }
        }
        pr(logFile, "\nOverriding non-bonded interaction energies for docking calculation;\n");
        break;
	}  /* end block for local storage allocation name space */



/******************************************************************************/

    default:
        break;

/******************************************************************************/

    } /* second switch */

} /* while: finished reading gpf */

#ifdef USE_BHTREE
// M Sanner 2015 BHTREE
// build BHTREE for receptor atoms	
 BHat = (BHpoint **)malloc(num_receptor_atoms*sizeof(BHpoint *));

 for (int ia=0;ia<num_receptor_atoms;ia++) {
   BHat[ia] = (BHpoint *)malloc(sizeof(BHpoint));
   BHat[ia]->x[0]=coord[ia][X];
   BHat[ia]->x[1]=coord[ia][Y];
   BHat[ia]->x[2]=coord[ia][Z];
   BHat[ia]->r=2.0;
   BHat[ia]->at=ia;
 }
 bht = generateBHtree(BHat, num_receptor_atoms, 10);
// M Sanner 2015 BHTREE END
	
for(int p=0;p<nthreads;p++) {
  closeAtomsIndices[p] = (int *)malloc(num_receptor_atoms*sizeof(int));
  closeAtomsDistances[p] = (float *)malloc(num_receptor_atoms*sizeof(float));
}
#endif

/* Map files checkpoint  (number of maps, desolv and elec maps )    SF  */

/* Number of maps defined for atom types*/
//DEBUG fprintf(stderr,"\n  num_atom_maps %d,  num_maps-2 %d\n", num_atom_maps, num_maps-2);fflush(stderr);
if ( num_atom_maps  > (num_maps - (elecPE>=0?1:0) - (dsolvPE>=0?1:0)) ) {   
             (void) fprintf( logFile, "Too few \"map\" keywords ;  the \"ligand_types\" command declares %d atom types.\nAdd a \"map\" keyword from the GPF.\n", num_atom_maps );
             (void) sprintf( message, "Too few \"map\" keywords found for the number of ligand atom types.\n" );
            print_error( logFile, FATAL_ERROR, message );
	} 

/* Electrostatic map (optional, but warn if not requested) */
if (( not use_vina_potential) && (elecPE<0 || (strlen( gridmap[elecPE].map_filename ) == 0 ))) {  
             (void) fprintf( logFile, "The electrostatic map file is not defined in the GPF.\n" );
             (void) sprintf( message, "No electrostatic map file requested.\n" );
            print_error( logFile, WARNING, message );
	}
/* Desolvation map (optional, but warn if not requested) */
if (( not use_vina_potential) && (dsolvPE<0 || (strlen( gridmap[dsolvPE].map_filename ) == 0 ))) {  
             (void) fprintf( logFile, "The desolvation map file is not defined in the GPF.\n" );
             (void) sprintf( message, "No desolvation map file requested.\n" );
            print_error( logFile, WARNING, message );
	}

/* End of map files checkpoint SF */


(void) fprintf( logFile, "\n>>> Closing the grid parameter file (GPF)... <<<\n\n");
(void) fprintf( logFile, UnderLine);
(void) fclose( GPF );

if ( ! floating_grid ) {
    // (void) fprintf( logFile, "\n\nNo Floating Grid was requested.\n");
}

if ( AVS_fld_fileptr == NULL ) {
             (void) fprintf( logFile, "Missing \"gridfld\" keyword;\nAdd a \"gridfld <file.fld>\" keyword to the GPF after the \"npts\" keyword.\n");
             (void) sprintf( message, "Missing \"gridfld\" keyword.");
            print_error( logFile, FATAL_ERROR, message );
	} 

(void) fprintf( AVS_fld_fileptr, "# AVS field file\n#\n");
(void) fprintf( AVS_fld_fileptr, "# AutoDock Atomic Affinity and Electrostatic Grids\n#\n");
(void) fprintf( AVS_fld_fileptr, "# Created by %s %s.\n#\n", programname, version_num);
(void) fprintf( AVS_fld_fileptr, "#SPACING %.3f\n", (float) spacing);
(void) fprintf( AVS_fld_fileptr, "#NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
(void) fprintf( AVS_fld_fileptr, "#CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
(void) fprintf( AVS_fld_fileptr, "#MACROMOLECULE %s\n", receptor_filename);
(void) fprintf( AVS_fld_fileptr, "#GRID_PARAMETER_FILE %s\n#\n", grid_param_fn );
(void) fprintf( AVS_fld_fileptr, "ndim=3\t\t\t# number of dimensions in the field\n");
(void) fprintf( AVS_fld_fileptr, "dim1=%d\t\t\t# number of x-elements\n", n1[X]);
(void) fprintf( AVS_fld_fileptr, "dim2=%d\t\t\t# number of y-elements\n", n1[Y]);
(void) fprintf( AVS_fld_fileptr, "dim3=%d\t\t\t# number of z-elements\n", n1[Z]);
(void) fprintf( AVS_fld_fileptr, "nspace=3\t\t# number of physical coordinates per point\n");
(void) fprintf( AVS_fld_fileptr, "veclen=%d\t\t# number of affinity values at each point\n", num_maps);
(void) fprintf( AVS_fld_fileptr, "data=float\t\t# data type (byte, integer, float, double)\n");
(void) fprintf( AVS_fld_fileptr, "field=uniform\t\t# field type (uniform, rectilinear, irregular)\n");
for (int i = 0;  i < XYZ;  i++) {
    (void) fprintf( AVS_fld_fileptr, "coord %d file=%s filetype=ascii offset=%d\n", (i + 1), xyz_filename, (i*2));
}
for (int i = 0;  i < num_atom_maps;  i++) {
    (void) fprintf( AVS_fld_fileptr, "label=%s-affinity\t# component label for variable %d\n", gridmap[i].type, (i + 1));
} /* i */
if(elecPE>=0)
(void) fprintf( AVS_fld_fileptr, "label=Electrostatics\t# component label for variable %d\n", elecPE+1);
if(dsolvPE>=0)
(void) fprintf( AVS_fld_fileptr, "label=Desolvation\t# component label for variable %d\n", dsolvPE+1);
if (floating_grid) {
    (void) fprintf( AVS_fld_fileptr, "label=Floating_Grid\t# component label for variable %d\n", num_maps);
}
(void) fprintf( AVS_fld_fileptr, "#\n# location of affinity grid files and how to read them\n#\n");
for (int i = 0;  i < num_atom_maps;  i++) {
    (void) fprintf( AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", (i + 1), gridmap[i].map_filename);
}
if(elecPE>=0)
(void) fprintf( AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_atom_maps + 1, gridmap[elecPE].map_filename);
if(dsolvPE>=0)
(void) fprintf( AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_atom_maps + 2, gridmap[dsolvPE].map_filename);
if (floating_grid) {
    (void) fprintf( AVS_fld_fileptr, "variable %d file=%s filetype=ascii skip=6\n", num_maps+1, floating_grid_filename);
}
(void) fclose( AVS_fld_fileptr);
AVS_fld_fileptr = NULL;

#ifdef BOINCCOMPOUND
 boinc_fraction_done(0.1);
#endif

(void) fprintf( logFile, "\n\nCalculating Pairwise Interaction Energies\n");
(void) fprintf( logFile,   "=========================================\n\n");

/**************************************************
 * do the map stuff here: 
 * set up xA, xB, npb_r, npb_eps and hbonder 
 * before this pt
 **************************************************/
if (use_vina_potential) {
    (void) fprintf( logFile,   "Use vina potential is %d\n\n", use_vina_potential);
}
    float map_Rij;
    //from autodock_vina_1_1_1/src/main.cpp,line 394
    float wt_gauss1 = -0.035579;
    float wt_gauss2 = -0.005156;
    float wt_repulsion = 0.840245;
    float wt_hydrogen = -0.587439;
    float wt_hydrophobic = -0.035069; //C_H,F_H,Cl_H,Br_H,I_H
    // xs_vdw_radii from atom_constants.h now in read_parameter_library.cc
    //float C_H = 1.9;//C_P
    //float N_P = 1.8;//N_D,N_A,N_DA
    //float O_P = 1.7;//O_D,O_A,O_DA
    //float S_P = 2.0;
    //float P_P = 2.1;
    //float F_H = 1.5;
    //float Cl_H = 1.8;
    //float Br_H = 2.0;
    //float I_H = 2.2;
    //float Met_D = 1.2; //metal_donor:Mg,Mn,Zn,Ca,Fe,Cl,Br
    //float Met_non_ad = 1.75;//metal_non_ad:Cu,Fe,Na,K,Hg,Co,U,Cd,Ni
    double rddist; 
    double delta_e = 0.0;
    //                                                               ia_dist
    //    interatom_distance                                   |.................|
    //    interatom_distance - xs_radius(t1) -xs_radius(t2)    .--|........|-----.  
    //                                                       at1               at2
    //vina distance from current gridpt to atom ia:        xs_rad1  rddist   xs_rad2
    //0. process receptor to setup type[ia], coords[ia], xs_rad[ia]
    //1. setup map types
    //2. setup the energy_lookup tables (NOTE: replaced by the "et" structure)
    //3. loop over all the maps 
    //4.     loop over all pts in current map_ia
    //5.        loop over all the receptor atoms adding to this pt
    //             rdist: interatom_distance from atom coords to current grid pt
    //                rddist based on types: ia_dist - (xs_rad1 + xs_rad2)
    //             e_attractive:
    //             delta_e = rgauss1*exp(-((rddist)/0.5)**2) + rgauss2*exp(-((rddist-3.)/2.)**2);
    //             energy_lookup[i][indx_r][ia] += delta_e
    //             e_repulsive:
    //             if (rddist<0.0){
    //             delta_e = rrepulsive*(rddist**2);
    //             energy_lookup[i][indx_r][ia] += delta_e;
    //             };
    //             e_hbond:
    //             (1)set ihb from types; it is set to 1 if pair of types is suitable for hbond int.
    //             ihb = 0;
    //             if (ihb>0){
    //               if (rddist<0.7)
    //                  delta_e = 1*weight_hydrogen;
    //                  energy_lookup[i][indx_r][ia] += delta_e;
    //               if ((-0.7<rddist) && (rddist<0.))
    //                  delta_e =(rddist/0.7)*weight_hydrogen;
    //                  energy_lookup[i][indx_r][ia] -= delta_e;
    //             }
    //             e_hydrophobic:
    //             //energy_lookup[atom_type[ia]][indx_r][map_ia]+= e_hphob
    //             (1)TODO: set ihb from types; it is set to 1 if pair of types is suitable for hydrophobic int.
    //             ihb = 0
    //             if (rddist<0.5){
    //                  energy_lookup[i][indx_r][ia]+= 1*weight_hydrophobic;}
    //             else if ((0.5<rddist)&& (rdist<1.5)) {
    //                  energy_lookup[i][indx_r][ia]+=(0.5-rddist)*weight_hydrophobic;
    //             };
    // to use in filling out the gridmap
    //????      gridmap[map_index].energy += energy_lookup[i][indx_r][i];
    //6.     output this map
    //

for (int ia=0; ia<num_atom_maps; ia++){
    if (gridmap[ia].is_covalent == FALSE) {
        /* i is the index of the receptor atom type, that 
         * the ia type ligand probe will interact with. */ /* GPF_MAP */
#ifdef DEBUG
        printf("receptor_types_ct=%d\n", receptor_types_ct);
#endif
        ParameterEntry * lig_parm = apm_find(ligand_types[ia]);
        for (int i = 0;  i < receptor_types_ct;  i++) {
	    double cA,cB,Rij,epsij; int xA,xB; /* local copy of gridmap info, for code simplicity */
            /*for each receptor_type*/
            cA = gridmap[ia].cA[i];
            cB = gridmap[ia].cB[i];
            xA = gridmap[ia].xA[i];
            xB = gridmap[ia].xB[i];
            Rij = gridmap[ia].nbp_r[i];
            epsij = gridmap[ia].nbp_eps[i];
            ParameterEntry * rec_parm = apm_find(receptor_types[i]);
            if (use_vina_potential){//from vina: atom_constants.h
#ifdef DEBUG
                printf("in use_vina_potential loop ia=%d, i=%d\n", ia,i);
                fprintf(logFile, "in use_vina_potential loop ia=%d, i=%d\n", ia,i);
#endif
                // get xs_radius for this probe type from parameter_library
                Rij = lig_parm->xs_radius; //see read_parameter_library
                map_Rij = rec_parm->xs_radius; //see read_parameter_library
                //TODO: add SER-OG,THR-OG, TYR_OH: X(1.2) Cl_H(1.8),Br_H(2.0),I_H(2.2),Met_D(1.2)
                /* loop over distance index, indx_r, from 0 to NEINT (scaled NBC non-bond cutoff) */ /* GPF_MAP */
#ifdef DEBUG
                printf("%d-%d-building  Rij=%6.3lf, map_Rij=%10.8f for %s %s\n",ia,i, Rij, map_Rij, gridmap[ia].type, ligand_types[ia]);
                (void) fprintf( logFile, "Calculating vina energies for %s-%s interactions (%d, %d).\n", gridmap[ia].type, receptor_types[i], ia, i );
#endif
                for (int indx_r = 1;  indx_r < NEINT;  indx_r++) {
                    r  = angstrom(indx_r);
                    // compute rddist:                                   map_Rij  rddist   Rij
                    //  interatom_distance - xs_radius(t1) -xs_radius(t2)    .--|........|-----.  
                    rddist =  r - (map_Rij + Rij);
                    //use rddist for computing the vina component energies
                    //TODO: replace with functions from vina..
                    //attraction:
                    delta_e = wt_gauss1 * exp(-pow(((rddist)/0.5),2)) + wt_gauss2 * exp(-pow(((rddist-3.)/2.),2));
                    //at distance 'indx_r': interaction of receptor atomtype 'ia' - ligand atomtype 'i'
                    et.e_vdW_Hb[indx_r][i][ia] += delta_e; 
                    //repulsion
                    if (rddist<0){
                        delta_e = wt_repulsion*pow(rddist,2);
                        et.e_vdW_Hb[indx_r][i][ia] += delta_e;
                    }
                    //hbond
                    if (gridmap[ia].hbonder[i]>0){ //check that ia-i must be hbonder
#ifdef DEBUG
                        printf(" processing gridmap= %d-hbonder i= %d\n",ia, i);
#endif
                        if (rddist<=0.7) { //what about EXACTLY 0.7?
                           delta_e = 1*wt_hydrogen;
                           et.e_vdW_Hb[indx_r][i][ia] += delta_e;
                        }
                        if ((-0.7<rddist) && (rddist<=0.)){
                           delta_e =(rddist/0.7)*wt_hydrogen;
                           et.e_vdW_Hb[indx_r][i][ia] -= delta_e;
                        }
                    }
                    // hydrophobic: check using index 'i' compared with 'carbon',
                    // carbon/aromatic_carbon       to non-hbonder
                    //TODO: add support for these other hydrophobic interactions:
                    //if (((i==carbon)||(i==arom_carbon)||(i==fluorine)||(i==chlorine)||(i==bromine)||(i==iodine))
                    //&& ((ia==carbon)||(ia==arom_carbon)||(ia==fluorine)||(ia==chlorine)||(ia==bromine)||(ia==iodine))) 
                    //if (((i==carbon)||(i==arom_carbon)) && ((ia==carbon)||(ia==arom_carbon)))
                    if ((gridmap[ia].hbonder[i]==0)&&(gridmap[i].hbonder[ia]==0)) { //4-3-12:hydrophobic update
                        delta_e = 0.;
                        if (rddist<0.5) {
                           delta_e = 1*wt_hydrophobic;
                        } else if (rddist<1.5){
                           delta_e = (0.5-rddist)*wt_hydrophobic;
                        }
                        et.e_vdW_Hb[indx_r][i][ia] += delta_e;
                    }
                } /*for each distance*/ 
#ifdef DEBUG
                printf("END USE_VINA_POTENTIAL\n");
#endif
             } /* END use_vina_potential*/
            else { //use regular autogrid potential
            /*for each receptor_type get its parms and fill in tables*/
            if (cA==0 && cB==0){
		double tmpconst;
                cA = (tmpconst = epsij / (double)(xA - xB)) * pow( Rij, (double)xA ) * (double)xB;
                cB = tmpconst * pow( Rij, (double)xB ) * (double)xA;
            /*printf("tmpconst = %6.4f, cA = %6.4f, cB = %6.4f\n",tmpconst, cA, cB);*/
            }
            if ( isnan( cA ) ) {
                print_error( logFile, FATAL_ERROR, "Van der Waals coefficient cA is not a number.  AutoGrid must exit." );
            }
            if ( isnan( cB ) ) {
                print_error( logFile, FATAL_ERROR, "Van der Waals coefficient cB is not a number.  AutoGrid must exit." );
            }
            if ( xA == 0 ) {
                print_error( logFile, FATAL_ERROR, "Van der Waals exponent xA is 0.  AutoGrid must exit." );
            }
            if ( xB == 0 ) {
                print_error( logFile, FATAL_ERROR, "Van der Waals exponent xB is 0.  AutoGrid must exit." );
            }
            (void) fprintf( logFile, "\n             %9.1lf       %9.1lf \n", cA, cB);
            (void) fprintf( logFile, "    E    =  -----------  -  -----------\n");
            (void) fprintf( logFile, "     %s, %s         %2d              %2d\n", gridmap[ia].type, receptor_types[i], xA, xB);
            (void) fprintf( logFile, "                r               r \n\n");
            /* loop over distance index, indx_r, from 0 to max(NEINT,NDIEL) */ /* GPF_MAP */
            (void) fprintf( logFile, "Calculating energies for %s-%s interactions.\n", gridmap[ia].type, receptor_types[i] );

	    // do up to non-bond cutoff distance
            for (int indx_r = 1;  indx_r < NEINT;  indx_r++) {
		double r, rA, rB;
                r  = angstrom(indx_r);
                rA = pow( r, (double) xA);
                rB = pow( r, (double) xB);
		// these should probably be assert()s - MP TODO
		if(i>=NUM_RECEPTOR_TYPES) printf("i>=%d %d\n", NUM_RECEPTOR_TYPES,i);
		if(ia>=MAX_MAPS) printf("ia>=%d %d\n", MAX_MAPS, ia);
                et.e_vdW_Hb[indx_r][i][ia] = min(EINTCLAMP, (cA/rA - cB/rB));
            } /*for each distance*/ 
            et.e_vdW_Hb[0][i][ia] = EINTCLAMP;
            et.e_vdW_Hb[NEINT-1][i][ia] = 0.;

#ifdef PRINT_BEFORE_SMOOTHING
            /*PRINT OUT INITIAL VALUES before smoothing here */
            (void) fprintf( logFile, "before smoothing\n  r ");
            for (int iat = 0;  iat < receptor_types_ct;  iat++) {
                (void) fprintf( logFile, "    %s    ", receptor_types[iat]);
            } 
            (void) fprintf( logFile, "\n ___");
            for (int iat = 0;  iat < receptor_types_ct;  iat++) {
                (void) fprintf( logFile, " ________");
            }
            (void) fprintf( logFile, "\n");
	    for (int j = 0;  j <= min(500,NEINT);  j += 10) {
                (void) fprintf( logFile, "%4.1lf", angstrom(j));
                for (int iat = 0;  iat < receptor_types_ct;  iat++) {
                    (void) fprintf( logFile, (et.e_vdW_Hb[j][iat][ia]<100000.)?"%9.2lf":"%9.2lg", et.e_vdW_Hb[j][iat][ia]);
				} 
                (void) fprintf( logFile, "\n");
            } 
            (void) fprintf( logFile, "\n");
#endif

            /* smooth with min function */ /* GPF_MAP */

	    /* Angstrom is divided by A_DIV in look-up table. */
	    /* Typical value of r_smooth is 0.5 Angstroms  */
	    /* so i_smooth = 0.5 * 100. / 2 = 25 */
	    i_smooth = (int) (r_smooth*A_DIV/2.);
            if (i_smooth > 0) {
                for (int indx_r = 0;  indx_r < NEINT;  indx_r++) {
                    energy_smooth[indx_r] = 100000.;
                    for (int j = max(0, indx_r - i_smooth);  j < min(NEINT, indx_r + i_smooth + 1);  j++) {
                      energy_smooth[indx_r] = min(energy_smooth[indx_r], et.e_vdW_Hb[j][i][ia]);
                    }
                }
                for (int indx_r = 0;  indx_r < NEINT;  indx_r++) {
                    et.e_vdW_Hb[indx_r][i][ia] = energy_smooth[indx_r];
                }
            } /* endif smoothing */
        } /* end regular autogrid potential */
        } /* for i in receptor types: build energy table for this map */

       /*
        * Print out a table, of distance versus energy...
        */ /* GPF_MAP */
        (void) fprintf( logFile, "\n\nFinding the lowest pairwise interaction energy within %.1f Angstrom (\"smoothing\").\n\n  r ", r_smooth);
        for (int iat = 0;  iat < receptor_types_ct;  iat++) {
            (void) fprintf( logFile, "    %s    ", receptor_types[iat]);
            /*(void) fprintf( logFile, "    %c    ", receptor_atom_type_string[iat]);*/
        } /* iat */
        (void) fprintf( logFile, "\n ___");
        for (int iat = 0;  iat < receptor_types_ct;  iat++) {
            (void) fprintf( logFile, " ________");
        } /* iat */
        (void) fprintf( logFile, "\n");
        for (int j = 0;  j <= min(500,NEINT);  j += 10) {
            (void) fprintf( logFile, "%4.1lf", angstrom(j));
            for (int iat = 0;  iat < receptor_types_ct;  iat++) {
                (void) fprintf( logFile, (et.e_vdW_Hb[j][iat][ia]<100000.)?"%9.2lf":"%9.2lg", et.e_vdW_Hb[j][iat][ia]);
		} /* iat */
            (void) fprintf( logFile, "\n");
        } /* j */
        (void) fprintf( logFile, "\n");
        (void) fprintf( logFile, "\n\nEnergyTable:\n");
        (void) fprintf( logFile, "Finding the lowest pairwise interaction energy within %.1f Angstrom (\"smoothing\").\n\n  r ", r_smooth);
        for (int iat = 0;  iat < receptor_types_ct;  iat++) {
            (void) fprintf( logFile, "    %s    ", receptor_types[iat]);
            /*(void) fprintf( logFile, "    %c    ", receptor_atom_type_string[iat]);*/
        } /* iat */
        (void) fprintf( logFile, "\n ___");
        for (int iat = 0;  iat < receptor_types_ct;  iat++) {
            (void) fprintf( logFile, " ________");
        } /* iat */
        (void) fprintf( logFile, "\n");
        for (int j = 0;  j <= min(500,NEINT);  j += 10) {
            (void) fprintf( logFile, "%4.1lf", angstrom(j));
            for (int iat = 0;  iat < receptor_types_ct;  iat++) {
                (void) fprintf( logFile, (et.e_vdW_Hb[j][iat][ia]<100000.)?"%9.2lf":"%9.2lg", et.e_vdW_Hb[j][iat][ia]);
		} /* iat */
            (void) fprintf( logFile, "\n");
        } /* j */
        (void) fprintf( logFile, "\n");
    } else {
        /* parsing for intnbp not needed for covalent maps */
        (void) fprintf( logFile, "\nAny internal non-bonded parameters will be ignored for this map, since this is a covalent map.\n");
    } /*end of else parsing intnbp*/
} /*end of loop over all the atom maps*/

/* exponential function for receptor and ligand desolvation */
/* note: the solvation term ranges beyond the non-bond cutoff 
 * and will not be smoothed 
 */
double sigma = 3.6;
for (int indx_r = 1;  indx_r < NDIEL;  indx_r++) {
     r  = angstrom(indx_r);
     et.sol_fn[indx_r] = AD4.coeff_desolv * exp(-sq(r)/(2.*sq(sigma)));
}

/**************************************************
 * Loop over all RECEPTOR atoms to
 * calculate bond vectors for directional H-bonds
 **************************************************/
 //setup the canned atom types here....
//at this point set up hydrogen, carbon, oxygen and nitrogen
hydrogen = get_rec_index("HD");
nonHB_hydrogen = get_rec_index("H");
carbon = get_rec_index("C");
arom_carbon = get_rec_index("A");
oxygen = get_rec_index("OA");
nitrogen = get_rec_index("NA");
nonHB_nitrogen = get_rec_index("N");
sulphur = get_rec_index("SA");
nonHB_sulphur = get_rec_index("S");


if (not use_vina_potential){
/********************************************
 * Start bond vector loop
 ********************************************/
int i1=0, i2=0, i3=0;
for (int ia=0; ia<num_receptor_atoms; ia++) {  /*** ia = i_receptor_atom_a ***/
    double d[XYZ], dc[XYZ];
    disorder[ia] = FALSE;  /* initialize disorder flag. */
    warned = 'F';
    /*
     * Set scan limits looking for bonded atoms
     */
    int from = max(ia-20, 0);
    int to   = min(ia + 20, num_receptor_atoms-1);

    /*
     * If 'ia' is a hydrogen atom, it could be a
     * RECEPTOR hydrogen-BOND DONOR,
     */
    /*8:CHANGE HERE: fix the atom_type vs atom_types problem in following*/
    if ((int)hbond[ia] == 2) { /*D1 hydrogen bond donor*/

        for ( int ib = from; ib <= to; ib++) if (ib != ia) { /* ib = i_receptor_atom_b */
                /*
                 * =>  NH-> or OH->
                 */
                /*if ((atom_type[ib] == nitrogen) || (atom_type[ib]==nonHB_nitrogen) ||(atom_type[ib] == oxygen)||(atom_type[ib] == sulphur)||(atom_type[ib]==nonHB_sulphur)) {*/
        
                    /*
                     * Calculate the square of the N-H or O-H bond distance, rd2,
                     *                            ib-ia  ib-ia
                     */
                    for (int i = 0;  i < XYZ;  i++) {
                        d[i] = coord[ia][i] - coord[ib][i];
                    }
                    rd2 = sq( d[X] ) + sq( d[Y] ) + sq( d[Z]);
                    /*
                     * If ia & ib are less than 1.3 A apart -- they are covalently bonded,
                     */
                    if (rd2 < sq(1.378)) { /*INCREASED for H-S bonds*/
                        if (rd2 < APPROX_ZERO) {
                            if (rd2 == 0.) {
                                (void) sprintf ( message, "While calculating an H-O or H-N bond vector...\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                                print_error( logFile, WARNING, message );
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
                        for (int i = 0;  i < XYZ;  i++) {
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
        for ( ib = from; ib <= to; ib++) if ( ib != ia ) {
                float rd2 = 0.;

                for (int i = 0;  i < XYZ;  i++) {
                    dc[i] = coord[ia][i] - coord[ib][i];
                    rd2 += sq( dc[i]);
                }

                if (((rd2 < sq(1.9)) && ((atom_type[ib] != hydrogen)&&(atom_type[ib]!=nonHB_hydrogen))) ||
                    ((rd2 < sq(1.3)) && ((atom_type[ib] == hydrogen)||(atom_type[ib]==nonHB_hydrogen)))) {
                    if (nbond == 2) {
                        (void) sprintf( message, "Found an H-bonding atom with three bonded atoms, atom serial number %d\n", ia + 1);
                        print_error( logFile, WARNING, message );
                    }
                    else if (nbond == 1) {
                        nbond = 2;
                        i2 = ib;
                    }
                    else if (nbond == 0) {
                        nbond = 1;
                        i1 = ib;
                    }
                }
        } /*ib-loop*/

        /* if no bonds, something is wrong */

        if (nbond == 0) {
            (void) sprintf( message, "Oxygen atom found with no bonded atoms, atom serial number %d, atom_type %d\n", ia + 1, atom_type[ia]);
            print_error( logFile, WARNING, message );
        }

        /* one bond: Carbonyl Oxygen O=C-X */

        if (nbond == 1) {

            /* calculate normalized carbonyl bond vector rvector[ia][] */

            float rd2 = 0.;
            for (int i = 0;  i < XYZ;  i++) {
                rvector[ia][i] = coord[ia][i]-coord[i1][i];
                rd2 += sq(rvector[ia][i]);
            }
            if (rd2 < APPROX_ZERO) {
                if ((rd2 == 0.) && (warned == 'F')) {
                    (void) sprintf ( message, "At receptor carbonyl oxygen i1:\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, i1 + 1);
                    print_error( logFile, WARNING, message );
                    warned = 'T';
                }
                rd2 = APPROX_ZERO;
            }
            inv_rd = 1./sqrt(rd2);
            for (int i = 0;  i < XYZ;  i++) {
                rvector[ia][i] *= inv_rd;
            }

            /* find a second atom (i2) bonded to carbonyl carbon (i1) */
            for ( i2 = from; i2 <= to; i2++) if (i2 != i1 && i2 != ia) {
                    rd2 = 0.;
                    for (int i = 0;  i < XYZ;  i++) {
                        dc[i] = coord[i1][i] - coord[i2][i]; /*NEW*/
                        rd2 += sq( dc[i]);
                    }
                    if ( (atom_type[i2] == hydrogen && rd2 < sq(1.3)) ||
                         (atom_type[i2] != hydrogen && rd2 < sq(1.7))) {

                        /* found one */
                        /* d[i] vector from carbon to second atom */
                        rd2 = 0.;
                        for (int i = 0;  i < XYZ;  i++) {
                            d[i] = coord[i2][i]-coord[i1][i];
                            rd2 += sq( d[i]);
                        }
                        if (rd2 < APPROX_ZERO) {
                            if ((rd2 == 0.) && (warned == 'F')) {
                                (void) sprintf ( message, "At receptor carbonyl oxygen i2:\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", i1 + 1, i2 + 1);
                                print_error( logFile, WARNING, message );
                                warned = 'T';
                            }
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1./sqrt(rd2);
                        for (int i = 0;  i < XYZ;  i++) {
                            d[i] *= inv_rd;
                        }

                        /* C=O cross C-X gives the lone pair plane normal */
                        rvector2[ia][0] = rvector[ia][1]*d[2] - rvector[ia][2]*d[1];
                        rvector2[ia][1] = rvector[ia][2]*d[0] - rvector[ia][0]*d[2];
                        rvector2[ia][2] = rvector[ia][0]*d[1] - rvector[ia][1]*d[0];
                        rd2 = 0.;
                        for (int i = 0;  i < XYZ;  i++) {
                            rd2 += sq(rvector2[ia][i]);
                        }
                        if (rd2 < APPROX_ZERO) {
                            if ((rd2 == 0.) && (warned == 'F')) {
                                (void) sprintf ( message, "At receptor carbonyl oxygen C==O cross C-X:\n:Attempt to divide by zero was just prevented.\n\n" );
                                print_error( logFile, WARNING, message );
                                warned = 'T';
                            }
                            rd2 = APPROX_ZERO;
                        }
                        inv_rd = 1./sqrt(rd2);
                        for (int i = 0;  i < XYZ;  i++) {
                            rvector2[ia][i] *= inv_rd;
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
                for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i] - coord[ib][i];
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) sprintf ( message, "At disordered hydroxyl:\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia + 1, ib + 1);
                        print_error( logFile, WARNING, message );
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1./sqrt(rd2);
                for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] *= inv_rd;
                }

            } else {

                /* not a disordered hydroxyl */
                /* normalized X1 to X2 vector, defines lone pair plane */

                rd2 = 0.;
                for (int i = 0;  i < XYZ;  i++) {
                    rvector2[ia][i] = coord[i2][i] - coord[i1][i];
                    rd2 += sq(rvector2[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) sprintf ( message, "At not disordered hydroxyl:\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", i1 + 1, i2 + 1);
                        print_error( logFile, WARNING, message );
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1./sqrt(rd2);
                for (int i = 0;  i < XYZ;  i++) {
                    rvector2[ia][i] *= inv_rd;
                }

                /* vector pointing between the lone pairs:
                ** front of the vector is the oxygen atom,
                ** X1->O vector dotted with normalized X1->X2 vector plus
                ** coords of X1 gives the point on the X1-X2 line for the
                ** back of the vector.
                */
                double rdot = 0.;
                for (int i = 0;  i < XYZ;  i++) {
                    rdot += (coord[ia][i] - coord[i1][i]) * rvector2[ia][i] ;
                }
                rd2 = 0.;
                for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i] - ( (rdot*rvector2[ia][i]) + coord[i1][i] ) ;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) sprintf ( message, "At disordered hydroxyl lone pair vector:\nAttempt to divide by zero was just prevented.\n\n" );
                        print_error( logFile, WARNING, message );
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1./sqrt(rd2);
                for (int i = 0;  i < XYZ;  i++) {
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
        for ( ib = from; ib <= to; ib++) if ( ib != ia) {
                float rd2 = 0.;

                for (int i = 0;  i < XYZ;  i++) {
                    dc[i] = coord[ia][i] - coord[ib][i];
                    rd2 += sq( dc[i] );
                }

                /*
                    for (int i = 0;  i < XYZ;  i++) {
                        rd2 += sq(coord[ia][i] - coord[ib][i]);
                    }
                */
                if (((rd2 < sq(1.7)) && ((atom_type[ib] != hydrogen)&&(atom_type[ib]!=nonHB_hydrogen))) || 
                    ((rd2 < sq(1.3)) && ((atom_type[ib] == hydrogen)||(atom_type[ib]==nonHB_hydrogen)))) {
                    if (nbond == 2) {
                        nbond = 3;
                        i3 = ib;
                    }
                    else if (nbond == 1) {
                        nbond = 2;
                        i2 = ib;
                    }
                    else if (nbond == 0) {
                        nbond = 1;
                        i1 = ib;
                    }
                }
        } /*ib-loop*/

        /* if no bonds, something is wrong */

        if (nbond == 0) {
            (void) sprintf( message, "Nitrogen atom found with no bonded atoms, atom serial number %d\n",ia);
            print_error( logFile, WARNING, message );
        }

        /* one bond: Azide Nitrogen :N=C-X */

        if (nbond == 1) {

            /* calculate normalized N=C bond vector rvector[ia][] */

            rd2 = 0.;
            for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i]-coord[i1][i];
                    rd2 += sq(rvector[ia][i]);
            }
            if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) sprintf ( message, "At azide nitrogen 1 bond:\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, i1);
                        print_error( logFile, WARNING, message );
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
            }
            inv_rd = 1./sqrt(rd2);
            for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] *= inv_rd;
            }
        } /* endif nbond==1 */

        /* two bonds: X1-N=X2 */
        if (nbond == 2) {
                /* normalized vector from Nitrogen to midpoint between X1 and X2 */

                rd2 = 0.;
                for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i]-(coord[i2][i]+coord[i1][i])/2.;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) sprintf ( message, "At azide nitrogen 2 bonds:\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        print_error( logFile, WARNING, message );
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1./sqrt(rd2);
                for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] *= inv_rd;
                }

        }  /* end two bonds for nitrogen*/

        /* three bonds: X1,X2,X3 */
        if (nbond == 3) {
                /* normalized vector from Nitrogen to midpoint between X1, X2, and X3 */

                rd2 = 0.;
                for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] = coord[ia][i]-(coord[i1][i]+coord[i2][i]+coord[i3][i])/3.;
                    rd2 += sq(rvector[ia][i]);
                }
                if (rd2 < APPROX_ZERO) {
                    if ((rd2 == 0.) && (warned == 'F')) {
                        (void) sprintf ( message, "At azide nitrogen 3 bonds:\nAttempt to divide by zero was just prevented.\nAre the coordinates of atoms %d and %d the same?\n\n", ia, ib);
                        print_error( logFile, WARNING, message );
                        warned = 'T';
                    }
                    rd2 = APPROX_ZERO;
                }
                inv_rd = 1./sqrt(rd2);
                for (int i = 0;  i < XYZ;  i++) {
                    rvector[ia][i] *= inv_rd;
                }

        }  /* end three bonds for Nitrogen */
    /* endNEW directional N Acceptor */
    } /* end test for atom type */

} /* Do Next receptor atom... */

/********************************************
 * End bond vector loop
 ********************************************/
} /* not use_vina_potential*/

for (int k = 0;  k < num_atom_maps + 1;  k++) {
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
for (int k = 0;  k < num_maps;  k++) {
    (void) fprintf( gridmap[k].map_fileptr, "GRID_PARAMETER_FILE %s\n", grid_param_fn );
    (void) fprintf( gridmap[k].map_fileptr, "GRID_DATA_FILE %s\n", AVS_fld_filename);
    (void) fprintf( gridmap[k].map_fileptr, "MACROMOLECULE %s\n", receptor_filename);
    (void) fprintf( gridmap[k].map_fileptr, "SPACING %.3lf\n", spacing);
    (void) fprintf( gridmap[k].map_fileptr, "NELEMENTS %d %d %d\n", nelements[X], nelements[Y], nelements[Z]);
    (void) fprintf( gridmap[k].map_fileptr, "CENTER %.3lf %.3lf %.3lf\n", center[X], center[Y], center[Z]);
    fflush(gridmap[k].map_fileptr);
}
if (floating_grid) {
    (void) fprintf( floating_grid_fileptr, "GRID_PARAMETER_FILE %s\n", grid_param_fn );
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
 fflush( logFile);

/*
 * Iterate over all grid points, Z( Y ( X ) ) (X is fastest)...
 */
#pragma omp parallel for ordered schedule(dynamic)
for (iz=0;iz<n1[Z];iz++) {
	Clock      grd_start;
	Clock      grd_end;
	struct tms tms_grd_start;
	struct tms tms_grd_end;
    grd_start = times( &tms_grd_start); // this is per-plane timing 

    double r_min=BIG; /* for floating_grid only MP TODO NEEDS to be plane-size array */
    int icoord[XYZ];
    icoord[Z] = iz - ne[Z];

    double c[XYZ];
    /*
     *  c[0:2] contains the current grid point coordinates (angstroms)
     */
    c[Z] = ((double)icoord[Z]) * spacing;
    int thread = omp_get_thread_num(); // this thread's working storage

if(outlev>LOGRUNV)
fprintf(logFile, "Starting plane iz=%d icoord=%d z=%8.2f thread=%d\n", iz,icoord[Z],c[Z],thread);fflush(logFile);

    /* Zero-out non-covalent maps' energy for this plane */
    // MP TODO not necessary?
    for (int j = 0;  j < num_maps ;  j++) {
                if (!gridmap[j].is_covalent) {
			bzero(gridmap[j].energy[thread], n1[X]*n1[Y]*sizeof(double)); 
		}
    }

if(outlev>LOGRUNVV)
fprintf(logFile, "Zeroed plane/thread=%d iz=%d icoord=%d z=%8.2f\n", thread,iz,icoord[Z],c[Z]);fflush(logFile);
    for (int iy=0; iy<n1[Y]; iy++) {
        icoord[Y] = iy - ne[Y];
        c[Y] = ((double)icoord[Y]) * spacing;

        for (int ix=0; ix<n1[X]; ix++) {
	    // local storage for one grid point calculation:
	    float fcc[XYZ]; // for USE_BHTREE
	    bool hbondflag[MAX_MAPS];
            double hbondmin[MAX_MAPS], hbondmax[MAX_MAPS];
	    double rminH; int closestH;
	    double d[XYZ];
double inv_rd, rd2, r; /* re, r2, rd, */
double racc, rdon;
double t0, ti;

            icoord[X] = ix - ne[X];
            c[X] = ((double)icoord[X]) * spacing;
	    for ( int i = 0; i< XYZ; i++) fcc[i] = c[i]; // for USE_BHTREE
	    
	    if(outlev>=LOGRUNVVV) fprintf(logFile, "grid point %3d %3d %3d xyz(%5.2f, %5.2f, %5.2f)\n",
		icoord[X], icoord[Y], icoord[Z], c[X], c[Y], c[Z]);

	    /* handle covalent map(s). */
            for (int j = 0;  j < num_maps ;  j++)  if (gridmap[j].is_covalent) {
                    /* Calculate rcov, the distance from the current
                     * grid point, c, to the covalent attachment point, covpos */
		    double rcov; 
                    for (int ii = 0;  ii < XYZ;  ii++) {
                        d[ii]  = covpos[ii] - c[ii];
                    }
                    rcov = hypotenuse( d[X], d[Y], d[Z] ) / covhalfwidth;
                    if (rcov < APPROX_ZERO) {
                        rcov = APPROX_ZERO;
                    }
                    gridmap[j].energy[thread][mapindex(icoord[X],icoord[Y])] 
		      = covbarrier * (1. - exp(ln_half * rcov * rcov));
            }

            /* Initialize Min Hbond variables  for each new point*/
            if (not use_vina_potential) //if vina_potential skip hbond stuff
            for (int map_index = 0; map_index < num_atom_maps; map_index++){
                hbondmin[map_index] = 999999.;
                hbondmax[map_index] = -999999.;
                hbondflag[map_index] = FALSE;
	    }
            
// M Sanner 2015 use BHTree to find atoms very close to grid point: write all atom grids' energies as "HIGH"
//  write electrostatic and desolvation maps as 0
#ifdef USE_BHTREE
	    if (! map_receptor_interior) {
	      // check for collision with any receptor atoms
 	      // if collision, modify energy, write all maps' values, skip to next grid point.
	      int bhTreeNbIndices = findBHcloseAtomsdist(bht, fcc, 
		BH_collision_dist, 
		closeAtomsIndices[thread], closeAtomsDistances[thread],
		num_receptor_atoms);

	      if (bhTreeNbIndices > 0) {
		  for (int j = 0;  j < num_maps;  j++) {
		    if(!gridmap[j].is_covalent)  
                      gridmap[j].energy[thread][mapindex(ix,iy)] =
                      (j==elecPE||j==dsolvPE)?0:INTERIOR_VALUE;
		    }
		  continue; // next icoord[X]
		  }
	    }

#endif

		

            /*
             *  Loop 1 of 2:
	     *  Performed only if electrostatic or floating maps are requested.
             *  Do all Receptor (protein, DNA, etc.) atoms, regardless of distance cutoff
             */
            if(elecPE>=0 || floating_grid) for (int ia = 0;  ia < num_receptor_atoms;  ia++) {
                /*
                 *  Get distance, r, from current grid point, c, to this receptor atom, coord,
                 */
                for (int i = 0;  i < XYZ;  i++) {
                    d[i]  = coord[ia][i] - c[i];
                }
                r = hypotenuse( d[X], d[Y], d[Z]);
                if (r < APPROX_ZERO)  r = APPROX_ZERO;
                inv_r  = 1./r;
                inv_rmax = 1./max(r, 0.5);

                /* make sure lookup index is in the tables */
                int indx_r = min(lookup(r), NDIEL-1);

                if (floating_grid) {
                    /* Calculate the so-called "Floating Grid"... */
                    r_min = min(r, r_min);
                }

                /* elecPE is the next-to-last last grid map, i.e. electrostatics */
                /* if use_vina_potential, electPE is -1 */
                if ((not use_vina_potential) && elecPE>=0) { 
                if (dddiel) {
                    /* Distance-dependent dielectric... */
                    /*apply the estat forcefield coefficient/weight here */
                    gridmap[elecPE].energy[thread][mapindex(ix,iy)]  
                     += charge[ia] *inv_rmax * et.r_epsilon_fn[indx_r] * AD4.coeff_estat;
                } else {
                    /* Constant dielectric... */
                    gridmap[elecPE].energy[thread][mapindex(ix,iy)]  
                     += charge[ia] * inv_rmax * invdielcal * AD4.coeff_estat;
                }
		}
            }/* ia loop, over all receptor atoms, no distance cutoff... */

	    // M Pique Oct 2015 TODO combine this loop with one above
	    /* NEW2: Find Closest Hbond */
            rminH=999999.;
            closestH= -1; // none found yet
            if(num_atom_maps>0) for (int inh = 0;  inh < nhbond_12;  inh++) {
		int ia = hbond_12[inh];
                // if ((hbond[ia]==1)||(hbond[ia]==2)||(hbond[ia]==6))  DS or D1 or AD // N3P: directionality for AD not required, right?
                //if ((hbond[ia]==1)||(hbond[ia]==2))  {/*DS or D1 or AD*/
		    bool breakout=FALSE;
		    double d[XYZ];
                    for (int i = 0;  i < XYZ;  i++) { 
                        d[i]  = coord[ia][i] - c[i]; 
			if(fabs(d[i])>rminH) {
				breakout=TRUE;
				break;
			}
                    }
		    if(breakout) continue;
                    r = hypotenuse( d[X],d[Y],d[Z] );
                    if (r < rminH) {
                        rminH = r;
                        closestH = ia; 
                    }
                //} /* Hydrogen test */
            } /* inh loop */

	    /* MPique declare error if no closestH was found */
	    if(closestH<0) print_error( logFile, FATAL_ERROR, "no closestH atom was found");
            /* END NEW2: Find Min Hbond */

	    if(num_atom_maps>0 || dsolvPE>=0) {   // huge block only invoked if atom affinity or desolvation maps requested...
	    /* Loop 2 of 2: consider only atoms within distance cutoff */
#ifdef USE_BHTREE
            int bhTreeNbIndices = findBHcloseAtomsdist(bht, fcc, BH_cutoff_dist, 
		closeAtomsIndices[thread], closeAtomsDistances[thread], num_receptor_atoms);

	    if(outlev>=LOGRUNVVV) fprintf(logFile, "  bhTreeNbIndices= %3d BH_cutoff_dist= %.2f\n",
		bhTreeNbIndices, BH_cutoff_dist);
            for (int ibh = 0;  ibh < bhTreeNbIndices;  ibh++) {
	      int ia = closeAtomsIndices[thread][ibh];
	      r = closeAtomsDistances[thread][ibh];

#else
            for (int ia = 0;  ia < num_receptor_atoms;  ia++) {
#endif
		double dnorm[XYZ]; // normalized grid-point- to - atom vector
		double cos_theta;
                /*
                 *  Get distance, r, from current grid point, c, to this receptor atom, coord,
                 */
                for (int i = 0;  i < XYZ;  i++) {
                    d[i]  = coord[ia][i] - c[i];
                }
#ifndef USE_BHTREE
                r = hypotenuse( d[X], d[Y], d[Z]);
#endif
                /*
                 * If distance from grid point to atom ia is too large,
                 * or if atom is a disordered hydrogen,
                 *   add nothing to the grid-point's non-bond energy;
                 *   just continue to next atom...
                 */
                if ( r > NBC) {
                    continue; /* onto the next atom... */
                }
                if ((atom_type[ia] == hydrogen) && disorder[ia] ) { /* N3P: add check for AD here too?*/
                    continue; /* on to the next atom... */
                }

                if (r < APPROX_ZERO) {
                    r = APPROX_ZERO;
                }
                inv_r  = 1./r;
                inv_rmax = 1./max(r, 0.5);

                for (int i = 0;  i < XYZ;  i++) {
                    dnorm[i] = d[i] * inv_r;
                }
                /* make sure both lookup indices are in the tables */
                int indx_r = min(lookup(r), NDIEL-1);
		int indx_n = min(lookup(r), NEINT-1);

                if (not use_vina_potential){

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
                     *  dnorm[] = Unit vector from current grid pt to ia_th m/m atom.
                     *  cos_theta = d dot rvector == cos(angle) subtended.
                     */
                    for (int i = 0;  i < XYZ;  i++) {
                        cos_theta -= dnorm[i] * rvector[ia][i];
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
			    double tmp;
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
			    double theta;
                            cos_theta = 0.;
                            for (int i = 0;  i < XYZ;  i++) {
                                cos_theta += rvector[closestH][i] * rvector[ia][i];
                            }
                            cos_theta = min(cos_theta, 1.0); 
                            cos_theta = max(cos_theta, -1.0);
                            theta = acos(cos_theta);
                            Hramp = 0.5-0.5*cos(theta * 120./90.);
                        } /* ia test for closestH */
/* END NEW2 calculate dot product of bond vector with bond vector of best hbond */
                    }
                    /* endif (atom_type[ia] == hydrogen) */
                } else if (hbond[ia] == 4) {/*A1*/
                    /* NEW Directional N acceptor */
                    /*
                    **  ia-th macromolecule atom = Nitrogen ( 4 = H )
                    **  calculate rdon for H-bond Donor PROBES at this grid pt.
                    **            ====     ======================
                    */
                    cos_theta = 0.;
                    /*
                    **  dnorm[] = Unit vector from current grid pt to ia_th m/m atom.
                    **  cos_theta = d dot rvector == cos(angle) subtended.
                    */
                    for (int i = 0;  i < XYZ;  i++) {
                        cos_theta -= dnorm[i] * rvector[ia][i];
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
                    for (int i = 0;  i < XYZ;  i++) {
                        cos_theta -= dnorm[i] * rvector[ia][i];
                    }
                    /*
                    ** t0 is the angle out of the lone pair plane, calculated
                    ** as 90 deg - acos (vector to grid point DOT lone pair
                    ** plane normal)
                    */
                    t0 = 0.;
                    for (int i = 0;  i < XYZ;  i++) {
                        t0 += dnorm[i] * rvector2[ia][i];
                    }
                    if (t0 > 1.) {
                        t0 = 1.;
                        /*(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", t0);
                        print_error( logFile, WARNING, message );Feb2012*/ 
                    } else if (t0 < -1.) {
                        t0 = -1.;
                        /*(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", t0);
                        print_error( logFile, WARNING, message );Feb2012*/
                    }
                    t0 = PI_halved - acos(t0);
                    /*
                    ** ti is the angle in the lone pair plane, away from the
                    ** vector between the lone pairs,
                    ** calculated as (grid vector CROSS lone pair plane normal)
                    ** DOT C=O vector - 90 deg
                    */
		    double cross[XYZ];
                    cross[0] = dnorm[1] * rvector2[ia][2] - dnorm[2] * rvector2[ia][1];
                    cross[1] = dnorm[2] * rvector2[ia][0] - dnorm[0] * rvector2[ia][2];
                    cross[2] = dnorm[0] * rvector2[ia][1] - dnorm[1] * rvector2[ia][0];
                    rd2 = sq(cross[0]) + sq(cross[1]) + sq(cross[2]);
                    if (rd2 < APPROX_ZERO) {
                        if ((rd2 == 0.) && (warned == 'F')) {
                            (void) sprintf ( message, "At receptor H-bond acceptor, oxygen:\nAttempt to divide by zero was just prevented.\n\n" );
                            print_error( logFile, WARNING, message );
                            warned = 'T';
                        }
                        rd2 = APPROX_ZERO;
                    }
                    inv_rd = 1./sqrt(rd2);
                    ti = 0.;
                    for (int i = 0;  i < XYZ;  i++) {
                        ti += cross[i] * inv_rd * rvector[ia][i];
                    }

                    /* rdon expressions from Goodford */
                    rdon = 0.;
                    if (cos_theta >= 0.) {
                        if (ti > 1.) {
                            ti = 1.;
                            /*(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", ti);
                            print_error( logFile, WARNING, message );Feb2012*/
                        } else if (ti < -1.) {
                            ti = -1.;
                            /*(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", ti);
                            print_error( logFile, WARNING, message );Feb2012*/
                        }
                        ti = acos(ti) - PI_halved;
                        if (ti < 0.) {
                            ti = -ti;
                        }
                        /* the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);*/
                        rdon = (0.9 + 0.1*sin(ti + ti))*cos(t0);
                    } else if (cos_theta >= -0.34202) {
                        /* 0.34202 = cos (100 deg) */
                        rdon = 562.25*pow(0.116978 - sq(cos_theta), 3.)*cos(t0);
                    }
                /* endif atom_type == oxygen, not disordered */
                } else if ((hbond[ia] == 5) && (disorder[ia] == TRUE)) {/*A2*/
                    /* cylindrically disordered hydroxyl */
		    double theta;
                    cos_theta = 0.;
                    for (int i = 0;  i < XYZ;  i++) {
                        cos_theta -= dnorm[i] * rvector[ia][i];
                    }
                    if (cos_theta > 1.) {
                        cos_theta = 1.;
                        /*(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value greater than 1.\n", cos_theta);
                        print_error( logFile, WARNING, message );Feb2012*/
                    } else if (cos_theta < -1.) {
                        cos_theta = -1.;
                        /*(void) sprintf( message, "I just prevented an attempt to take the arccosine of %f, a value less than -1.\n", cos_theta);
                        print_error( logFile, WARNING, message );Feb2012*/
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
                } /* end atom_type tests used to set rdon and racc */
                }//end  not use_vina_potential

                /*
                 * For each probe atom-type,
                 * Sum pairwise interactions between each probe
                 * at this grid point (c[0:2])
                 * and the current receptor atom, ia...
                 */
                for (int map_index = 0;  map_index < num_atom_maps;  map_index++) {
                    /* We do not want to change the current energy value 
                     * for any covalent maps, make sure iscovalent is
                     * false... */                    

                    if (gridmap[map_index].is_covalent) continue;
                    if ((not use_vina_potential) && (gridmap[map_index].is_hbonder == TRUE)) {
                            /*  current map_index PROBE forms H-bonds... */
			    double rsph;
                            /* rsph ramps in angular dependence for distances with negative energy */
                            rsph = et.e_vdW_Hb[indx_n][atom_type[ia]][map_index]/100.;
                            rsph = max(rsph, 0.);
                            rsph = min(rsph, 1.);
                            if ((gridmap[map_index].hbond==3||gridmap[map_index].hbond==5||gridmap[map_index].hbond==6) /*AS or A2 or AD N3P:modified*/ 
                                      &&(hbond[ia]==1||hbond[ia]==2||hbond[ia]==6)){/*DS or D1 or AD N3P:modified*/
                                  /* PROBE can be an H-BOND ACCEPTOR, */
                                if (disorder[ia] == FALSE ) {
				    gridmap[map_index].energy[thread][mapindex(ix,iy)] 
                                     += et.e_vdW_Hb[indx_n][atom_type[ia]][map_index] * Hramp * (racc + (1. - racc)*rsph);
                                } else {
				    gridmap[map_index].energy[thread][mapindex(ix,iy)] 
                                     += et.e_vdW_Hb[max(0, indx_n - 110)][hydrogen][map_index] * Hramp * (racc + (1. - racc)*rsph);

                                }
                            } else if ((gridmap[map_index].hbond==4 || gridmap[map_index].hbond==6 ) /*A1, AD: N3P: modified*/
                                             &&(hbond[ia]==1||hbond[ia]==2||hbond[ia]==6)) { /*DS,D1, AD: N3P: modified*/
                                hbondmin[map_index] = min( hbondmin[map_index],et.e_vdW_Hb[indx_n][atom_type[ia]][map_index] * (racc+(1.-racc)*rsph));
                                hbondmax[map_index] = max( hbondmax[map_index],et.e_vdW_Hb[indx_n][atom_type[ia]][map_index] * (racc+(1.-racc)*rsph));
                                hbondflag[map_index] = TRUE;
                            } else if ((gridmap[map_index].hbond==1||gridmap[map_index].hbond==2||gridmap[map_index].hbond==6)&& (hbond[ia]>2||hbond[ia]==6)){/*DS,D1 vs AS,A1,A2,AD N3P:modified*/

				double temp_hbond_enrg;
                                /*  PROBE is H-BOND DONOR, */
                                temp_hbond_enrg = et.e_vdW_Hb[indx_n][atom_type[ia]][map_index] * (rdon + (1. - rdon)*rsph);
                                hbondmin[map_index] = min( hbondmin[map_index], temp_hbond_enrg);
                                hbondmax[map_index] = max( hbondmax[map_index], temp_hbond_enrg);
                                hbondflag[map_index] = TRUE;
                            } else {
                                /*  hbonder PROBE-ia cannot form a H-bond..., */
			        gridmap[map_index].energy[thread][mapindex(ix,iy)] 
                                 += et.e_vdW_Hb[indx_n][atom_type[ia]][map_index];
                            }
                        } else { /*end of is_hbonder*/
                            /*  PROBE does not form H-bonds..., */
			    gridmap[map_index].energy[thread][mapindex(ix,iy)] 
                             += et.e_vdW_Hb[indx_n][atom_type[ia]][map_index];
                        }/* end hbonder tests */

                        if (not use_vina_potential){

                        /* add desolvation energy  */
                        /* forcefield desolv coefficient/weight in sol_fn*/
			gridmap[map_index].energy[thread][mapindex(ix,iy)] 
				+= gridmap[map_index].solpar_probe * vol[ia]*et.sol_fn[indx_r] + 
				(solpar[ia]+solpar_q*fabs(charge[ia]))*gridmap[map_index].vol_probe*et.sol_fn[indx_r];
                       }
                } /* end of loop over all num_atom_maps index values */

	    /* compute desolvation map value */
            if ((not use_vina_potential) && dsolvPE>=0){
	        if(outlev>=LOGRUNVVV) fprintf(logFile, " r=%8.3f ", r);
	           if(outlev>=LOGRUNVVV) fprintf(logFile, 
			"  dsolvPE += solpar_q=%8.3f * vol[ia=%3d]=%8.3f * et.sol_fn[indx_r=%4d]=%9.4f  %9.4f += %9.4f",
                      solpar_q, ia, vol[ia], indx_r,  et.sol_fn[indx_r],
                      gridmap[dsolvPE].energy[thread][mapindex(ix,iy)],
                      solpar_q * vol[ia] * et.sol_fn[indx_r]);
                gridmap[dsolvPE].energy[thread][mapindex(ix,iy)]
                 += solpar_q * vol[ia] * et.sol_fn[indx_r];
	        if(outlev>=LOGRUNVVV) fprintf(logFile, 
			" -> %9.4f\n", gridmap[dsolvPE].energy[thread][mapindex(ix,iy)]);
                }

#ifdef USE_BHTREE
            }/* ia loop, over all receptor atoms within distance cutoff */
#else
	    } /* ibh loop */
#endif

            /* adjust maps of hydrogen-bonding atoms by adding largest and
             * smallest interaction of all 'pair-wise' interactions with receptor atoms
             * M Pique TODO can this damage covalent maps?
             */
            if (not use_vina_potential) 
            for (int map_index = 0; map_index < num_atom_maps; map_index++) {
                if (hbondflag[map_index]) {
                    gridmap[map_index].energy[thread][mapindex(ix,iy)] 
                     += ( hbondmin[map_index] + hbondmax[map_index] );
                }
            }
	    } // end if num_atom_maps>0 or dsolvPE>=0
        } /* ix/icoord[X] loop */
    } /* iy/icoord[Y] loop */

#pragma omp ordered

if(outlev>LOGRUNV)
fprintf(logFile, "Writing iz=%d icoord=%d z=%8.2f plane/thread=%d\n", iz,icoord[Z],c[Z],thread);fflush(logFile);

            /*
             * O U T P U T . . .
             *
             * Now output one plane's grid point energies to the map files:
             *
             */
            for (int j = 0;  j < num_maps;  j++) {
		for( int iy=0; iy<n1[Y]; iy++) for(int ix=0;ix<n1[X];ix++) {
		double e = gridmap[j].energy[thread][mapindex(ix,iy)] ;
                if (!problem_wrt) {
                    if (fabs(e) < PRECISION) {
                        fprintf_retval = fprintf(gridmap[j].map_fileptr, "0.\n");
                    } else {
                        fprintf_retval = fprintf(gridmap[j].map_fileptr, "%.3f\n", (float)round3dp(e));
                    }
                    if (fprintf_retval < 0) {
                        problem_wrt = TRUE;
                    }
                }
                gridmap[j].energy_max = max(gridmap[j].energy_max, e);
                gridmap[j].energy_min = min(gridmap[j].energy_min, e);
	        }
            }
            if (floating_grid) {
                if ((!problem_wrt)&&(fprintf(floating_grid_fileptr, "%.3f\n", (float)round3dp(r_min)) < 0)) {
                    problem_wrt = TRUE;
                }
            }

    if (problem_wrt) {
        (void) sprintf( message, "Problems writing grid maps - there may not be enough disk space.\n");
        print_error( logFile, FATAL_ERROR, message );
    }
    grd_end = times( &tms_grd_end);
    ++nDone;
    timeRemaining = (float)(grd_end - grd_start) * idct * (float)(n1[Z] - nDone);
    (void) fprintf( logFile, " %6d   %8.3lf   %5.1lf%%   ", icoord[Z], cgridmin[Z] + c[Z], percentdone*(double)(iz+1));
    prHMSfixed( timeRemaining);
    (void) fprintf( logFile, "  ");
    timesys( grd_end - grd_start, &tms_grd_start, &tms_grd_end, logFile);
    (void) fflush( logFile);
} /* icoord[Z] loop */

#ifdef BOINCCOMPOUND
 boinc_fraction_done(0.9);
#endif

/*____________________________________________________________________________
 * Print a summary of extrema-values from the atomic-affinity,
 * electrostatics, and desolvation grid-maps.
 *____________________________________________________________________________*/

(void) fprintf(logFile, "\nGrid\tAtom\tMinimum   \tMaximum\n");
(void) fprintf(logFile, "Map \tType\tEnergy    \tEnergy \n");
(void) fprintf(logFile, "\t\t(kcal/mol)\t(kcal/mol)\n");
(void) fprintf(logFile, "____\t____\t_____________\t_____________\n");

for (int i = 0;  i < num_atom_maps;  i++) {
    (void) fprintf( logFile, " %d\t %s\t  %6.2lf\t%9.2le\n", i + 1, gridmap[i].type, gridmap[i].energy_min, gridmap[i].energy_max);
}

if (not use_vina_potential){
if(elecPE>=0)
(void) fprintf( logFile, " %d\t %c\t  %6.2lf\t%9.2le\tElectrostatic Potential\n", elecPE, 'e', gridmap[elecPE].energy_min, gridmap[elecPE].energy_max);

if(dsolvPE>=0)
(void) fprintf( logFile, " %d\t %c\t  %6.2lf\t%9.2le\tDesolvation Potential\n", dsolvPE, 'd', gridmap[dsolvPE].energy_min, gridmap[dsolvPE].energy_max);
(void) fprintf( logFile, "\n\n * Note:  Every pairwise-atomic interaction was clamped at %.2f\n\n", EINTCLAMP);
}
/*
 * Close all files, ************************************************************
 */

for (int i = 0;  i < num_maps;  i++) {
    (void) fclose( gridmap[i].map_fileptr);
}
if (floating_grid) {
    (void) fclose(floating_grid_fileptr);
}
/* Free up the memory allocated to the gridmap objects... */
// dont bother    free(gridmap);

(void) fprintf( logFile, "\n%s: Successful Completion.\n", programname); // do not tinker with this, used by ADT and by automated tests

job_end = times( &tms_job_end);
timesyshms( job_end - job_start, &tms_job_start, &tms_job_end, logFile);

(void) fclose( logFile);

#ifdef BOINCCOMPOUND
 boinc_fraction_done(1.);
#endif

#ifdef BOINC        
   
    boinc_finish(0);       /* should not return */
#endif

return EXIT_SUCCESS; // POSIX, defined in stdlib.h
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


static int get_map_index(const char key[]) {
    ParameterEntry * found_parm;
    found_parm = apm_find(key);
    if (found_parm != NULL)
        return found_parm->map_index;
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

#ifdef NOTNEEDED
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
#endif



/*
 * EOF
 */ 
