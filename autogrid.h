/* autogrid.h */

#include "autocomm.h"
#include "gpftoken.h"

/******************************************************************************/
/*      Name: autogrid.h                                                      */
/*  Function: Header file for Autogrid.                                       */
/* Copyright: (C) 1995, TSRI                                                  */
/*----------------------------------------------------------------------------*/
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
/*      Date: 02/06/95  6-FEB-1995                                            */
/*----------------------------------------------------------------------------*/
/*    Inputs: None.                                                           */
/*   Returns: Parameters, Macro substitutions, Prototyped functions.          */
/*   Globals: (see 'autoglobal.h')                                            */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 04/01/93 GMM     Created for use in makefile.                              */
/******************************************************************************/

#define MAX_DIST     16384   /* Maximum distance in 100ths of an Angstrom.    */
                             /*  = 163.84 Angstroms                           */
#define MAX_ATOMS    32768   /* Maximum number of atoms in macromolecule.     */
/*    32768 = 2^15	*/
/*    int 16-bit two's complement ranges 0-32767, 0 to (2^15 - 1)	*/

#define ORDERED 	0
#define CYLINDRICAL 	1
#define SPHERICAL 	2

#define A_DIVISOR    100.    /* Angstrom is divided by this in look-up table. */

#define	NBCUTOFF     8.      /* non-bond cutoff = 8 Angstroms.                */

#define PRECISION 0.0001 /* fabs(Energies) less than this will be written as '0.' */

/*----------------------------------------------------------------------------*/
/* Macros,                                                                    */
/*----------------------------------------------------------------------------*/

#define sq(a)               ( (a) * (a) )
#define hypotenuse(x,y,z)   ( sqrt( (x)*(x) + (y)*(y) + (z)*(z) )  )
#define sq_hyp(x,y,z)       ( (x)*(x) + (y)*(y) + (z)*(z) )
#define max(x,y)            ( ((x) > (y)) ? (x) : (y) )
#define min(x,y)            ( ((x) < (y)) ? (x) : (y) )
#define angstrom(i)         ( ( (double) (i) ) / A_DIVISOR )
#define lookup(r)           ( (int) ( (r) * A_DIVISOR ) )
#define equal(a,b,n)        ( strncmp(a,b,(n)) == 0 )

/*----------------------------------------------------------------------------*/
/* Prototypes,                                                                */
/*----------------------------------------------------------------------------*/

#include "prototypes.h"


#define MAX_NUM_AUTOGRID_TYPES 100
#define NUM_ALL_TYPES 20   /*??? IS THIS REASONABLE???*/
#define MAX_LEN_AUTOGRID_TYPE 7
/*enum babel_type { C3, C2, C1, Cac, Cpl, 
                  N3pl, Nox, N3, Ntr, Npl, N1, Nam, 
                  O3, O2, Om, 
                  S3pl, S3, S2, Sac, Sox, S,
                  HC, H, 
                  P3, Pac, Pox, 
                  B, Bac, Box, 
                  Al, As , Be, Br, Ca, Cl, Cu, Fl, 
                  Fe, Ge, I, K, Li, Mg, Mn, Na, Ni, 
                  Pb, Si, Zn};*/



struct parm_info {
        char  autogrid_type[MAX_LEN_AUTOGRID_TYPE + 1]; /*KEY:  autogrid_type
                                                    is based on babel_types
                                                    assigned by PyBabel*/
        int   num;              /* ct of autogrid_type as above*/
        int   map_number;        /*index of autogrid_type in ligand_types*/
        int   rec_index;        /*index of autogrid_type in receptor_types*/
        double Rij;             /*Lennard-Jones equilibrium separation*/
        double epsij;           /*Lennard-Jones energy well-depth*/
        double vol;             /*solvation volume*/
        double solpar;          /*solvation parameter*/
        double constant;       /*OLDSTYLE: constant parameter*/
        /*lists of parameters for this type vs each receptor type*/
        double rij[NUM_ALL_TYPES];       /*rij per thistype-receptor_type*/
        double eps[NUM_ALL_TYPES];       /*epsij per thistype-receptor_type*/
        double xA[NUM_ALL_TYPES];       /*exponent1 per thistype-receptor_type*/
        double xB[NUM_ALL_TYPES];       /*exponent2 per thistype-receptor_type*/
        int hb[NUM_ALL_TYPES];       /*hb per thistype-receptor_type*/
        int is_metal;                /*1 if FE,ZN,MN,MG,CA*/
};


/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
