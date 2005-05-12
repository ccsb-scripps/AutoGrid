/* structs.h */
#ifndef _PARAMETERS_H
#define _PARAMETERS_H

//#include "constants.h"
//#include "typedefs.h"

/* *****************************************************************************
 *      Name: parameters.h                                                       *
 *  Function: Defines structures used in Molecular Applications.              *
 * Copyright: (C) Garrett Matthew Morris, TSRI                                *
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: SEP/07/1995                                                     *
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 02/28/95 GMM     This header added                                         *
 ***************************************************************************** */

/* ____________________________________________________________________________ */

/* ______________________________________________________________________________
** Parameter Dictionary */

#include <search.h>

#define MAX_NUM_AUTOGRID_TYPES 100
#define MAX_LEN_AUTOGRID_TYPE 7

enum hbond_type
{ NON, DS, D1, AS, A1, A2 };	/* hbonding character: */

typedef struct parameter_entry
{				// was "parm_info" in earlier AutoGrid 4 code
  char autogrid_type[MAX_LEN_AUTOGRID_TYPE + 1];	/* autogrid_type is based on babel_types assigned by PyBabel */
  double Rij;			/* Lennard-Jones equilibrium separation */
  double epsij;			/* Lennard-Jones energy well-depth */
  double vol;			/* solvation volume */
  double solpar;		/* solvation parameter */
  hbond_type hbond;		/* hbonding character: 
				   NON: none, 
				   DS: spherical donor 
				   D1: directional donor
				   AS: spherical acceptor
				   A1: acceptor of 1 directional hbond
				   A2: acceptor of 2 directional hbonds */
  double Rij_hb;		/* 12-10 Lennard-Jones equilibrium separation */
  double epsij_hb;		/* 12-10 Lennard-Jones energy well-depth */
  int rec_index;		/* used to set up receptor atom_types */
  int map_index;		/* used to set up map atom_types */
  int bond_index;		/* used to set up bonds; corresponds to the enum in mdist.h */
} ParameterEntry;

typedef struct linear_FE_model
{
    double coeff_vdW;                 // Free energy coefficient for van der Waals term
    double coeff_hbond;               // Free energy coefficient for H-bonding term
    double coeff_estat;               // Free energy coefficient for electrostatics term
    double coeff_desolv;              // Free energy coefficient for desolvation term
    double coeff_tors;                // Free energy coefficient for torsional term

    double stderr_vdW;                // Free energy standard error for van der Waals term
    double stderr_hbond;              // Free energy standard error for H-bonding term
    double stderr_estat;              // Free energy standard error for electrostatics term
    double stderr_desolv;             // Free energy standard error for desolvation term
    double stderr_tors;               // Free energy standard error for torsional term
} Linear_FE_Model;

#endif
/* EOF */
