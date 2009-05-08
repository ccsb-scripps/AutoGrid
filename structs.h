/*

 $Id: structs.h,v 1.7 2009/05/08 23:17:35 rhuey Exp $

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

#ifndef _STRUCTS_H
#define _STRUCTS_H

#include "constants.h"
#include "typedefs.h"

/* *****************************************************************************
 *      Name: structs.h                                                       *
 *  Function: Defines structures used in Molecular Applications.              *
 *Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
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

#ifdef USE_INT_AS_LONG
typedef int FourByteLong;
typedef unsigned int UnsignedFourByteLong;
#else
typedef long FourByteLong;
typedef unsigned long UnsignedFourByteLong;
#endif

/* ____________________________________________________________________________ */

typedef double Vector [3];  // for vectors and points

/* ____________________________________________________________________________ */

typedef struct coord
{
  double x;			/* Cartesian x-coordinate */
  double y;			/* Cartesian y-coordinate */
  double z;			/* Cartesian z-coordinate */
} Coord;

/* ____________________________________________________________________________ */

typedef struct quat
{
  double nx;			/* unit vector's x-component */
  double ny;			/* unit vector's y-component */
  double nz;			/* unit vector's z-component */
  double ang;			/* angle of rotation about unit-vector */
  double w;			/* quaternion's w-component */
  double x;			/* quaternion's x-component */
  double y;			/* quaternion's y-component */
  double z;			/* quaternion's z-component */
  double qmag;			/* quaternion's 4-D magnitude */
} Quat;

/* ____________________________________________________________________________ */

typedef struct energy
{
  double total;			/* total energy */
  double intra;			/* intramolecular energy, a.k.a. "internal" energy */
  double inter;			/* intermolecular energy */
  double FE;			/* estimated Free Energy of binding */
} Energy;

/* ____________________________________________________________________________ */

typedef struct state
{
  Coord T;			/* coordinates of center of molecule */
  Quat Q;			/* rigid-body orientation */
  double tor[MAX_TORS];		/* torsion angles in radians */
  int ntor;			/* number of torsions in molecule */
  int hasEnergy;		/* if 0, this state has an undefined energy */
  Energy e;			/* energy structure */
} State;

/* ____________________________________________________________________________ */

typedef struct molecule
{
  Real     crdpdb[MAX_ATOMS][SPACE];	    /* original coordinates of atoms */
  Real     crd[MAX_ATOMS][SPACE];      	/* current coordinates of atoms */
  char              atomstr[MAX_ATOMS][MAX_CHARS];	/* strings describing atoms, from PDB file, cols,1-30. */
  int               natom;			                /* number of atoms in molecule */
  Real     vt[MAX_TORS][SPACE];        	/* vectors  of torsions */
  int               tlist[MAX_TORS][MAX_ATOMS];	    /* torsion list of movable atoms */
  State             S;		                    	/* state of molecule */
} Molecule;

/* ____________________________________________________________________________ */

typedef struct rotamer
{
  double tor[MAX_TORS_IN_ROTAMER];	/* torsion angles in radians */
  int ntor;			/* number of torsions */
} Rotamer;

/* ____________________________________________________________________________ */

typedef struct chargestruct
{
    double charge;
    double abs_charge;
    double qsp_abs_charge;
} Charge;

/* ____________________________________________________________________________ */

typedef struct atom
{
  double    coords[3];			    /* transformed point */
  Boole     has_charge;			    /* TRUE if the atom has a charge */

  double    charge;
  double    abs_charge;
  double    qsp_abs_charge;

  int       type;			        /* atom type as integer */
  char      type_string[MAX_CHARS]; /* atom type as string */
  Boole     is_hydrogen;		    /* TRUE if atom is a hydrogen */

  int       serial;			        /* serial ID */
  char      name[5];			    /* PDB atom name; formerly "pdbaname" */
  char      stuff[MAX_CHARS];       /* PDB atom string; formerly "atomstuff" */

  int       nnb;			        /* number of non-bonds for this atom */
} Atom;
/* ____________________________________________________________________________ */

typedef struct bond
{
  Atom *atom1;
  Atom *atom2;
  double bondLength;
  Coord bondVector;
} Bond;
/* ____________________________________________________________________________ */

typedef struct pair_id
{
  Atom *atom1;			/* pointer to one atom in pair */
  Atom *atom2;			/* pointer to other atom */
} PairID;
/* ____________________________________________________________________________ */

typedef struct dist_constraint
{
  PairID bond;			/* two atoms defining distance constraint */
  double length;		/* current bond length */
  double lower;			/* lower bound on distance */
  double upper;			/* upper bound on distance */
} DisCon;

/* ____________________________________________________________________________ */

typedef struct torsion
{
  PairID rotbnd;		/* atom serial-IDs of rotatable bond */
  int nmoved;			/* number of atoms moved by this */
  int IDmove[MAX_ATOMS];	/* atom serial-IDs of atoms moved by this */
  Coord vt;			/* bond-vector of rotatable bond */
} Torsion;

/* ______________________________________________________________________________
** Molecular fragments, side-chains, even entire ligands */

typedef struct group
{
  int natom;			/* Number of atoms in fragment */
  Atom atm[MAX_ATOMS];		/* Atom data */
  int ntor;			/* Number of torsions in fragment */
  Torsion tors[MAX_TORS];	/* Torsion data */
  int nnb;			/* Number of non-bonds in fragment */
  PairID nbs[MAX_NONBONDS];	/* Non-bond data */
  Boole B_constrain;		/* TRUE if any distance constraints */
  DisCon distcon;		/* Distance constraint data */
  char pdbqfilnam[MAX_CHARS];	/* PDBQ filename holding these data */
} Group;

#include "grid.h"

/* ______________________________________________________________________________
** Parameter Dictionary */

#include "parameters.h"
/* ______________________________________________________________________________ */

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

/* ______________________________________________________________________________ */
/* Energy Lookup Tables */

typedef struct energy_tables
{
    Real e_vdW_Hb[NEINT][MAX_ATOM_TYPES][MAX_ATOM_TYPES];  // vdW & Hb energies
    Real sol_fn[NEINT];                            // distance-dependent desolvation function
    Real epsilon_fn[NDIEL];                        // distance-dependent dielectric function
    Real r_epsilon_fn[NDIEL];                      // r * distance-dependent dielectric function
} EnergyTables;


#endif
/* EOF */
