/* structs.h */
#ifndef _PARAMETERS_H
#define _PARAMETERS_H

/* ______________________________________________________________________________
** Parameter Dictionary */

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

#endif
/* EOF */
