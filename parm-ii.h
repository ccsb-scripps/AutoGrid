#ifndef PARM_II
#define PARM_II

/* gpf3gen.awk parameters */

/*
 *  AutoDock Binding Free Energy Model 140n
 *  New well-depths using new Solvation Model,
 *  multiplied by van der Waals coefficient, 0.1485
 */

/*
 *  Well-depths prior to Solvation Model, multiplied by FE_vdW_coeff
 */

#define FE_vdW_coeff    0.1485
#define FE_estat_coeff  0.1146
#define FE_hbond_coeff  0.0656
#define FE_tors_coeff   0.3113
#define FE_desol_coeff  0.1711


/*  C-C,  A-A,  N-N,  O-O,  S-S,  P-P,  H-H,  Fe-Fe
 *  F-F,  Cl-Cl,  Br-Br,  I-I,  Mg-Mg,  Zn-Zn,  Ca-Ca,  Xx-Xx */

/*
C    A    N    O    S    P    H    f
F    c    b    I    M    Z    L    X
*/

/*  Recognizable atom types are C,A,N,O,S,P & H
 *  and...                      f,F,c,b,I (Fe, F, Cl, Br, I)
 *  and...                      X,M (Any, Mg)
 *  and...                      Z,L (Zn, Ca)
 *  (note: A is aromatic carbon) */

#define ALL_TYPES "CANOSPHfFcbIMZLX"
/*                 0123456789012345 */
#define NUM_ALL_TYPES 16                 /* must be the length of ALL_TYPES */

#define CARBON		0 /*  C  */
#define AROMATIC_CARBON 1 /* A */
#define NITROGEN	2 /*  N  */
#define OXYGEN		3 /*  O  */
#define SULPHUR		4 /*  S  */
#define PHOSPHORUS  5 /*  P  */
#define HYDROGEN	6 /*  H  */
#define IRON        7 /*  Fe  */
#define FLUORINE    8 /*  F  */
#define CHLORINE    9 /*  Cl  */
#define BROMINE     10 /* Br  */
#define IODINE      11 /* I  */
#define MAGNESIUM	12 /*  Mg  */
#define METAL		MAGNESIUM /*  M  */
#define ZINC		13 /*  Z  */
#define CALCIUM		14 /*  L  */
#define UNKNOWN		15 /*  X  */

#define COVALENT    16 /*  COVALENTTYPE  */
#define COVALENTTYPE  'z'
#define COVALENT2   17 /*  COVALENTTYPE2  */
#define COVALENTTYPE2 'y'
#define NATOMTYPES	7 /* Number of atom types for atomic interactions  */
#define MAX_TYPES   8 /* Maximum number of atom types used. */

/* Prototypes */
float get_Rij(int type1, int type2);
float get_epsij(int type1, int type2);
float get_Rij_Hbond(int type1, int type2);
float get_epsij_Hbond(int type1, int type2);
float get_SolVol(int type1);
float get_SolPar(int type1);
float get_SolCon(int type1);

#endif
/* EOF */
