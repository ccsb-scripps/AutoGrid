/* autodock.h */

#include <sys/types.h>
#include "constants.h" /* autocomm.h, dpftoken.h included in constants.h */

/*----------------------------------------------------------------------------* 
 * Prototypes,                                                                * 
 *----------------------------------------------------------------------------*/

/******************************************************************************
 *      Name: autodock.h                                                      *
 *  Function: Defines prototypes for AUTODOCK.                                *
 * Copyright: (C) Garrett Matthew Morris, TSRI                                *
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: 02/28/1995                                                      *
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 02/28/95 GMM     This header added                                         *
 ******************************************************************************/


void  analysis( int Nnb, 
                char atomstuff[MAX_ATOMS][MAX_CHARS], 
                float charge[MAX_ATOMS], 
		Boole B_calcIntElec,
		float q1q2[MAX_NONBONDS],
                float clus_rms_tol, 
                float crdpdb[MAX_ATOMS][XYZ], 
                float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
                float inv_spacing, 
                float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                float econf[MAX_RUNS], 
                int irunmax, 
                float xlo, 
                float ylo, 
                float zlo, 
                int natom, 
                int nonbondlist[MAX_NONBONDS][2], 
                int nconf, 
                int ntor, 
                float qtnHist[MAX_RUNS][QUAT], 
                char smFileName[MAX_CHARS], 
                float sml_center[XYZ], 
                Boole B_symmetry_flag, 
                int torList[MAX_TORS][MAX_ATOMS], 
                float torHist[MAX_RUNS][MAX_TORS], 
                int type[MAX_ATOMS], 
                float vt[MAX_TORS][XYZ],
		char rms_ref_crds[MAX_CHARS]);

void  banner( float version_num );

void  bestpdb( int ncluster, 
               int num_in_clu[MAX_RUNS], 
               int cluster[MAX_RUNS][MAX_RUNS], 
               float econf[MAX_RUNS], 
               float crd[MAX_RUNS][MAX_ATOMS][XYZ], 
               char atomstuff[MAX_ATOMS][MAX_CHARS], 
               int natom, 
               Boole B_write_all_clusmem, 
               float clu_rms[MAX_RUNS][MAX_RUNS] );

void  check_header_float( float f1, 
			  float f2, 
			  char keyword[], 
			  char filename[] );

void  check_header_int( int i1, 
                        int i2, 
                        char axis, 
                        char filename[MAX_CHARS] );

void  check_header_line( char s1[], 
                         char s2[]
			 );

void  clmode( char atm_tyP_str[ATOM_MAPS], 
              int num_atm_maps, 
              float clus_rms_tol, 
              char hostnm[MAX_CHARS], 
              Clock jobStart, 
	      struct tms tms_jobStart, 
              Boole B_write_all_clusmem, 
              char clusFN[MAX_CHARS], 
              float crdpdb[MAX_ATOMS][XYZ], 
              float sml_center[XYZ], 
              Boole B_symmetry_flag,
	      char rms_ref_crds[MAX_CHARS]);

int  cluster_analysis( float clus_rms_tol, 
                       int cluster[MAX_RUNS][MAX_RUNS], 
                       int num_in_clus[MAX_RUNS], 
                       int isort[MAX_RUNS], 
                       int nconf, 
                       int natom, 
                       int type[MAX_ATOMS], 
                       float crd[MAX_RUNS][MAX_ATOMS][XYZ], 
		       float crdpdb[MAX_ATOMS][XYZ], 
		       float sml_center[XYZ], 
		       float clu_rms[MAX_RUNS][MAX_RUNS], 
                       Boole B_symmetry_flag,
		       float ref_crds[MAX_ATOMS][XYZ],
		       int ref_natoms);

int cmdmode( int natom,
             Clock jobStart,
             struct tms tms_jobStart,
             float xlo,
             float ylo,
             float zlo,
             float xhi,
             float yhi,
             float zhi,
             float inv_spacing,
             float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
             float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
             float WallEnergy,
             float vt[MAX_TORS][XYZ],
             int torList[MAX_TORS][MAX_ATOMS],
             int ntor,
             int Nnb,
             int nonbondlist[MAX_NONBONDS][2],
             char atomstuff[MAX_ATOMS][MAX_CHARS],
             float crdpdb[MAX_ATOMS][XYZ],
             char hostnm[MAX_CHARS],
             int type[MAX_ATOMS],
             float charge[MAX_ATOMS],
	     Boole B_calcIntElec,
	     float q1q2[MAX_NONBONDS],
	     char atm_tyP_str[ATOM_MAPS]);

void    dpffld( float *P__inv_spacing, 
                float *P__spacing, 
                char gdfldFileName[MAX_CHARS], 
                char gpfFileName[MAX_CHARS], 
                int gridpts1[XYZ], 
                int gridpts[XYZ], 
		float *xhi,
		float *yhi,
		float *zhi,
                Clock jobStart, 
                char line[LINE_LEN], 
                float *xlo, 
                float *ylo, 
                float *zlo, 
                char macromolFileName[MAX_CHARS], 
                float maP_center[XYZ], 
		struct tms tms_jobStart );

void    readmap( Boole *P_B_HaveMap, 
                int *P_Imap,
                int *P_NumAtmMaps, 
                float *P_ExtSpacing, 
                char AtmTypStr[ATOM_MAPS], 
                char ExtFldFileName[MAX_CHARS], 
                int ExtGridPts1[XYZ], 
                int ExtGridPts[XYZ], 
                Clock jobStart, 
                char line[LINE_LEN], 
                char ExtMacromolFileName[MAX_CHARS], 
                float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                float MapCenter[XYZ], 
                float MapMax[MAX_MAPS], 
                float MapMin[MAX_MAPS], 
		struct tms tmsJobStart );


void dpfmove( char  line[LINE_LEN],

	      char  atm_typ_str[ATOM_MAPS],
	      int   num_atm_maps,

	      int   *P_natom,

	      float crdpdb[MAX_ATOMS][NTRN],
	      float charge[MAX_ATOMS],
	      Boole *P_B_haveCharges,
	      int   type[MAX_ATOMS],
	      char  pdbaname[MAX_ATOMS][5],
	      char  pdbqFileName[MAX_CHARS],
	      char  atomstuff[MAX_ATOMS][MAX_CHARS],
	      int   Htype,

	      Boole *P_B_constrain,
	      int   *P_atomC1,
	      int   *P_atomC2,
	      float *P_sqlower,
	      float *P_squpper,

	      int   *P_ntor1,
	      int   *P_ntor,
	      int   tortree[MAX_TORS][MAX_ATOMS],
	      float vt[MAX_TORS][NTRN],

	      int   *P_Nnb,
	      int   Nnbonds[MAX_ATOMS],
	      int   nonbondlist[MAX_NONBONDS][2],

	      Clock jobStart,
	      struct tms tms_jobStart,
	      char  hostnm[MAX_CHARS]);

void    dpftypes( int *P__Htype, 
                  int *P__num_all_maps,
                  int *P__num_atm_maps, 
                  char atm_tyP_str[ATOM_MAPS], 
                  char line[LINE_LEN] );

float  eintcal( int nonbondlist[MAX_NONBONDS][2], 
                float eint_table[NEINT][ATOM_MAPS][ATOM_MAPS], 
                float tcoord[MAX_ATOMS][XYZ], 
                int atmtyp[MAX_ATOMS], 
                int Nnb,
		Boole B_calcIntElec,
		float q1q2[MAX_NONBONDS]);

float  eintcalPrint( int nonbondlist[MAX_NONBONDS][2], 
                float eint_table[NEINT][ATOM_MAPS][ATOM_MAPS], 
                float tcoord[MAX_ATOMS][XYZ], 
                int atmtyp[MAX_ATOMS], 
                int Nnb,
		Boole B_calcIntElec,
		float q1q2[MAX_NONBONDS]);

int getpdbcrds( char rms_ref_crds_FN[MAX_CHARS],
		float ref_crds[MAX_ATOMS][XYZ] );

float  getrms( float Crd[MAX_ATOMS][XYZ], 
               float CrdRef[MAX_ATOMS][XYZ], 
               Boole B_symmetry_flag, 
               int natom, 
               int type[MAX_ATOMS] );

void initautodock(float *P_e0, 
                  int Nnb, 
                  char atomstuff[MAX_ATOMS][MAX_CHARS], 
                  float charge[MAX_ATOMS], 
		  Boole B_calcIntElec,
		  float q1q2[MAX_NONBONDS],
                  float crd[MAX_ATOMS][XYZ], 
                  float crdpdb[MAX_ATOMS][XYZ], 
                  float elec[MAX_ATOMS], 
                  float emap[MAX_ATOMS], 
                  float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
                  float xhi, 
                  float yhi, 
                  float zhi, 
                  float inv_spacing, 
                  float xlo, 
                  float ylo, 
                  float zlo, 
                  float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                  int natom, 
                  int nonbondlist[MAX_NONBONDS][2], 
                  int ntor, 
                  float qtn0[QUAT], 
                  int torList[MAX_TORS][MAX_ATOMS], 
                  float tor0[MAX_TORS], 
                  int type[MAX_ATOMS], 
                  float vt[MAX_TORS][XYZ],
		  unsigned short US_torProfile[MAX_TORS][NTORDIVS],
		  Boole B_isTorConstrained[MAX_TORS]);

int  input_state( FILE *fp, 
                  char line[LINE_LEN], 
                  float qm[QUAT], 
                  float tor[MAX_TORS], 
                  int ntor, 
		  int *P_istep, 
                  float *P_energy, 
		  float *P_eint, 
                  char *P_lastmove );

void intnbtable(Boole *P_B_havenbp,
                int *P_a1,
                int *P_a2,
		int num_atm_maps,
                char atm_tyP_str[ATOM_MAPS],
                float cA,
                float cB,
                int xA,
                int xB,
                float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS] );

int  main( int  argc, 
           char **argv );

void  mkNewState( float qtn[NQTN],
                float tor[MAX_TORS],
                float qtnLast[NQTN],
                float torLast[MAX_TORS],
                float qtnChange[NQTN],
                float torChange[MAX_TORS],
                float vt[MAX_TORS][NTRN],
                int torList[MAX_TORS][MAX_ATOMS],
                int ntor,
                float crd[MAX_ATOMS][NTRN],
                float crdpdb[MAX_ATOMS][NTRN],
                int natom,
                float trnStep,
                float qtwStep,
                float torStep );
                /* float trnStepHi,
                ** float qtwStepHi,
                ** float torStepHi,
		** float eLast,
		** float xloTrn,
		** float xhiTrn,
		** float yloTrn,
		** float yhiTrn,
		** float zloTrn,
		** float zhiTrn);
		*/

void  mkTorTree( int atomnumber[MAX_RECORDS],
                char record[MAX_RECORDS][LINE_LEN],
                int nrecord,
                int torList[MAX_TORS][MAX_ATOMS],
                int *P_ntor,
                char smFileName[MAX_CHARS],
                char pdbaname[MAX_ATOMS][5],
                Boole *P_B_constrain,
                int *P_atomC1,
                int *P_atomC2,
                float *P_sqlower,
                float *P_squpper );

void  nbe( char atm_tyP_str[ATOM_MAPS],
	   float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
           int num_atm_maps );

void  nonbonds( float crdpdb[MAX_ATOMS][XYZ],
                int nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
                int natom,
                int atomnumber[MAX_RECORDS], 
                int nrecord,
                char record[MAX_RECORDS][LINE_LEN],
                int piece[MAX_ATOMS],
                int Htype,
                int type[MAX_ATOMS] );

int  openfile( char filename[MAX_CHARS],
               char mode[],
               FILE **fp );

void  output_state( FILE *fp,
                    float qm[QUAT],
                    float tor[MAX_TORS],
                    int ntor,
                    int istep,
                    float energy,
                    float eint,
                    char lastmove,
                    Boole B_watch,
                    char FN_watch[MAX_CHARS],
                    char atomstuff[MAX_ATOMS][MAX_CHARS],
                    int natom,
                    float crd[MAX_ATOMS][XYZ]);

int  parse_com_line( char line[LINE_LEN] );

int  parse_dpf_line( char line[LINE_LEN] );

int  parse_pdbq_line( char line[LINE_LEN] );

int  parse_trj_line( char line[LINE_LEN] );

void  print_2x( FILE *stream1,
                FILE *stream2,
                char *string );

void  print_atomic_energies( int natom,
                            char atomstuff[MAX_ATOMS][MAX_CHARS],
                            int type[MAX_ATOMS],
                            float emap[MAX_ATOMS],
                            float elec[MAX_ATOMS],
                            float coord[MAX_ATOMS][XYZ],
                            float charge[MAX_ATOMS] );

void  print_avsfld( FILE *logFile,
                   int veclen,
                   int natom,
                   int nframe,
                   int offset[VECLENMAX],
                   int stride,
                   char label[MAX_CHARS],
                   char filename[MAX_CHARS] );

void  printdate( FILE *fp, int flag );

void  print_docked( int irun,
                    char smFileName[MAX_CHARS],
		    char dpfFN[MAX_CHARS],
                    float sml_center[XYZ],
                    float qtnold[QUAT],
                    int ntor,
                    float torold[MAX_TORS],
                    float eintra,
                    float einter,
                    int natom,
                    char atomstuff[MAX_ATOMS][MAX_CHARS],
                    float crd[MAX_ATOMS][XYZ],
                    float emap[MAX_ATOMS],
                    float elec[MAX_ATOMS],
                    float charge[MAX_ATOMS] );

void  printhms( float t );

void  prClusterHist( int ncluster,
                    int irunmax,
                    float clus_rms_tol,
                    int num_in_clu[MAX_RUNS],
                    int cluster[MAX_RUNS][MAX_RUNS],
                    float econf[MAX_RUNS],
                    float clu_rms[MAX_RUNS][MAX_RUNS] );

void  print_rem(  FILE *outFile,
                  int Rank,
                  int NumMem,
                  int Run,
                  float Energy,
                  float rms );

void  print_time_since( time_t then );

void  qmultiply( float *P_qm, 
                 float qml[QUAT], 
                 float qmr[QUAT]);

void  qtransform( float q[QUAT],
                  float tcoord[MAX_ATOMS][XYZ],
                  int natom );

void  quicksort( float e[], 
                 int isort[],
                 int left,
                 int right );

void  readpdbq( char line[LINE_LEN],
                float crd[XYZ], 
                float *P_q );

void  set_cmd_io_std( void );

int  setflags( int argc, 
               char **argv );

void simanneal( int *P__nconf, 
                int Nnb, 
                float WallEnergy, 
                char atomstuff[MAX_ATOMS][MAX_CHARS], 
                float charge[MAX_ATOMS], 
		Boole B_calcIntElec,
		float q1q2[MAX_NONBONDS],
                float crd[MAX_ATOMS][XYZ], 
                float crdpdb[MAX_ATOMS][XYZ], 
		char dpfFN[MAX_CHARS],
                float e0, 
                float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
                float econf[MAX_RUNS], 
                Boole B_either, 
                float elec[MAX_ATOMS], 
                float emap[MAX_ATOMS], 
                float xhi, 
                float yhi, 
                float zhi, 
                int icyclemax, 
                float inv_spacing, 
                int irunmax, 
                Clock jobStart, 
                float xlo, 
                float ylo, 
                float zlo, 
                float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                int naccmax, 
                int natom, 
                int nonbondlist[MAX_NONBONDS][2], 
                int nrejmax, 
                int ntor1, 
                int ntor, 
                int outlev, 
                float qtn0[QUAT], 
                float qtnHist[MAX_RUNS][QUAT], 
                float qtwFac, 
                Boole B_qtwReduc, 
                float qtwStep0, 
                float qtwStep, 
                Boole B_selectmin, 
                char smFileName[MAX_CHARS], 
                float sml_center[XYZ], 
                float RT0, 
                Boole B_RTChange, 
                float RTFac, 
                struct tms tms_jobStart, 
                int torList[MAX_TORS][MAX_ATOMS], 
                float tor0[MAX_TORS], 
                float torFac, 
                float torHist[MAX_RUNS][MAX_TORS], 
                Boole B_torReduc, 
                float torStep0, 
                float torStep, 
                char trjFileName[MAX_CHARS], 
                int trj_cyc_max, 
                int trj_cyc_min, 
                int trj_freq, 
                float trnFac, 
                Boole B_trnReduc, 
                float trnStep0, 
                float trnStep, 
                int type[MAX_ATOMS], 
                float vt[MAX_TORS][XYZ], 
                Boole B_write_trj, 
                Boole B_constrain, 
                int atomC1, 
                int atomC2, 
                float sqlower, 
                float squpper,
		Boole B_linear_schedule,
		float RTreduc,
		float maxrad,
		Boole B_watch,
		char FN_watch[MAX_CHARS],
		unsigned short US_torProfile[MAX_TORS][NTORDIVS],
		Boole B_isTorConstrained[MAX_TORS] );

void sort_enrg( float econf[MAX_RUNS],
                int isort[MAX_RUNS],
                int nconf );

void  stop( char reason[LINE_LEN] );

int  strindex( char s[], 
               char t[] );

void  success( char hostnm[MAX_CHARS],
               Clock jobStart,
               struct tms tms_jobStart );

void  summarizegrids( char atm_tyP_str[ATOM_MAPS],
                      float mapmax[MAX_MAPS],
                      float mapmin[MAX_MAPS],
                      int num_all_maps,
                      int num_atm_maps );

void swap( int v[], int i, int j );

time_t  timenow( void );

void  timesys( Clock duration,
               struct tms *start,
               struct tms *end );

void  timesyshms( Clock duration,
                  struct tms *start,
                  struct tms *end );

void  torNorVec( float crdpdb[MAX_ATOMS][XYZ],
                 int ntor,
                 int torList[MAX_TORS][MAX_ATOMS],
                 float vt[MAX_TORS][XYZ] );

void  torsion( float v[MAX_TORS][XYZ], 
               float tor[MAX_TORS], 
               float crd[MAX_ATOMS][XYZ], 
               int torlis[MAX_TORS][MAX_ATOMS], 
               int ntor );

float  trilinterp( float tcoord[MAX_ATOMS][XYZ], 
                   float charge[MAX_ATOMS], 
                   int type[MAX_ATOMS], 
                   int natom, 
                   float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                   float inv_spacing, 
                   float elec[MAX_ATOMS], 
                   float evdW[MAX_ATOMS], 
                   float xlo, 
                   float ylo, 
                   float zlo );
		   /* float lo[XYZ] ); */

int  usage( void );

void  warn_bad_file( char filename[MAX_CHARS],
                     char message[LINE_LEN] );

void  weedbonds( int natom,
                 char pdbaname[MAX_ATOMS][5],
                 int piece[MAX_ATOMS],
                 int ntor,
                 int torList[MAX_TORS][MAX_ATOMS],
                 int *P__Nnb,
                 int Nnbonds[MAX_ATOMS],
                 int nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
                 int nonbondlist[MAX_NONBONDS][2] );


/*----------------------------------------------------------------------------* 
 * End of file                                                                * 
 *----------------------------------------------------------------------------*/
