#include <sys/types.h>
#include <string.h>
#include <stdio.h>

/******************************************************************************/
/*      Name: autoglobal.h                                                    */
/*  Function: Global variables for Autodock modules.                          */
/* Copyright: (C) TSRI                                                        */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett M. Morris                                               */
/*                                                                            */
/*            e-mail: garrett@scripps.edu				      */
/*                                                                            */
/*            The Scripps Research Institute                                  */
/*            Department of Molecular Biology, MB5                            */
/*            10666 North Torrey Pines Road                                   */
/*            La Jolla, CA 92037.                                             */
/*                                                                            */
/*      Date: 03/18/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: None.                                                           */
/*   Returns: Nothing.                                                        */
/*   Globals: programname, AutoGridHelp, AutoDockHelp, command_mode,          */
/*            command_in_fp, command_out_fp, GPF, logFile.                    */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 03/18/93 GMM     Created.                                                  */
/******************************************************************************/


/*----------------------------------------------------------------------------*/
/* Globals,                                                                   */
/*----------------------------------------------------------------------------*/

char    *programname;
char    AutoDockHelp[] = "           -p parameter_filename\n           -l log_filename\n           -o (Use old PDBQ format, charge q in columns 55-61)\n           -k (Keep original residue numbers)\n           -i (Ignore header-checking)\n           -t (Parse the PDBQ file to check torsions, then stop.)\n           -c < command_file (Command mode, by file)\n           -c | control_program (Command mode, by control_program)\n\n";

char    AutoGridHelp[] = "-p parameter_filename\n-l log_filename\n-o (old PDBQ format)\n-d (increment debug level)\n-u (display this message)\n";

char    dock_param_fn[MAX_CHARS];
char    GPF_filename[MAX_CHARS];

int     command_mode = FALSE;
int     debug = 0;
int	ElecMap = 0;
int     ignore_errors = FALSE;
int     keepresnum = 0;
int     oldpdbq = FALSE;
int     parse_tors_mode = FALSE;

float	idct = 1.;

FILE    *command_in_fp;
FILE    *command_out_fp;
FILE    *parFile;
FILE    *GPF_fileptr;
FILE    *logFile;

/*
** struct  Quat {
**             float angle;
**             float vec[XYZ];
**             };
*/

/* EOF */
