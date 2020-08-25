/*

 $Id: gpfparser.cpp,v 1.18 2020/08/25 20:26:43 mp Exp $

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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "gpftoken.h"
#include "autogrid.h"
#include "constants.h"

int gpfparser( char line[LINE_LEN] )

/******************************************************************************/
/*      Name: gpfparser                                                       */
/*  Function: Parse the AutoGrid parameter file line                          */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 02/01/95 (1-feb-1995)                                           */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 02/01/95 GMM     Entered code.                                             */
/******************************************************************************/

{
#define streq(a,b) (0==strcasecmp(a,b)) // case-independent string match
    int token = -1 ;	       /* return -1 if nothing is recognized. */
    char c[LINE_LEN];  // keyword token (first white-space-terminated string)
    char * tptr = &line[0]; // start of keyword token in "line"

    /* skip any leading whitespace */
    while( isascii(*tptr) && isspace(*tptr) ) tptr++;

    /* copy non-whitespace chars from "line" to "c" */
    int i=0;
    for ( i=0; isascii(tptr[i]) && !isspace(tptr[i]); i++ ) c[i] = tptr[i];

    c[i] = '\0'; // NULL-terminate c

    if ((c[0]=='\n')||(c[0]=='\0')) {
        token = GPF_NULL;

    } else if (c[0]=='#') {
        token = GPF_COMMENT;

    } else if (streq(c,"receptor_types")) {
        token = GPF_RECEPTOR_TYPES;

    } else if (streq(c,"receptor")) {
        token = GPF_RECEPTOR;

    } else if (streq(c,"gridfld")) {
        token = GPF_GRIDFLD;

    } else if (streq(c,"npts")) {
        token = GPF_NPTS;

    } else if (streq(c,"spacing")) {
        token = GPF_SPACING;

    } else if (streq(c,"gridcenter")) {
        token = GPF_GRIDCENTER;

    } else if (streq(c,"types")) {
        token = GPF_LIGAND_TYPES;

    } else if (streq(c,"ligand_types")) {
        token = GPF_LIGAND_TYPES;


    } else if (streq(c,"map")) {
        token = GPF_MAP;

    } else if (streq(c,"elecmap")) {
        token = GPF_ELECMAP;

    } else if (streq(c,"dsolvmap") || streq(c,"desolvmap")) {
        token = GPF_DSOLVMAP;

    } else if (streq(c,"covalentmap")) {
        token = GPF_COVALENTMAP;

    } else if (streq(c,"coeff_vdw")) {
        token = GPF_COEFF_VDW;

    } else if (streq(c,"coeff_hbond")) {
        token = GPF_COEFF_HBOND;

    } else if (streq(c,"coeff_estat")) {
        token = GPF_COEFF_ESTAT;

    } else if (streq(c,"coeff_desolv")) {
        token = GPF_COEFF_DESOLV;

    } else if (streq(c,"coeff_tors")) {
        token = GPF_COEFF_TORS;

    } else if (streq(c,"nbp_coeffs")) {
        token = GPF_NBP_COEFFS;

    } else if (streq(c,"nbp_r_eps")) {
        token = GPF_NBP_R_EPS;

    } else if (streq(c,"dielectric")) {
        token = GPF_DIEL;

    } else if (streq(c,"qasp")) {
        token = GPF_QASP;

    } else if (streq(c,"fmap")) {
        token = GPF_FMAP;

    } else if (streq(c,"disorder_h")) {
        token = GPF_DISORDER;

    } else if (streq(c,"smooth")) {
        token = GPF_SMOOTH;

    } else if (streq(c,"sol_par")) {
        token = GPF_SOL_PAR;

    } else if (streq(c,"constant")) {
        token = GPF_CONSTANT;

    } else if (streq(c,"parameter_file")) {
        token = GPF_PARAM_FILE;

    } else if (streq(c,"use_vina_potential")) {
        token = GPF_USE_VINA_POTENTIAL;

    } else if (streq(c,"outlev")) {
        token = GPF_OUTLEV;

    } else if (streq(c,"map_receptor_interior")) {
        token = GPF_MAP_RECEPTOR_INTERIOR;

    } else if (streq(c,"cmap") || streq(c,"constriction_map")) {
        token = GPF_CMAP;

    } else if (streq(c,"constriction_distance_cutoff")) {
        token = GPF_CONSTRICTION_DISTANCE_CUTOFF;

    } else if (streq(c,"separate_desolvation_maps")) {
        token = GPF_SEPARATE_DESOLVATION_MAPS;
    }
    return(token);
}
/* EOF */
