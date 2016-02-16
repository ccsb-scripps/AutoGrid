/*

 $Id: gpfparser.cpp,v 1.14 2016/02/16 23:49:27 mp Exp $

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
    int l, i, token = -1 ;	       /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    l = (int)strindex(line, " ");
    if (l == -1) {
        l = (int)strindex(line, "\t");
        if (l == -1) {
            l = (int)strlen(line);
	}
    }
    for(i=0; i<l; i++) {
        c[i] = (char)tolower( (int)line[i] );
    }

    if ((c[0]=='\n')||(c[0]=='\0')) {
        token = GPF_NULL;

    } else if (c[0]=='#') {
        token = GPF_COMMENT;

    } else if (equal(c,"receptor_types")) {
        token = GPF_RECEPTOR_TYPES;

    } else if (equal(c,"receptor")) {
        token = GPF_RECEPTOR;

    } else if (equal(c,"gridfld")) {
        token = GPF_GRIDFLD;

    } else if (equal(c,"npts")) {
        token = GPF_NPTS;

    } else if (equal(c,"spacing")) {
        token = GPF_SPACING;

    } else if (equal(c,"gridcenter")) {
        token = GPF_GRIDCENTER;

    } else if (equal(c,"types")) {
        token = GPF_LIGAND_TYPES;

    } else if (equal(c,"ligand_types")) {
        token = GPF_LIGAND_TYPES;


    } else if (equal(c,"map")) {
        token = GPF_MAP;

    } else if (equal(c,"elecmap")) {
        token = GPF_ELECMAP;

    } else if (equal(c,"dsolvmap") || equal(c,"desolvmap")) {
        token = GPF_DSOLVMAP;

    } else if (equal(c,"covalentmap")) {
        token = GPF_COVALENTMAP;

    } else if (equal(c,"nbp_coeffs")) {
        token = GPF_NBP_COEFFS;

    } else if (equal(c,"nbp_r_eps")) {
        token = GPF_NBP_R_EPS;

    } else if (equal(c,"dielectric")) {
        token = GPF_DIEL;

    } else if (equal(c,"qasp")) {
        token = GPF_QASP;

    } else if (equal(c,"fmap")) {
        token = GPF_FMAP;

    } else if (equal(c,"disorder_h")) {
        token = GPF_DISORDER;

    } else if (equal(c,"smooth")) {
        token = GPF_SMOOTH;

    } else if (equal(c,"sol_par")) {
        token = GPF_SOL_PAR;

    } else if (equal(c,"constant")) {
        token = GPF_CONSTANT;

    } else if (equal(c,"parameter_file")) {
        token = GPF_PARAM_FILE;

    } else if (equal(c,"use_vina_potential")) {
        token = GPF_USE_VINA_POTENTIAL;

    } else if (equal(c,"outlev")) {
        token = GPF_OUTLEV;

    } else if (equal(c,"map_receptor_interior")) {
        token = GPF_MAP_RECEPTOR_INTERIOR;

    }
    return(token);
}
/* EOF */
