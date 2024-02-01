/*

 $Id: typedefs.h,v 1.14 2020/05/04 17:00:40 mp Exp $

 AutoDock  

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

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

#ifndef _TYPEDEFS_H
#define _TYPEDEFS_H

/******************************************************************************
 *      Name: typedefs.h                                                      *
 *  Function: Defines types used in Molecular Applications.                   *
 *Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: JAN/18/2003                                                     *
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 01/18/03 GMM     This header added                                         *
 ******************************************************************************/


#ifdef USE_INT_AS_FOURBYTELONG
    typedef int  FourByteLong;
    typedef unsigned int UnsignedFourByteLong;
#define FBL_FMT "%d"
#define UFBL_FMT "%u"
#else
    typedef long FourByteLong;
    typedef unsigned long UnsignedFourByteLong;
#define FBL_FMT "%ld"
#define UFBL_FMT "%lu"
#endif

#ifdef USE_DOUBLE
    typedef double Real;
#   define FDFMT "%lf"
#else
    typedef float Real;
#   define FDFMT "%f"
#endif
#define FDFMT2 FDFMT " " FDFMT
#define FDFMT3 FDFMT " " FDFMT " " FDFMT

// MP note: "const2" tests were "const Real", "const3" tests were "const Real&"
// "const4" tests (early November 2010) were "const Real&"
// MP note: this type is for scalar declarations only, use "const Real a[NN];"
//   for array definitions.
#define ConstReal	const Real&
#define ConstDouble	const double&

typedef float MapType; //the grid map type can be changed independent of other 'Real's
#   define MAPFMT "%f"
//
// the 1-D e_vdW_Hb internal energy table type - this table can range into gigabytes
//  when MAX_ATOM_TYPES is large and the 8-Angstrom NB cutoff is not used.
typedef float ETableType; // for_evdW_Hb


#ifdef USE_VELOCITY_ENGINE
typedef union
{
	vector float vec;
	float elements[4];
} Float4;
#endif

#endif
/* EOF */
