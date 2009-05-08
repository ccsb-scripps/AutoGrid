/*

 $Id: ranlib.h,v 1.4 2009/05/08 23:17:35 rhuey Exp $

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

/* Prototypes for all user accessible RANLIB routines */
#ifndef _RANLIB_H
#define _RANLIB_H

#include "typedefs.h"

extern void advnst(FourByteLong k);
extern Real genbet(Real aa,Real bb);
extern Real genchi(Real df);
extern Real genexp(Real av);
extern Real genf(Real dfn, Real dfd);
extern Real gengam(Real a,Real r);
extern void genmn(Real *parm,Real *x,Real *work);
extern void genmul(FourByteLong n,Real *p,FourByteLong ncat,FourByteLong *ix);
extern Real gennch(Real df,Real xnonc);
extern Real gennf(Real dfn, Real dfd, Real xnonc);
extern Real gennor(Real av,Real sd);
extern void genprm(FourByteLong *iarray,int larray);
extern Real genunf(Real low,Real high);
extern void getsd(FourByteLong *iseed1,FourByteLong *iseed2);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong ignbin(FourByteLong n,Real pp);
extern FourByteLong ignnbn(FourByteLong n,Real p);
extern FourByteLong ignlgi(void);
extern FourByteLong ignpoi(Real mu);
extern FourByteLong ignuin(FourByteLong low,FourByteLong high);
extern void initgn(FourByteLong isdtyp);
extern FourByteLong mltmod(FourByteLong a,FourByteLong s,FourByteLong m);
extern void phrtsd(char* phrase,FourByteLong* seed1,FourByteLong* seed2);
extern Real ranf(void);
extern void setall(FourByteLong iseed1,FourByteLong iseed2);
extern void setant(FourByteLong qvalue);
extern void setgmn(Real *meanv,Real *covm,FourByteLong p,Real *parm);
extern void setsd(FourByteLong iseed1,FourByteLong iseed2);
extern Real sexpo(void);
extern Real sgamma(Real a);
extern Real snorm(void);
extern Real rcauchy(Real, Real);
extern Real scauchy1(void);

#endif
