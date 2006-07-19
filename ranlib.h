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
