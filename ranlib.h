/* Prototypes for all user accessible RANLIB routines */
#ifndef _RANLIB_H
#define _RANLIB_H

#include "typedefs.h"

extern void advnst(FourByteLong k);
extern FloatOrDouble genbet(FloatOrDouble aa,FloatOrDouble bb);
extern FloatOrDouble genchi(FloatOrDouble df);
extern FloatOrDouble genexp(FloatOrDouble av);
extern FloatOrDouble genf(FloatOrDouble dfn, FloatOrDouble dfd);
extern FloatOrDouble gengam(FloatOrDouble a,FloatOrDouble r);
extern void genmn(FloatOrDouble *parm,FloatOrDouble *x,FloatOrDouble *work);
extern void genmul(FourByteLong n,FloatOrDouble *p,FourByteLong ncat,FourByteLong *ix);
extern FloatOrDouble gennch(FloatOrDouble df,FloatOrDouble xnonc);
extern FloatOrDouble gennf(FloatOrDouble dfn, FloatOrDouble dfd, FloatOrDouble xnonc);
extern FloatOrDouble gennor(FloatOrDouble av,FloatOrDouble sd);
extern void genprm(FourByteLong *iarray,int larray);
extern FloatOrDouble genunf(FloatOrDouble low,FloatOrDouble high);
extern void getsd(FourByteLong *iseed1,FourByteLong *iseed2);
extern void gscgn(FourByteLong getset,FourByteLong *g);
extern FourByteLong ignbin(FourByteLong n,FloatOrDouble pp);
extern FourByteLong ignnbn(FourByteLong n,FloatOrDouble p);
extern FourByteLong ignlgi(void);
extern FourByteLong ignpoi(FloatOrDouble mu);
extern FourByteLong ignuin(FourByteLong low,FourByteLong high);
extern void initgn(FourByteLong isdtyp);
extern FourByteLong mltmod(FourByteLong a,FourByteLong s,FourByteLong m);
extern void phrtsd(char* phrase,FourByteLong* seed1,FourByteLong* seed2);
extern FloatOrDouble ranf(void);
extern void setall(FourByteLong iseed1,FourByteLong iseed2);
extern void setant(FourByteLong qvalue);
extern void setgmn(FloatOrDouble *meanv,FloatOrDouble *covm,FourByteLong p,FloatOrDouble *parm);
extern void setsd(FourByteLong iseed1,FourByteLong iseed2);
extern FloatOrDouble sexpo(void);
extern FloatOrDouble sgamma(FloatOrDouble a);
extern FloatOrDouble snorm(void);
extern FloatOrDouble rcauchy(FloatOrDouble, FloatOrDouble);
extern FloatOrDouble scauchy1(void);

#endif
