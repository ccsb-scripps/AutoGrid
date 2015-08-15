/******************************************************************************
**                                                                           **
** MSMS v2.5                                                                 **
**                                                                           **
** Copyright: (C) 1997 Michel F. Sanner and TSRI                             **
**		                                                             **
** ALL RIGHTS RESERVED. THIS SOURCE CODE CAN ONLY BE USED TO PORT MSMS       **
** TO A NEW PLATFORM, AFTER A WRITTEN AGREEMENT WITH THE AUTHOR              **
**                                                                           **
**___________________________________________________________________________**
**                                                                           **
**   Authors: Michel F. Sanner                               Nov. 1999       **
**                                                                           **
**   e-mail: sanner@scripps.edu                                              **
**   www:    http://www.scripps.edu/sanner                                   **
**                                                                           **
**            The Scripps Research Institute                                 **
**            Department of Molecular Biology, MB5                           **
**            10666 North Torrey Pines Road                                  **
**            La Jolla, CA 92037.                                            **
**                                                                           **
**      Date: 11/3/99                                                        **
**__________________________________________________________________________**/

/* $Header: /Users/mp/facil/autodock/git-luna/autogrid-cvstar/bhtree.h,v 1.1 2015/08/15 00:55:21 sanner Exp $
 *
 * $Id: bhtree.h,v 1.1 2015/08/15 00:55:21 sanner Exp $
 *
 * $Log: bhtree.h,v $
 * Revision 1.1  2015/08/15 00:55:21  sanner
 * - first commit
 *
 * Revision 1.1  2003/07/03 01:36:52  sanner
 * Initial revision
 *
 * Revision 0.2  2000/08/14 18:00:28  sanner
 * removed copyright text
 *
 * Revision 0.1  2000/08/14 17:45:35  sanner
 * added copyright text
 *
 * Revision 0.0  1999/10/27 17:51:28  sanner
 * *** empty log message ***
 *
 * Revision 1.3  1998/03/19  22:10:14  sanner
 * fixed the RCS Header Id and Log in source files
 *
 */
/* 
   provided by Armin Widmer
*/
#ifndef BHTREEDEF
#define BHTREEDEF

#include <stdio.h>

typedef struct BHpoint {
  float x[3];
  float r;
  int   at;
} BHpoint;

typedef struct BHnode {
  struct BHnode *left,*right;
  struct BHpoint **atom;
  float  cut;
  int    dim,n;
} BHnode;  

typedef struct BHtree {
  struct BHnode *root;
  struct BHpoint **atom;
  float xmin[3];
  float xmax[3];
  float rm;
#ifdef STATBHTREE
  long tot;    /* total number of neighbors returned by findBHclose */
  int max,min; /* min and max of these numbers */
  int nbr;     /* number of calls to findBHclose */
#endif
  short bfl;
} BHtree;

BHtree *generateBHtree(BHpoint **atoms, int nbat, int granularity);
BHnode *findBHnode(BHtree *tree,float *x);
int    findBHcloseAtoms(BHtree *tree,float *x,float cutoff,
		        int *atom,int maxn);
int    findBHcloseAtomsdist(BHtree *tree,float *x,float cutoff,
		            int *atom,float *d,int maxn);
void   freeBHtree(BHtree *tree);
void   divideBHnode(BHnode *node,float *xmin,float *xmax,int granularity);
void   freeBHnode(BHnode *node);

extern BHtree *bht;
extern BHpoint **BHat;

#endif




