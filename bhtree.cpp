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

/* $Header: /Users/mp/facil/autodock/git-luna/autogrid-cvstar/bhtree.cpp,v 1.2 2015/10/02 20:03:29 mp Exp $
 *
 * $Id: bhtree.cpp,v 1.2 2015/10/02 20:03:29 mp Exp $
 *
 * $Log: bhtree.cpp,v $
 * Revision 1.2  2015/10/02 20:03:29  mp
 * Changed 'include "mystdlib.h"' to 'include <stdlib.h>'
 *
 * Revision 1.1  2015/08/15 00:55:21  sanner
 * - first commit
 *
 * Revision 1.1  2003/07/03 01:36:52  sanner
 * Initial revision
 *
 * Revision 0.2  2000/08/14 18:00:34  sanner
 * removed copyright text
 *
 * Revision 0.1  2000/08/14 17:45:19  sanner
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bhtree.h"

/* Barnes - Hut Trees ? */

#define NSTEPS 128

static int findBHcloseAtomsInNode(BHnode *node,float *x,float cutoff,
				  int *atom,int maxn);
static int findBHcloseAtomsInNodedist(BHnode *node,float *x,float cutoff,
				      int *atom,float *dist,int maxn);

/*-------------------- generateBHtree --------------------*/

BHtree *generateBHtree(BHpoint **atoms,int nbat,int granularity)
{
  /* 3D space-sort atoms into Barnes-Hut tree with given granularity */
  /* (max number of atoms per leaf-node) */

  BHtree *r;
  BHpoint **p;
  int i,k;

  /* allocate tree data structure */

  r=(BHtree *)malloc(sizeof(BHtree));
  if (r==NULL) return(r);
  r->atom=(BHpoint **)NULL;
  r->bfl=0;
  r->rm=0.0;
  for (i=0;i<nbat;i++)
    if (r->rm<atoms[i]->r) r->rm = atoms[i]->r;
  r->rm += 0.1;
#ifdef STATBHTREE
  r->tot=0;    /* total number of neighbors returned by findBHclose */
  r->max=0; r->min=9999999; /* min and max of these numbers */
  r->nbr=0;     /* number of calls to findBHclose */
#endif

  /* allocate root node data structure */

  r->root=(BHnode *)malloc(sizeof(BHnode));
  if (r->root==NULL) {
    freeBHtree(r);
    return((BHtree *)NULL);
  }

  /* initialize root node data structure */

  r->root->atom=(BHpoint **)NULL;
  r->root->n=0;
  r->root->dim= -1;
  r->root->left=(BHnode *)NULL;
  r->root->right=(BHnode *)NULL;

  /* count atoms in chain */

  if (nbat==0) {
    freeBHtree(r);
    return((BHtree *)NULL);
  }
    
  /* allocate atom pointer array */

  r->atom=atoms;
  if (r->atom == NULL) { 
    freeBHtree(r);
    return((BHtree *)NULL);
  }
  r->root->atom=r->atom;

  /* fill atom pointer array */

  p=r->root->atom;
  r->root->n=nbat;

  /* determine box dimension */

  p=r->root->atom;

  for (k=0;k<3;k++) {
    r->xmin[k]=p[0]->x[k];
    r->xmax[k]=r->xmin[k];
  }

  for (i=1;i<r->root->n;i++) {
    for (k=0;k<3;k++) {
      if (r->xmin[k] > p[i]->x[k]) r->xmin[k] = p[i]->x[k];
      if (r->xmax[k] < p[i]->x[k]) r->xmax[k] = p[i]->x[k];
    }
  }

  /* start recursive 3D space sort */

  divideBHnode(r->root,r->xmin,r->xmax,granularity);

  /* done... */

  return(r);
}
/*-------------------- find_BHnode --------------------*/

BHnode *findBHnode(BHtree *tree,float *x)
{
  /* find leaf in BH tree that contains 3D point x */

  BHnode *r;
  int k;

  if (tree==NULL) return((BHnode *)NULL);

  /* first check if point is in tree */

  for (k=0;k<3;k++) {
    if (x[k]<tree->xmin[k]) return((BHnode *)NULL);
    if (x[k]>tree->xmax[k]) return((BHnode *)NULL);
  }

  /* if point is in tree, locate the leaf */

  r=tree->root;
  while (r!=NULL) {
    if (r->dim<0) break;
    if (x[r->dim]<r->cut) r=r->left;
    else r=r->right;
  }

  return(r);
}
/*-------------------- freeBHtree --------------------*/

void freeBHtree(BHtree *tree)
{
  int i;

  if (tree->atom!=NULL) {
    for (i=0;i<tree->root->n;i++) free(tree->atom[i]);
    free(tree->atom);
  }
  freeBHnode(tree->root);
  free(tree);
}
/*-------------------- freeBHnode --------------------*/

void freeBHnode(BHnode *node)
{
  if (node!=NULL) {
    freeBHnode(node->left);
    freeBHnode(node->right);
    free(node);
  }
}
/*-------------------- divideBHtree --------------------*/

void divideBHnode(BHnode *node,float *xmin,float *xmax,int granularity)
{
  float cut,dx,xminl[3],xmaxl[3],xminr[3],xmaxr[3];
  int dim,i,j,k,n[NSTEPS],lm,rm;
  BHpoint *a;

  /* if nothing left to divide just return */

  if (node==NULL) return;
  if (granularity<1 || node->n <= granularity) return;
  if (node->atom==NULL) return;

  /* determine dimension along which to cut */

  dim=0;
  if (xmax[1]-xmin[1] > xmax[dim]-xmin[dim]) dim=1;
  if (xmax[2]-xmin[2] > xmax[dim]-xmin[dim]) dim=2;

  /* determine position of cutting plane */

  dx=(xmax[dim]-xmin[dim])/NSTEPS;
  if (dx<0.0001) return;
  for (i=0;i<NSTEPS;i++) n[i]=0;
  for (j=0;j<node->n;j++) {
    i=(node->atom[j]->x[dim]-xmin[dim])/dx;
    if (i>=0 && i<NSTEPS) n[i]++;
  }
  for (i=1;i<NSTEPS;i++) {
    n[i]+=n[i-1];
    if (n[i]>node->n/2) break;
  }
  cut=xmin[dim]+i*dx;
  if (n[i]>=node->n) return;

  /* create left/right descendants */

  node->left=(BHnode *)malloc(sizeof(BHnode));
  if (node->left==NULL) return;
  node->left->dim= -1;
  node->left->left=(BHnode *)NULL;
  node->left->right=(BHnode *)NULL;
        
  node->right=(BHnode *)malloc(sizeof(BHnode));
  if (node->right==NULL) {
    freeBHnode(node->left);
    return;
  }
  node->right->dim= -1;
  node->right->left=(BHnode *)NULL;
  node->right->right=(BHnode *)NULL;

  node->cut=cut;
  node->dim=dim;

  /* sort atoms into left/right descendants */

  lm=0;
  rm=node->n-1;

  while (lm<rm) {
    for(;lm<node->n;lm++) if (node->atom[lm]->x[dim]>=cut) break;
    for(;rm>=0;rm--)      if (node->atom[rm]->x[dim]<cut) break;
    if (lm<rm) {
      a=node->atom[rm];
      node->atom[rm]=node->atom[lm];
      node->atom[lm]=a;
      rm--;
      lm++;
    }
  }

  if (lm==rm) {
    if (node->atom[rm]->x[dim]<cut) lm++;
    else rm--;
  }

  node->left->n=rm+1;
  node->left->atom=node->atom;

  node->right->n=node->n-(rm+1);
  node->right->atom=node->atom+lm;

  /* if descendants are coarse, cut them up... */

  if (node->left->n > granularity) {
    for (k=0;k<3;k++) {
      xminl[k]=xmin[k];
      xmaxl[k]=xmax[k];
    }
    xmaxl[dim]=cut;
    divideBHnode(node->left,xminl,xmaxl,granularity);
  }

  if (node->right->n > granularity) {
    for (k=0;k<3;k++) {
      xminr[k]=xmin[k];
      xmaxr[k]=xmax[k];
    }
    xminr[dim]=cut;
    divideBHnode(node->right,xminr,xmaxr,granularity);
  }

  /* done... */

  return;

}
/*-------------------- findBHcloseAtoms --------------------*/

int findBHcloseAtomsdist(BHtree *tree,float *x,float cutoff,
		     int *atom,float *dist,int maxn)
{
  int i;

  if (maxn<1 || tree==NULL || cutoff <=0.0) return(0);
  if (tree->root==NULL) return(0);

  for (i=0;i<3;i++) {
    if (x[i] < tree->xmin[i] - cutoff) break;
    if (x[i] > tree->xmax[i] + cutoff) break;
  }
  if (i<3) return(0);

  return(findBHcloseAtomsInNodedist(tree->root,x,cutoff,atom,dist,maxn));
}
/*-------------------- findBHcloseAtomsInNode --------------------*/

static int findBHcloseAtomsInNodedist(BHnode *node,float *x,float cutoff,
				  int *atom,float *dist,int maxn)
{
  int j,n;
  float d[3],D,C;

  if (node==NULL) return(0);
  if (maxn<=0) return(0);
  if (node->n < 1) return(0);
  
  if (node->dim<0) {
    n=0;
    C=cutoff*cutoff;
    for (j=0;j<node->n;j++) {
      d[0]=x[0]-node->atom[j]->x[0];
      if (d[0]>cutoff || d[0]< -cutoff) continue;
      d[1]=x[1]-node->atom[j]->x[1];
      if (d[1]>cutoff || d[1]< -cutoff) continue;
      d[2]=x[2]-node->atom[j]->x[2];
      if (d[2]>cutoff || d[2]< -cutoff) continue;
      D=d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
      if (D>C) continue;
      if (n<maxn) {
	atom[n]=node->atom[j]->at;
	dist[n]=sqrt(D);
	n++;
      }else{
	n++;
	break;
      }
    }
  }else{
    n=0;
    if (x[node->dim]<node->cut+cutoff) {
      n+=findBHcloseAtomsInNodedist(node->left,x,cutoff,atom,dist,maxn);
    }
    if (x[node->dim]>=node->cut-cutoff) {
      n+=findBHcloseAtomsInNodedist(node->right,x,cutoff,atom+n,dist+n,maxn-n);
    }
  }
  return(n);
}

/*-------------------- findBHcloseAtoms --------------------*/

int findBHcloseAtoms(BHtree *tree,float *x,float cutoff,
		     int *atom,int maxn)
{
  int i,n;

  if (maxn<1 || tree==NULL || cutoff <=0.0) return(0);
  if (tree->root==NULL) return(0);

  for (i=0;i<3;i++) {
    if (x[i] < tree->xmin[i] - cutoff) break;
    if (x[i] > tree->xmax[i] + cutoff) break;
  }
  if (i<3) return(0);

#ifdef STATBHTREE

  n = findBHcloseAtomsInNode(tree->root,x,cutoff,atom,maxn);
  tree->nbr++;
  tree->tot += n;
  if (n>tree->max) tree->max=n;
  if (n<tree->min) tree->min=n;
  return(n);

#else
  return(findBHcloseAtomsInNode(tree->root,x,cutoff,atom,maxn));
#endif
}
/*-------------------- findBHcloseAtomsInNode --------------------*/

static int findBHcloseAtomsInNode(BHnode *node,float *x,float cutoff,
				  int *atom, int maxn)
{
  int j,n;
  float C,D;
  double d1,d2,d3;
  BHpoint *p;
/*
  if (node==NULL) return(0);
  if (maxn<=0) return(0);
  if (node->n < 1) return(0);
*/  
  if (node->dim<0) {

    C=cutoff*cutoff;

    for (n=j=0;j<node->n;j++) {
      p = node->atom[j];
      d1 = x[0] - p->x[0];
      if (d1 > cutoff || d1 < -cutoff) continue;
      d2 = x[1] - p->x[1];
      if (d2 > cutoff || d2 < -cutoff) continue;
      d3 = x[2] - p->x[2];
      if (d3 > cutoff || d3 < -cutoff) continue;

      D=d1*d1 + d2*d2 + d3*d3;

      if (D>C) continue;
      if (n<maxn) {
	atom[n] = p->at;
	n++;
      }else{
	printf("ERROR: findBHcloseAtomsInNode: result array too small\n");
	break;
      }
    }
  } else {
    n=0;
    if (x[node->dim]<node->cut+cutoff) {
      n+=findBHcloseAtomsInNode(node->left,x,cutoff,atom,maxn);
    }
    if (x[node->dim]>=node->cut-cutoff) {
      n+=findBHcloseAtomsInNode(node->right,x,cutoff,atom+n,maxn-n);
    }
  }
  return(n);
}
