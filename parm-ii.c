#include "parm-ii.h"

float reqm[NUM_ALL_TYPES] = {
4.00, 4.00, 3.50, 3.20, 4.00, 4.20, 2.00, 1.30, 
3.09, 4.09, 4.33, 4.72, 1.30, 1.48, 1.98, 2.00};

float eps[NUM_ALL_TYPES] = {
0.1500, 0.1500, 0.1600, 0.2000, 0.2000, 0.2000, 0.0200, 0.0100,
0.0800, 0.2760, 0.3890, 0.5520, 0.8750, 0.5500, 0.5500, 0.0200 }; /* *FE_vdW_coeff */

float SolVol[NUM_ALL_TYPES] = {
12.77, 10.80, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
 0.0,   0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

float SolPar[NUM_ALL_TYPES] = {
4.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

float SolCon[NUM_ALL_TYPES] = {
0.0, 0.0, 0.0, 0.236, 0.0, 0.0, 0.118, 0.0, 
0.0, 0.0, 0.0, 0.0,   0.0, 0.0, 0.0,   0.0};

float reqm_Hbond[NUM_ALL_TYPES] = {
0.0, 0.0, 1.9, 1.9, 2.5, 0.0, 0.0, 0.0, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

float eps_Hbond[NUM_ALL_TYPES] = {
0.0, 0.0, 5.0, 5.0, 1.0, 0.0, 0.0, 0.0, 
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; /*  *FE_hbond_coeff  */


/* Functions */
float get_Rij(int type1, int type2)
{
    return (reqm[type1] + reqm[type2])/2.;
}

float get_epsij(int type1, int type2)
{
    return sqrt(eps[type1] * eps[type2]) * FE_vdW_coeff;
}

float get_Rij_Hbond(int type1, int type2)
{
    return (reqm_Hbond[type1] + reqm_Hbond[type2])/2.;
}

float get_epsij_Hbond(int type1, int type2)
{
    return sqrt(eps_Hbond[type1] * eps_Hbond[type2]) * FE_vdW_coeff;
}

float get_SolVol(int type1)
{
    return SolVol[type1];
}

float get_SolPar(int type1)
{
    return SolPar[type1];
}

float get_SolCon(int type1)
{
    return SolCon[type1];
}

/* EOF */
