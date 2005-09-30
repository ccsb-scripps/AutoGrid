/* distdepdiel.cpp */
/*
  $Id: distdepdiel.cpp,v 1.1.6.1 2005/09/30 22:43:35 alther Exp $
*/

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

#include "distdepdiel.h"

double calc_ddd_Mehler_Solmajer( double distance, double approx_zero ) {
    /*____________________________________________________________________________
     * Distance-dependent dielectric ewds: Mehler and Solmajer, Prot Eng 4, 903-910.
     *____________________________________________________________________________*/
    double epsilon = 1.0L;
    double lambda = 0.003627L;
    double epsilon0 = 78.4L;
    double A = -8.5525L;
    double B;
    B = epsilon0 - A;
    double rk= 7.7839L;
    double lambda_B;
    lambda_B = -lambda * B;

    epsilon = A + B / (1.0L + rk*exp(lambda_B * distance));

    if (epsilon < approx_zero) {
        epsilon = 1.0L;
    }
    return epsilon;
}
 /* EOF */
