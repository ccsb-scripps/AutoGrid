/*

 $Id: distdepdiel.cpp,v 1.4 2015/10/02 20:06:03 mp Exp $

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

#include <math.h>
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
        
    // Programming note 2015-10
    //  I suspected a problme in the next line
    //  where I saw  l = 0.003627  not "1" according to manual for AD 3.0.5
    //  as downloaded from http://www.csb.yale.edu/userguides/datamanip/autodock/html/Using_AutoDock_305.c.html
    // Checking out the cited publication, I now believe the "1" is correct.
    // I have not, however, yet located a local copy of the 3.0.5 guide.
    // M Pique
    epsilon = A + B / (1.0L + rk*exp(lambda_B * distance));
    
    if (epsilon < approx_zero) {
        epsilon = 1.0L;
    }
    return epsilon;
}
 /* EOF */
