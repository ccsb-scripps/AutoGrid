#include <stdio.h> 
#include <math.h>
#include "parm.h"

/*indices:
         0    1    2    3    4    5    6    7    8    9   10   11   12   13   14  
        C3   C2   C1  Cac   Cpl N3pl Nox   N3  Ntr  Npl   N1  Nam   O3   O2   Om 
 
        15   16   17   18   19   20   21   22   23   24   25   26   27   28   29  
      S3pl   S3   S2  Sac  Sox    S   HC    H   P3  Pac  Pox    B  Bac  Box   Al  

        30   31   32   33   34   35   36   37   38   39   40   41   42   43   44  
        As   Be   Br   Ca   Cl   Cu   Fl   Fe   Ge    I    K   Li   Mg   Mn   Na   

        45   46   47   48
        Ni   Pb   Si   Zn
*/

float FE_vdW_coeff   = 0.1485;
float FE_hbond_coeff = 0.656 ;

float Rij[num_atom_types] = { 
        4.,  4.,  4.,  4.,  4., 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.2, 3.2, 3.2, 
        4.,  4.,  4.,  4.,  4.,  4.,  2.,  2., 4.2, 4.2, 4.2,  0.,  0.,  0., 0., 
        0.,  0.,4.33,  0.,4.09,  0.,3.09, 1.3,  0.,  0.,  0.,  0.,  0., 0., 0., 
        0.,  0.,  0.,  0., };


/* to be multiplied by FE_vdW_coeff =  0.1485*/
float eps[num_atom_types] = {
       .15, .15, .15, .15, .15, .16, .16, .16, .16, .16, .16, .16, .20, .20,.20,
       .20, .20, .20, .20, .20, .20, .02, .02, .20, .20, .20,  0.,  0., 0., 0.,
        0.,  0.,.389,  0.,.276,  0., .08, .01,  0.,  0.,  0.,  0.,  0., 0., 0.,
        0.,  0.,  0.,  0.,};

/*indices:
            0    1    2    3    4    5    6    7    8    9   
          N3+  Nox   N3  Ntr  Npl   N1  Nam   O3   O2   O-  
 
           10   11   12   13   14  15  
          S3+   S3   S2  Sac  Sox   S 
*/

float Rij_h[num_hbnd_types] = {
          1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9,
          2.5, 2.5, 2.5, 2.5, 2.5, 2.5};

/* to be multiplied by FE_hbond_coeff =  0.0656*/
float eps_h[num_hbnd_types]  = {
          5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
          1.0, 1.0, 1.0, 1.0, 1.0, 1.0};


float get_Rij(int type1,  int type2)
    {
        float val1,  val2;
        val1 = Rij[type1];
        val2 = Rij[type2];
        return (val1+val2)/2.;
    }


float get_epsij(int type1,  int type2)
    {
        float val1,  val2;
        val1 = eps[type1] * FE_vdW_coeff;
        val2 = eps[type2] * FE_vdW_coeff;
        return sqrt(val1*val2);
    }


/*for testing purposes*/
float main(int argc, char **argv)
{
        int val1, val2;
        printf("in main\n");
        if(argc<2){
            printf("sorry, need two input values\n");
            return 1;
        } 
        sscanf(argv[1],"%d", &val1);
        printf("%d ", val1);
        sscanf(argv[2],"%d", &val2);
        printf("%d\n", val2);
        printf("Rij=%f\n", get_Rij(val1,val2));
        printf("epsij=%f\n", get_epsij(val1,val2));
        /*printf("->->->%d<-<-<-\n", Cac);*/
        printf("get_Rij(C3,C3) = %f\n", get_Rij(0,0));
        printf("get_epsij(C3,C3) = %f\n", get_epsij(0,0));
        return 1;
}
