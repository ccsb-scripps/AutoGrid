#define C3      0
#define C2      1
#define C1      2
#define Cac     3
#define Cpl     4
    
#define N3pl    5
#define Nox     6
#define N3      7
#define Ntr     8
#define Npl     9
#define N1     10
#define Nam    11

#define O3     12
#define O2     13
#define Om     14

#define S3pl   15
#define S3     16
#define S2     17
#define Sac    18
#define Sox    19
#define S      20

#define HC     21 
#define H      22 

#define P3     23
#define Pac    24 
#define Pox    25

#define B      26
#define Bac    27
#define Box    28

#define Al     29
#define As     30
#define Be     31
#define Br     32
#define Ca     33
#define Cl     34
#define Cu     35  
#define Fl     36  
#define Fe     37 
#define Ge     38  
#define I      39 
#define K      40 
#define Li     41 
#define Mg     42  
#define Mn     43
#define Na     44
#define Ni     45
#define Pb     46
#define Si     47
#define Zn     48

#define num_atom_types     49
#define num_hbnd_types     16


extern float get_Rij(int type1, int type2);
extern float get_epsij(int type1, int type2);

