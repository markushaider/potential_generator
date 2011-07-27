#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "variablen.h"
#include "proto.h"

extern char* spath;

void speicher_pot_fits(int erwin, int n)
{
   char  filename1[80];
   kkanow = (float) header1.time;
   FILE *fp;    
   int ii;
   int nx=ZA, ny=ZA, nz=ZA;
   
     
   sprintf(filename1, "%sphi_level_%i_time_%i.bin", spath, n, erwin);
    
   fp=fopen(filename1, "w");

      
   /* do potentials */

   for (ii = 0; ii < nx*ny*nz; ii++)
    {
     Phi[ii]=(Phi[ii])*4.337E14;
     fwrite(&Phi[ii], sizeof(double), 1, fp);
    }

   fclose(fp);

     
   /* anow, total mass */

}
