#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "variablen.h"
#include "proto.h"

char *lpath, *spath;

int main(int argc, char** argv)
{
  N22 = (2*ZA)*(2*ZA)*(2*ZA);
  N = (ZA)*(ZA)*(ZA);     /*sizes of the arrays*/
  TotNumPart=128*128*128;
   
  int n, otto, erwin, ii;
  lpath=argv[1];
  spath=argv[2];
 
   for (ii=0; ii < 200; ii++ )
   {

   erwin = ii;
        
   dichterechner(erwin);
   /*hier wird mallocated und 1/r gerechnet*/
   getPoten();
   /*hier wird 1/r Fourier transformed*/
   FFTW_dist();
      
   /*hier wird das Dichtefeld Fourier transformed*/
   for (n = 0; n<4; n++) 
       {
	FFTW_rho(n, erwin);

        if (n!=0)
	  {		
	   for (otto=0; otto < n; otto++)
	      {
		    /*der Verschachtler fuer die Nested Grids*/
	            printf("otto: %i --- n: %i\n", otto, n);
		    getbackPoten(otto,n);			  		 
	       }
	 }

             speicher_pot_fits(erwin, n);   
	 
  	     free(Phi);

             printf("juhu\n") ;
        }

           
        free(POSITION_X);
        free(POSITION_Y);
        free(POSITION_Z);
        free(MASSE);
        free(Rho);
        fftw_free(rhoft);
        fftw_free(rho);
        fftw_free(distft);
        fftw_free(backpot);
        fftw_free(dist);
  
      }

 

 
  return 0;
}
