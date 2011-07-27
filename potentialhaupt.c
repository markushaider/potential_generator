#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "variablen.h"
#include "proto.h"
#include <fftw3.h>


/*------------------------------------------------------------------------------------------*/
//----------------------------------------------------------------------------------
   /*getPoten: top routine to calculate potential on nested grid*/
//-----------------------------------------------------------------------------------

void getPoten(void)
{
     int qq, nn, mm, qm, nm, m2, z=0;
     int nx, ny, nz;
     N22 = (2*ZA)*(2*ZA)*(2*ZA);
     N = (ZA)*(ZA)*(ZA);     /*sizes of the arrays*/    
    
    dist = (double*)fftw_malloc(sizeof(double)*N22); /*allocate arrays*/

    if(dist == NULL)
        { 
            printf("Not enough memory!!!\n");
            exit(0);
        }
        
    backpot = (double*)fftw_malloc(sizeof(double)*N);
       if(backpot == NULL){ 
        printf("Not enough memory!!!\n");
        exit(1);
      } 
    
    distft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N22);
      if(distft == NULL){ 
        printf("Not enough memory!!!\n");
        exit(2);
      }
    
    rho = (double*)fftw_malloc(sizeof(double)*N22);  /*allocate memory*/
      if(rho == NULL){ 
        printf("Not enough memory!!!\n");
        exit(0);
      }
        
    rhoft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N22);
      if(rhoft == NULL){ 
        printf("Not enough memory!!!\n");
        exit(2);
      }
 
 nx=ny=nz=ZA;
    
 

            for (qq=0; qq<2*nz; qq++)           /*set the 1/r kernel*/
               for (nn=0; nn<2*ny;nn++)
                for (mm=0;mm<2*nx;mm++){
                    if (qq<=(2*nz) - qq) qm = qq;
                    else qm = ((2*nz) - qq);
                    if (nn<=(2*ny) - nn) nm = nn;
                    else nm = ((2*ny) - nn);
                    if (mm<=(2*nx) - mm) m2 = mm;
                    else m2 = ((2*nx) - mm);
                        dist[z]=(double)(-1/(sqrt((double)(qm*qm+nm*nm+m2*m2)+(1.0/2.0)*(1.0/2.0))));
                        z++;        
                    }
        
        
     
    return;
            
}
