#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "variablen.h"
#include "proto.h"
#include <complex.h>
#include <fftw3.h>


/*-------------------------------------------------------------------------------------*/
void FFTW_rho(int n, int erwin)
{
    kkanow = (float) header1.time;
    
    int x;
    fftw_plan plan1, iplan1;
        
    int i, j, k;
    int l, g, o;
    
    //FILE *fd;
    //char filename[80];
       
    int N = (ZA)*(ZA)*(ZA), N22=N*8;

    int nx, ny, nz, pnx, pny, pnz;
    double Lx, Ly, Lz;
    double sx, sy, sz;
    
    Lx=Ly=Lz=HYDROBOX;
    nx=ny=nz=ZA;
    pnx=pny=pnz=2*ZA;
    
    
    sx=((Lx/(double)nx)/(double)pow(2,n));                    
    sy=((Ly/(double)ny)/(double)pow(2,n));
    sz=((Lz/(double)nz)/(double)pow(2,n));

     //printf("Servusoi\n"); 

    Phi=(double*) malloc(sizeof(double)*N);

     //printf("Servusio\n");   

     /* create the forward and backward plans: */
     plan1 = fftw_plan_dft_r2c_3d(pnx, pny, pnz, rho, rhoft, FFTW_ESTIMATE);
     iplan1 = fftw_plan_dft_c2r_3d(pnx, pny, pnz, rhoft, rho, FFTW_ESTIMATE);
    
    
        o = n*nx*ny*nz;
        for (k=0;k<2*nz;k++)            /*set density field to array*/
            for(j=0;j<2*ny;j++)
                for(i=0;i<2*nx;i++){
                l = o+(i-nx/2)+nx*(j-ny/2)+nx*ny*(k-nz/2);
                g = i+2*nx*j+4*nx*ny*k;
                if(k>=2*nz/4 && k<6*nz/4 && j>=2*ny/4 && j<6*ny/4
                    && i>=2*nx/4 && i<6*nx/4) rho[g]= (Rho[l]*sx*sy*sz);
                else rho[g]= 0.0;
            }
      
     printf("Servus3\n");

     /* sprintf(filename, "nbody_rho_level_%i_time_%i.bin", n, erwin);         

        fd = fopen(filename,"w");

   	for (k=0;k<2*nz;k++)            
            for(j=0;j<2*ny;j++)
                for(i=0;i<2*nx;i++){
   	    		 if(k>=2*nz/4 && k<6*nz/4 && j>=2*ny/4 && j<6*ny/4
                    	 && i>=2*nx/4 && i<6*nx/4) 
				{
				rho[i+2*nx*j+4*nx*ny*k]=rho[i+2*nx*j+4*nx*ny*k]*kkanow*kkanow*kkanow;
				fwrite(&rho[i+2*nx*j+4*nx*ny*k], sizeof(double), 1, fd);
				}
   	     }
        fclose(fd);  */
    
        


    /******************************************************************/
          
     
    /* Now, compute the forward transform: */
    fftw_execute(plan1);
           
      printf("Servus4\n");      
    /*Multiplikation im Fourierraum, zwecks Konvulotion*/

    for (x=0; x < pnx*pny*pnz; x++)
	{

         rhoft[x] = distft[x] * rhoft[x];
         
        }  
   
   printf("Servus5\n");   

   fftw_execute(iplan1);
     
   for (k=0;k<2*nz;k++)          /*set density field to array*/
       for(j=0;j<2*ny;j++)
           for(i=0;i<2*nx;i++){
           l = (i-nx/2)+nx*(j-ny/2)+nx*ny*(k-nz/2);
           g = i+2*nx*j+4*nx*ny*k;
           if(k>=2*nz/4 && k<6*nz/4 && j>=2*ny/4 && j<6*ny/4
              && i>=2*nx/4 && i<6*nx/4) Phi[l]=rho[g]/(sx*(double)N22);
           } 

     printf("Servus6\n");   

   fftw_destroy_plan(plan1);
   fftw_destroy_plan(iplan1); 

     printf("Servus7\n");   
      
   return;
}
