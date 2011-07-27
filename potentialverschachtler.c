#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "variablen.h"
#include "proto.h"
#include <complex.h>
#include <fftw3.h>


#define   MAXI(x,y)     (((x) > (y)) ? (x) : (y))
#define   MIN(x,y)     (((x) < (y)) ? (x) : (y))

/*------------------------------------------------------------------------------------------*/
void getbackPoten(int m, int q)
{
    int  i, j, k;
    int ll, g, o;
    double Lx=HYDROBOX, sx;
    int nx, ny, nz, pnx, pny, pnz;
    int x;
    kkanow = (float) header1.time;   
    nx=ny=nz=ZA;
    pnx=pny=pnz=2*ZA;
        
    sx=((Lx/(double)ZA)/(double)pow(2,m));    
    
    if (ThisTask == 0)
    {   
    
    o = m*nx*ny*nz;
    for (k=0;k<2*nz;k++)            /*set density field to array*/
        for(j=0;j<2*ny;j++)
            for(i=0;i<2*nx;i++){
            ll = o+(i-nx/2)+nx*(j-ny/2)+nx*ny*(k-nz/2);
            g = i+2*nx*j+4*nx*ny*k; 
            if(k>=nz/2 && k<3*nz/2 && j>=ny/2 && j<3*ny/2
               && i>=nx/2 && i<3*nx/2){
               if(k>=3*nz/4 && k<5*nz/4 && j>=3*ny/4 && j<5*ny/4
               && i>=3*nx/4 && i<5*nx/4) rho[g]=0.0;
                   else rho[g]=(Rho[ll]*sx*sx*sx);
            }
            else rho[g]= 0.0; 
            }
    }
    
    printf("juhu1\n");
    
 /*jetzt nun die das rho Feld Fouriertransformieren*/
 
    
    fftw_plan plan2, iplan2;
                        
    /* next create the plans for FFTW*/
    
     /* create the forward and backward plans: */
     plan2 = fftw_plan_dft_r2c_3d(pnx, pny, pnz, rho, rhoft, FFTW_ESTIMATE);
     iplan2 = fftw_plan_dft_c2r_3d(pnx, pny, pnz, rhoft, rho, FFTW_ESTIMATE);
     
        
     printf("juhu2\n");  
  
           
    /* Now, compute the forward transform: */
    fftw_execute(plan2);
    
     printf("juhu4\n");

    /*Multiplikation im Fourierraum, zwecks Konvulotion*/
    
     for (x=0; x < pnx*pny*pnz; x++)
	{

         rhoft[x]=rhoft[x]*distft[x];
         
        }


      printf("juhu5\n");      
 
     fftw_execute(iplan2);
    
      printf("juhu6\n");
 
    /*hier nun schneiden wir das Potential aus*/
    
    for (k=0;k<2*nz;k++)            /*set density field to array*/
        for(j=0;j<2*ny;j++)
            for(i=0;i<2*nx;i++){
             ll = (i-nx/2)+nx*(j-ny/2)+nx*ny*(k-nz/2);           
             g = i+2*nx*j+4*nx*ny*k;
               if(k>=nz/2 && k<3*nz/2 && j>=ny/2 && j<3*ny/2
               && i>=nx/2 && i<3*nx/2) backpot[ll]=rho[g]/(sx*(double)N22);
            }
   
   /*Interpolationsgetüdel jetzt ----------------------------------*/
    
    printf("juhu7\n");
    int *l, p, t;
    float u,v,w;
    double *intpnt;
          
    intpnt = (double*) malloc(12*sizeof(double));

    l = (int*)malloc(27*sizeof(int));
    
    for (p = 0; p < N; p++) {
        
        i = (int)(floor((posgitter[q*nx*nx*nx+p][0]/sx) + 0.5*(float)(nx-1))); /* get nearest neighbour*/
        j = (int)(floor((posgitter[q*nx*nx*nx+p][1]/sx) + 0.5*(float)(ny-1)));
        k = (int)(floor((posgitter[q*nx*nx*nx+p][2]/sx) + 0.5*(float)(nz-1)));
        
        //printf("%i, %i, %i\n", i,j,k);
                                
        l[0] = i+nx*j+nx*ny*k;
        l[1] = (i-1)+nx*(j+1)+nx*ny*k;
        l[2] = i+nx*(j+1)+nx*ny*k;
        l[3] = (i+1)+nx*(j+1)+nx*ny*k;
        l[4] = (i-1)+nx*j+nx*ny*k;
        l[5] = (i+1)+nx*j+nx*ny*k;
        l[6] = (i-1)+nx*(j-1)+nx*ny*k;
        l[7] = i+nx*(j-1)+nx*ny*k;
        l[8] = (i+1)+nx*(j-1)+nx*ny*k;
        
            
        l[9] = i+nx*j+nx*ny*(k+1);
        l[10] = (i-1)+nx*(j+1)+nx*ny*(k+1);
        l[11] = i+nx*(j+1)+nx*ny*(k+1);
        l[12] = (i+1)+nx*(j+1)+nx*ny*(k+1);
        l[13] = (i-1)+nx*j+nx*ny*(k+1);
        l[14] = (i+1)+nx*j+nx*ny*(k+1);
        l[15] = (i-1)+nx*(j-1)+nx*ny*(k+1);
        l[16] = i+nx*(j-1)+nx*ny*(k+1);
        l[17] = (i+1)+nx*(j-1)+nx*ny*(k+1);
        
        l[18] = i+nx*j+nx*ny*(k-1);
        l[19] = (i-1)+nx*(j+1)+nx*ny*(k-1);
        l[20] = i+nx*(j+1)+nx*ny*(k-1);
        l[21] = (i+1)+nx*(j+1)+nx*ny*(k-1);
        l[22] = (i-1)+nx*j+nx*ny*(k-1);
        l[23] = (i+1)+nx*j+nx*ny*(k-1);
        l[24] = (i-1)+nx*(j-1)+nx*ny*(k-1);
        l[25] = i+nx*(j-1)+nx*ny*(k-1);
        l[26] = (i+1)+nx*(j-1)+nx*ny*(k-1);
        
        u = (float)(abs((posgitter[q*nx*nx*nx+p][2]-posgitter[m*N+l[0]][2])/sx));
        v = (float)(abs((posgitter[q*nx*nx*nx+p][1]-posgitter[m*N+l[0]][1])/sx));
        w = (float)(abs((posgitter[q*nx*nx*nx+p][0]-posgitter[m*N+l[0]][0])/sx));
        
        /*do trilinear interpolation using the procedure of van Leer */
            
        if (posgitter[q*nx*nx*nx+p][2]<posgitter[m*N+l[0]][2]){
        
            for (t=0;t<9;t++) 
            intpnt[t] = backpot[l[t]] - 2*u*(MAXI((backpot[l[t]]-backpot[l[t+9]])*
              (backpot[l[t+18]]-backpot[l[t]]),0)/(backpot[l[t+18]]-backpot[l[t+9]]));
            
            if ((intpnt[11]-intpnt[9])==0)
                    {
                        printf("HILFE\n");
                    }
                                
            if (posgitter[q*nx*nx*nx+p][1]<posgitter[m*N+l[0]][1]){
            intpnt[9] = intpnt[4] - 2*v*(MAXI((intpnt[4]-intpnt[1])*
                    (intpnt[6]-intpnt[4]),0)/(intpnt[6]-intpnt[1]));
            intpnt[10] = intpnt[0] - 2*v*(MAXI((intpnt[0]-intpnt[2])*
                    (intpnt[7]-intpnt[0]),0)/(intpnt[7]-intpnt[2]));
            intpnt[11] = intpnt[5] - 2*v*(MAXI((intpnt[5]-intpnt[3])*
                    (intpnt[8]-intpnt[5]),0)/(intpnt[8]-intpnt[3]));
                    
            }
            
            if (posgitter[q*nx*nx*nx+p][1]>posgitter[m*N+l[0]][1]){
            intpnt[9] = intpnt[4] + 2*v*(MAXI((intpnt[4]-intpnt[1])*
                    (intpnt[6]-intpnt[4]),0)/(intpnt[6]-intpnt[1]));
            intpnt[10] = intpnt[0] + 2*v*(MAXI((intpnt[0]-intpnt[2])*
                    (intpnt[7]-intpnt[0]),0)/(intpnt[7]-intpnt[2]));
            intpnt[11] = intpnt[5] + 2*v*(MAXI((intpnt[5]-intpnt[3])*
                    (intpnt[8]-intpnt[5]),0)/(intpnt[8]-intpnt[3]));
            }
            if (posgitter[q*nx*nx*nx+p][0]<posgitter[m*N+l[0]][0])
            Phi[p] += intpnt[10] - 2*w*(MAXI((intpnt[10]-intpnt[9])*
                    (intpnt[11]-intpnt[10]),0)/(intpnt[11]-intpnt[9]));
            if (posgitter[q*nx*nx*nx+p][0]>posgitter[m*N+l[0]][0])
            Phi[p] += intpnt[10] + 2*w*(MAXI((intpnt[10]-intpnt[9])*
                    (intpnt[11]-intpnt[10]),0)/(intpnt[11]-intpnt[9]));
                    
                    
        }
            if (posgitter[q*nx*nx*nx+p][2]>posgitter[m*N+l[0]][2]){
               
            for (t=0;t<9;t++) 
            intpnt[t] = backpot[l[t]] + 2*u*(MAXI((backpot[l[t]]-backpot[l[t+9]])*
                (backpot[l[t+18]]-backpot[l[t]]),0)/(backpot[l[t+18]]-backpot[l[t+9]]));
                
            if (posgitter[q*nx*nx*nx+p][1]<posgitter[m*N+l[0]][1]){
            intpnt[9] = intpnt[4] - 2*v*(MAXI((intpnt[4]-intpnt[1])*
                    (intpnt[6]-intpnt[4]),0)/(intpnt[6]-intpnt[1]));
            intpnt[10] = intpnt[0] - 2*v*(MAXI((intpnt[0]-intpnt[2])*
                    (intpnt[7]-intpnt[0]),0)/(intpnt[7]-intpnt[2]));
            intpnt[11] = intpnt[5] - 2*v*(MAXI((intpnt[5]-intpnt[3])*
                    (intpnt[8]-intpnt[5]),0)/(intpnt[8]-intpnt[3]));
            }
            if (posgitter[q*nx*nx*nx+p][1]>posgitter[m*N+l[0]][1]){
            intpnt[9] = intpnt[4] + 2*v*(MAXI((intpnt[4]-intpnt[1])*
                    (intpnt[6]-intpnt[4]),0)/(intpnt[6]-intpnt[1]));
            intpnt[10] = intpnt[0] + 2*v*(MAXI((intpnt[0]-intpnt[2])*
                    (intpnt[7]-intpnt[0]),0)/(intpnt[7]-intpnt[2]));
            intpnt[11] = intpnt[5] + 2*v*(MAXI((intpnt[5]-intpnt[3])*
                    (intpnt[8]-intpnt[5]),0)/(intpnt[8]-intpnt[3]));
            }
            
            if (posgitter[q*nx*nx*nx+p][0]<posgitter[m*N+l[0]][0])
            Phi[p] += intpnt[10] - 2*w*(MAXI((intpnt[10]-intpnt[9])*
                    (intpnt[11]-intpnt[10]),0)/(intpnt[11]-intpnt[9]));
            if (posgitter[q*nx*nx*nx+p][0]>posgitter[m*N+l[0]][0])
            Phi[p] += intpnt[10] + 2*w*(MAXI((intpnt[10]-intpnt[9])*
                    (intpnt[11]-intpnt[10]),0)/(intpnt[11]-intpnt[9]));
                    
        }
    
    }
    
 free(intpnt);
 free(l);
  
 printf("juhu8\n");

return;

}
