#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "variablen.h"
#include "proto.h"


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
/*------------------------------------------------------------------------------------------*/
int dichterechner(int kk)
{
  POSITION_X = (float*) malloc(TotNumPart*sizeof(float));
  POSITION_Y = (float*) malloc(TotNumPart*sizeof(float));
  POSITION_Z = (float*) malloc(TotNumPart*sizeof(float));
  
  VELOCITY_X = (float*) malloc(TotNumPart*sizeof(float));
  VELOCITY_Y = (float*) malloc(TotNumPart*sizeof(float));
  VELOCITY_Z = (float*) malloc(TotNumPart*sizeof(float));

  MASSE = (float*) malloc(TotNumPart*sizeof(float));
  load_snapshot(kk);
  
  free(Id);
  free(P2);
  
  CICmethod(kk);     /*   use cloud-in-cell method to calculate the density field*/

  /*velocity_rechner*/
  free(VELOCITY_X);
  free(VELOCITY_Y);
  free(VELOCITY_Z);
  return 0;
}

/*------------------------------------------------------------------------------------------*/



/*------------------------------------------------------------------------------------------*/
void CICmethod(int erwin)
{
    
    /*-------------------------------------------*/
    /*erzeugen der Gittereigenschaften*/
    
    int n, m, l, nlevels=4, k, j, i;
    int nx, ny, nz;
    float Lx, Ly, Lz;
    float sx, sy, sz;
    double *otto_tmp;
    kkanow = (float) header1.time;
    printf("kkanow: %f, header1.time: %g\n", kkanow, header1.time);
    Lx=Ly=Lz=HYDROBOX;
    
    nx=ny=nz=ZA;
        
    sx=(Lx/(float)nx);
    sy=(Ly/(float)ny);
    sz=(Lz/(float)nz);
    velox= (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
    veloy= (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
    veloz= (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
    gwichterl = (float*) malloc(ZA*ZA*ZA*nlevels*sizeof(float));

    for (n=0;n<nlevels;n++)
        {
        m = n*nx*ny*nz;
        for (k=0;k<nz;k++)
        for (j=0;j<ny;j++)
        for (i=0;i<nx;i++)
            {
            l=m+i+nx*j+nx*ny*k;
            posgitter[l][0] = sx*(((float)i)-0.5*(float)(nx-1));
            
            posgitter[l][1] = sy*(((float)j)-0.5*(float)(ny-1));
                           
            posgitter[l][2] = sz*(((float)k)-0.5*(float)(nz-1));
                            
            }
       sx/=2.0;
       sy/=2.0;
       sz/=2.0;
        }
    /*-------------------------------------------*/
    
    int o,p;
    
    float *Weight;
    Weight = (float*) malloc((3 + 1) * sizeof(float));
    int nmp=ZA*ZA*ZA*nlevels;
    
    sx=(Lx/(float)nx);
    sy=(Ly/(float)ny);
    sz=(Lz/(float)nz);
    
    Rho = (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
    otto_tmp = (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));

    for (m = 0; m < nmp; m++)
        {
          Rho[m] = 0.0;  /* set potential to zero*/
          velox[m] = 0.0;  /* set to zero*/
          veloy[m] = 0.0;  /* set to zero*/
          veloz[m] = 0.0;  /* set to zero*/
          gwichterl[m] = 0.0;  /* set to zero*/
        }    
    
    for (n=0;n<nlevels;n++)        /*loop over levels*/
    {   
        for (p = 0; p < TotNumPart; p++)   /*loop over all bodies*/
         {
            if (MASSE[p] == 0.0) continue;     /*skip dead particles*/
    
            i = (int)(floor((((POSITION_X[p]/sx) + 0.5*(float)(nx-1))))); /* get nearest neighbour*/
            j = (int)(floor((((POSITION_Y[p]/sy) + 0.5*(float)(ny-1)))));
            k = (int)(floor((((POSITION_Z[p]/sz) + 0.5*(float)(nz-1)))));
            
            if ((i < 0) || (i >= (nx-1)) || (j < 0) || (j >= (ny-1)) 
                    || (k < 0) || (k >= (nz-1))) continue; /*grid point within grid?*/

            o = n*nx*ny*nz;       
            
            l = o+i+nx*j+nx*ny*k;
            getWeights( l, p, Weight, sx);  /* calculate weights*/
            Rho[l] += (double)((MASSE[p] * Weight[3])/(sx*sy*sz)); /* calculate density*/
            velox[l] += (double)VELOCITY_X[p]*Weight[3]*MASSE[p]; 
            veloy[l] += (double)VELOCITY_Y[p]*Weight[3]*MASSE[p];
            veloz[l] += (double)VELOCITY_Z[p]*Weight[3]*MASSE[p];
            gwichterl[l] += MASSE[p]*Weight[3];
        
            //printf("Weight:%f------Masse:%f------Velocity:%f-----POSITION:%f\n", Weight[3], MASSE[p], VELOCITY_X[p], POSITION_X[p]);

            l = o+(i+1)+nx*j+nx*ny*k;
            getWeights( l, p, Weight, sx);  /* calculate weights*/
            Rho[l] += (double)((MASSE[p] * Weight[3])/(sx*sy*sz)); /* calculate density*/ 
            velox[l] += (double)VELOCITY_X[p]*Weight[3]*MASSE[p]; 
            veloy[l] += (double)VELOCITY_Y[p]*Weight[3]*MASSE[p];
            veloz[l] += (double)VELOCITY_Z[p]*Weight[3]*MASSE[p];
            gwichterl[l] += MASSE[p]*Weight[3];
 
            l = o+i+nx*(j+1)+nx*ny*k;
            getWeights( l, p, Weight, sx);  /* calculate weights*/
            Rho[l] += (double)((MASSE[p] * Weight[3])/(sx*sy*sz)); /* calculate density*/         
            velox[l] += (double)VELOCITY_X[p]*Weight[3]*MASSE[p]; 
            veloy[l] += (double)VELOCITY_Y[p]*Weight[3]*MASSE[p];
            veloz[l] += (double)VELOCITY_Z[p]*Weight[3]*MASSE[p];
            gwichterl[l] += MASSE[p]*Weight[3];
 
            l = o+i+nx*j+nx*ny*(k+1);
            getWeights( l, p, Weight, sx);  /* calculate weights*/
            Rho[l] += (double)((MASSE[p] * Weight[3])/(sx*sy*sz)); /* calculate density*/
            velox[l] += (double)VELOCITY_X[p]*Weight[3]*MASSE[p]; 
            veloy[l] += (double)VELOCITY_Y[p]*Weight[3]*MASSE[p];
            veloz[l] += (double)VELOCITY_Z[p]*Weight[3]*MASSE[p];
            gwichterl[l] += MASSE[p]*Weight[3];

            l = o+(i+1)+nx*(j+1)+nx*ny*k;
            getWeights( l, p, Weight, sx);  /* calculate weights*/
            Rho[l] += (double)((MASSE[p] * Weight[3])/(sx*sy*sz)); /* calculate density*/
            velox[l] += (double)VELOCITY_X[p]*Weight[3]*MASSE[p]; 
            veloy[l] += (double)VELOCITY_Y[p]*Weight[3]*MASSE[p];
            veloz[l] += (double)VELOCITY_Z[p]*Weight[3]*MASSE[p];
            gwichterl[l] += MASSE[p]*Weight[3];
 
            l = o+(i+1)+nx*j+nx*ny*(k+1);
            getWeights( l, p, Weight, sx);  /* calculate weights*/
            Rho[l] += (double)((MASSE[p] * Weight[3])/(sx*sy*sz)); /* calculate density*/
            velox[l] += (double)VELOCITY_X[p]*Weight[3]*MASSE[p]; 
            veloy[l] += (double)VELOCITY_Y[p]*Weight[3]*MASSE[p];
            veloz[l] += (double)VELOCITY_Z[p]*Weight[3]*MASSE[p];
            gwichterl[l] += MASSE[p]*Weight[3];
 
            l = o+i+nx*(j+1)+nx*ny*(k+1);
            getWeights( l, p, Weight, sx);  /* calculate weights*/
            Rho[l] += (double)((MASSE[p] * Weight[3])/(sx*sy*sz)); /* calculate density*/
            velox[l] += (double)VELOCITY_X[p]*Weight[3]*MASSE[p]; 
            veloy[l] += (double)VELOCITY_Y[p]*Weight[3]*MASSE[p];
            veloz[l] += (double)VELOCITY_Z[p]*Weight[3]*MASSE[p];
            gwichterl[l] += MASSE[p]*Weight[3];

            l = o+(i+1)+nx*(j+1)+ nx*ny*(k+1);
            getWeights( l, p, Weight, sx);  /* calculate weights*/
            Rho[l] += (double)((MASSE[p] * Weight[3])/(sx*sy*sz)); /* calculate density*/
            velox[l] += (double)VELOCITY_X[p]*Weight[3]*MASSE[p]; 
            veloy[l] += (double)VELOCITY_Y[p]*Weight[3]*MASSE[p];
            veloz[l] += (double)VELOCITY_Z[p]*Weight[3]*MASSE[p];
            gwichterl[l] += MASSE[p]*Weight[3];
        }
        sx/=2.0;     /*half the spacing*/
        sy/=2.0;
        sz/=2.0;
      }

   for (l=0; l<nmp; l++)
      {
       if (gwichterl[l] == 0)  gwichterl[l]=1E-5;
       
       velox[l] /= (double)gwichterl[l];
       veloy[l] /= (double)gwichterl[l];
       veloz[l] /= (double)gwichterl[l];
 
      }
    

    /*schreiben der velocity arrays*/

    char filename_v[80];
    char filename_d[80];
    FILE *fp;
    
    for (n=0; n<nlevels; n++){

    sprintf(filename_v, "velo_level_%i_time_%i.bin", n, erwin);

    fp=fopen(filename_v, "wb");    
   
    o = n*nx*ny*nz;
    for (k=0;k<nz;k++)            
            for(j=0;j<ny;j++)
                for(i=0;i<nx;i++){
                l = o+i+nx*j+nx*ny*k;
                velox[l]=velox[l]*(double)sqrt(kkanow);
                veloy[l]=veloy[l]*(double)sqrt(kkanow);
                veloz[l]=veloz[l]*(double)sqrt(kkanow);
                //printf("velox:%g   kkanow:%f   posgitter[l][0]:%f, header1.HubbleParam:%g\n", velox[l], kkanow, posgitter[l][0], header1.HubbleParam);
   	        fwrite(&velox[l], sizeof(double), 1, fp);
                fwrite(&veloy[l], sizeof(double), 1, fp);
                fwrite(&veloz[l], sizeof(double), 1, fp);	
   	     }
    
    fclose(fp);

    sprintf(filename_d, "nbody_rho_level_%i_time_%i.bin", n, erwin);         

    fp = fopen(filename_d,"wb");

    o = n*nx*ny*nz;
   	for (k=0;k<nz;k++)            /*set density field to array*/
            for(j=0;j<ny;j++)
                for(i=0;i<nx;i++){
   	    	l = o+i+nx*j+nx*ny*k;
		otto_tmp[l]=Rho[l];
		fwrite(&otto_tmp[l], sizeof(double), 1, fp);
				
   	     }
    fclose(fp);  
    }

    free(Weight);
    free(velox);
    free(veloy);
    free(veloz);
    free(gwichterl);
    free(otto_tmp);
   
    return;
}

/*------------------------------------------------------------------------------------------*/

void getWeights( int m, int p, float *Weight, float sx)
    {
     
     float d[3], a1, a2, a3;
     int t; 
               
     a1 = (POSITION_X[p] - posgitter[m][0]);
     a2 = (POSITION_Y[p] - posgitter[m][1]);
     a3 = (POSITION_Z[p] - posgitter[m][2]);
     d[0] = fabs(a1);
     d[1] = fabs(a2);
     d[2] = fabs(a3);
     
     for (t=0; t<3; t++)
      {
      Weight[t] = (1.0 - (d[t]/sx)); /*assumtion: sx=sy=sz*/      
      }
    
       
     Weight[3] = Weight[0] * Weight[1] * Weight[2]; /* calculate total Weight*/
     if(Weight[3] < 0) Weight[3] = 0.0;   
     return;
}





/*------------------------------------------------------------------------------------------*/

     
