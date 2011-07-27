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



void getWeightTSC(int tscWidth, float dx, float * weight) {

  int ii;
  /* boundary cells of triangle */
  weight[0]          = 0.5*(0.5-dx)*(0.5-dx)/(tscWidth*tscWidth);
  weight[2*tscWidth] = 0.5*(0.5+dx)*(0.5+dx)/(tscWidth*tscWidth);
  /* central cell */
  weight[tscWidth]=-dx*dx/(tscWidth*tscWidth)-1./(4.*tscWidth*tscWidth)+1./tscWidth;
  /* left side */
  for (ii=1; ii!=tscWidth; ii++) {
    weight[ii]=1.0/(tscWidth*tscWidth)*(ii-tscWidth-dx)+1.0/tscWidth;
  }
  /* right side */
  for (ii=1; ii!=tscWidth; ii++) {
    weight[ii+tscWidth]=-1.0/(tscWidth*tscWidth)*(ii-dx)+1.0/tscWidth;
  }
  /* Check the weights */
  //printf("Check weights\n");
  float weightIntegral=0;
  for (ii=0; ii!=2*tscWidth+1; ii++) {
    //printf("ii: %i weight: %g\n",ii,weight[ii]);
    weightIntegral+=weight[ii];
  }
  //printf("\n");
  if(weightIntegral<0.999 || weightIntegral>1.0001) {
    printf("Achtung Achtung. Da hilft nur noch Hubshraubereinsatz\n");
    printf("Weight Integral: %g\n",weightIntegral);
}
 

}


/*------------------------------------------------------------------------------------------*/
void CICmethod(int erwin) {
    
  /*-------------------------------------------*/
  /*erzeugen der Gittereigenschaften*/
    
  int n, m, l, nlevels=4, k, j, i;
  int nx, ny, nz;
  float Lx, Ly, Lz;
  float sx, sy, sz;
  double *otto_tmp;
  kkanow = (float) header1.time;
  printf("kkanow: %f, header1.time: %g\n", kkanow, header1.time);
  Lx=Ly=Lz=HYDROBOX; //defined in makefile, currently 20000 kpc
    
  nx=ny=nz=ZA; // defined in Makefile, = resolution, e.g. 128
        
  sx=(Lx/(float)nx); // length of cell in kpc
  sy=(Ly/(float)ny);
  sz=(Lz/(float)nz);
  velox= (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
  veloy= (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
  veloz= (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
  float * massWeight = (float*) malloc(ZA*ZA*ZA*nlevels*sizeof(float));
    
  // create array of the physical cell positions in kpc
  for (n=0;n<nlevels;n++) {
    m = n*nx*ny*nz;
    for (k=0;k<nz;k++) {
      for (j=0;j<ny;j++) {
	for (i=0;i<nx;i++) {
	  l=m+i+nx*j+nx*ny*k;
	  // get the position of the center of a cell
	  /* in a centered cube */
	  posgitter[l][0] = sx*(((float)i)-0.5*(float)(nx-1));
	  posgitter[l][1] = sy*(((float)j)-0.5*(float)(ny-1));
          posgitter[l][2] = sz*(((float)k)-0.5*(float)(nz-1));
	}
      } 
    }
    sx/=2.0; // half the length for the next inner grid
    sy/=2.0;
    sz/=2.0;
  }
  /*-------------------------------------------*/
    
  int o,p;
    
  // old
  float *Weight;
  Weight = (float*) malloc((3 + 1) * sizeof(float));
  int nmp=ZA*ZA*ZA*nlevels;
    
  // set again to cell-length of outermost grid
  sx=(Lx/(float)nx);
  sy=(Ly/(float)ny);
  sz=(Lz/(float)nz);
  
  // 1 array for all for grids
  Rho = (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
  otto_tmp = (double*) malloc(ZA*ZA*ZA*nlevels*sizeof(double));
  
  for (m = 0; m < nmp; m++) {
    Rho[m] = 0.0;  /* set potential to zero*/
    velox[m] = 0.0;  /* set to zero*/
    veloy[m] = 0.0;  /* set to zero*/
    veloz[m] = 0.0;  /* set to zero*/
    massWeight[m] = 0.0;  /* set to zero*/
  }

  /* Compute Density for outermost grid */
  /* Loop over all particles */
  for (p=0; p<TotNumPart;p++) {

    if(p%10000==0) {
      printf("p: %i\n",p);
    }

    if(MASSE[p] == 0.0) {
      printf("Dead particle :-(");
      continue;
    }
    
    /* Get the index of the cell in which the particle lies */
    /* Position is running from -16000/h to +16000/h */
    i = (int)(floor(POSITION_X[p]/sx) + nx/2);
    j = (int)(floor(POSITION_Y[p]/sy) + ny/2);
    k = (int)(floor(POSITION_Z[p]/sz) + nz/2);
    //printf("p: %i, i,j,k: %i %i %i\n",p,i,j,k);
    //k = (int)(floor(POSITION_Z[p]/sz) + 0.5*(float)(nz-1));
    //printf("position X: %g sx: %g i: %i\n", POSITION_X[p], sx, i);
      
    /* Are we within our Hydro-Grid? */
    /* If not, continue to next particle */
    if( i < 0 || i > (nx-1)) {
      continue;
    }
    else if ( j < 0 || j > (ny-1)) {
      continue;
    }
    else if ( k < 0 || k > (nz-1)) {
      continue;
    }



    /* Loop over cells in outer grid, which are affected by the TSC of a particle */
    int tscWidth = 5; //in cell widths, arbitrary chosen, dependent on resolution


    /* Calculation of weights */
    int ii,jj,kk;
    /* dx,dy,dz are the distance from the cell center */
    /* it is important for calculating the weights*/
    float dx,dy,dz;
    dx=(POSITION_X[p]-posgitter[i+nx*j+nx*ny*k][0])/sx;
    //printf("Pos: %g posgitter: %g sx: %g dx: %g \n",POSITION_X[p],posgitter[i+nx*j+nx*ny*k][0],sx,dx);
    dy=(POSITION_Y[p]-posgitter[i+nx*j+nx*ny*k][1])/sx;
    dz=(POSITION_Z[p]-posgitter[i+nx*j+nx*ny*k][2])/sx;
    if(fabs(dx)>0.5 || fabs(dy)>0.5 || fabs(dz) >0.5) {
      printf("Attention, something is wrong:\n");
      printf("dx = %g, dy = %g, dz = %g \n",dx,dy,dz);
    }

    float weightX[2*tscWidth+1];
    float weightY[2*tscWidth+1];
    float weightZ[2*tscWidth+1];
    getWeightTSC(tscWidth,dx,weightX);
    getWeightTSC(tscWidth,dy,weightY);
    getWeightTSC(tscWidth,dz,weightZ);

    /* Filling density array */
    ii=0;
    int counter=0;
    for (ii=-tscWidth; ii<=tscWidth; ii++) {
      for (jj=-tscWidth; jj<=tscWidth; jj++) {
	for (kk=-tscWidth; kk<=tscWidth; kk++) {
	  l=(i+ii+nx)%ZA+nx*((j+jj+ny)%ZA)+nx*ny*((k+kk+nz)%ZA);
	  counter+=1;
	  //printf("%i laufindex: %i\n",counter, l);
	  Rho[l]+=(double)(MASSE[p]*weightX[ii+tscWidth]*weightY[jj+tscWidth]*weightZ[tscWidth+kk]/(sx*sy*sz));
	  velox[l]+=(double)(VELOCITY_X[p]*weightX[ii+tscWidth]*weightY[jj+tscWidth]*weightZ[tscWidth+kk]*MASSE[p]);
	  veloy[l]+=(double)(VELOCITY_Y[p]*weightX[ii+tscWidth]*weightY[jj+tscWidth]*weightZ[tscWidth+kk]*MASSE[p]);
	  veloz[l]+=(double)(VELOCITY_Z[p]*weightX[ii+tscWidth]*weightY[jj+tscWidth]*weightZ[tscWidth+kk]*MASSE[p]);
	  massWeight[l]+= MASSE[p]*weightX[ii+tscWidth]*weightY[jj+tscWidth]*weightZ[tscWidth+kk];
	}
      }
    }
    fflush(stdout);
  
  } //end of particle loop


  for (l=0; l<nx*ny*nz;l++) {
    if (massWeight[l] == 0.0) {
      printf("Achtung, massWeight: %g, rho: %g\n",massWeight[l],Rho[l]);
      massWeight[l]=1E-5;
    }
    /* Divide the velocity by the total weights */
    velox[l] /= (double)massWeight[l];
    veloy[l] /= (double)massWeight[l];
    veloz[l] /= (double)massWeight[l];
  }

  /* Fill smaller grids */
  int ii,jj,kk;
  int indX,indY,indZ;
  for(n=1; n<4; n++) {
    o = n*nx*ny*nz;             
    for(ii = nx/4; ii < 3*nx/4; ii++) {
      for(jj = ny/4; jj < 3*ny/4; jj++) {
	for(kk = nz/4; kk < 3*nz/4; kk++) {

	  l = (n-1)*nx*ny*nz+ii+nx*jj+nx*ny*kk;
	  //ll=o+ii+1+nx*jj+nx*ny*kk;
	  indX=(ii-nx/4)*2;
	  indY=(jj-ny/4)*2;
	  indZ=(kk-nz/4)*2;
	  Rho[o+(indX  )+nx*(indY  )+nx*ny*(indZ  )] = Rho[l];
	  Rho[o+(indX+1)+nx*(indY  )+nx*ny*(indZ  )] = Rho[l];
	  Rho[o+(indX+1)+nx*(indY+1)+nx*ny*(indZ  )] = Rho[l];
	  Rho[o+(indX+1)+nx*(indY+1)+nx*ny*(indZ+1)] = Rho[l];
	  Rho[o+(indX+1)+nx*(indY  )+nx*ny*(indZ+1)] = Rho[l];
	  Rho[o+(indX  )+nx*(indY+1)+nx*ny*(indZ  )] = Rho[l];
	  Rho[o+(indX  )+nx*(indY+1)+nx*ny*(indZ+1)] = Rho[l];
	  Rho[o+(indX  )+nx*(indY  )+nx*ny*(indZ+1)] = Rho[l];

	  velox[o+(indX  )+nx*(indY  )+nx*ny*(indZ  )] = velox[l];
	  velox[o+(indX+1)+nx*(indY  )+nx*ny*(indZ  )] = velox[l];
	  velox[o+(indX+1)+nx*(indY+1)+nx*ny*(indZ  )] = velox[l];
	  velox[o+(indX+1)+nx*(indY+1)+nx*ny*(indZ+1)] = velox[l];
	  velox[o+(indX+1)+nx*(indY  )+nx*ny*(indZ+1)] = velox[l];
	  velox[o+(indX  )+nx*(indY+1)+nx*ny*(indZ  )] = velox[l];
	  velox[o+(indX  )+nx*(indY+1)+nx*ny*(indZ+1)] = velox[l];
	  velox[o+(indX  )+nx*(indY  )+nx*ny*(indZ+1)] = velox[l];

	  veloy[o+(indX  )+nx*(indY  )+nx*ny*(indZ  )] = veloy[l];
	  veloy[o+(indX+1)+nx*(indY  )+nx*ny*(indZ  )] = veloy[l];
	  veloy[o+(indX+1)+nx*(indY+1)+nx*ny*(indZ  )] = veloy[l];
	  veloy[o+(indX+1)+nx*(indY+1)+nx*ny*(indZ+1)] = veloy[l];
	  veloy[o+(indX+1)+nx*(indY  )+nx*ny*(indZ+1)] = veloy[l];
	  veloy[o+(indX  )+nx*(indY+1)+nx*ny*(indZ  )] = veloy[l];
	  veloy[o+(indX  )+nx*(indY+1)+nx*ny*(indZ+1)] = veloy[l];
	  veloy[o+(indX  )+nx*(indY  )+nx*ny*(indZ+1)] = veloy[l];

	  veloz[o+(indX  )+nx*(indY  )+nx*ny*(indZ  )] = veloz[l];
	  veloz[o+(indX+1)+nx*(indY  )+nx*ny*(indZ  )] = veloz[l];
	  veloz[o+(indX+1)+nx*(indY+1)+nx*ny*(indZ  )] = veloz[l];
	  veloz[o+(indX+1)+nx*(indY+1)+nx*ny*(indZ+1)] = veloz[l];
	  veloz[o+(indX+1)+nx*(indY  )+nx*ny*(indZ+1)] = veloz[l];
	  veloz[o+(indX  )+nx*(indY+1)+nx*ny*(indZ  )] = veloz[l];
	  veloz[o+(indX  )+nx*(indY+1)+nx*ny*(indZ+1)] = veloz[l];
	  veloz[o+(indX  )+nx*(indY  )+nx*ny*(indZ+1)] = veloz[l];

	}
      }
    }
  }


  

  /*schreiben der velocity arrays*/

  char filename_v[80];
  char filename_d[80];
  FILE *fp;
  extern char* spath;
  for (n=0; n<nlevels; n++){

    sprintf(filename_v, "%svelo_level_%i_time_%i.bin",spath, n, erwin);

    fp=fopen(filename_v, "wb");    
   
    o = n*nx*ny*nz;
    for (k=0;k<nz;k++) {
      for(j=0;j<ny;j++) {
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
      }
    }   
    fclose(fp);

    sprintf(filename_d, "%snbody_rho_level_%i_time_%i.bin",spath, n, erwin);         

    fp = fopen(filename_d,"wb");

    o = n*nx*ny*nz;
    for (k=0;k<nz;k++)            /*set density field to array*/
      for(j=0;j<ny;j++)
	for(i=0;i<nx;i++) {
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

     
