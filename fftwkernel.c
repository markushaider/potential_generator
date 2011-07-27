#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "variablen.h"
#include "proto.h"
#include <fftw3.h>

/*------------------------------------------------------------------------------------------*/
void FFTW_dist(void)
{
        
    fftw_plan plan;
    
    int pnx=2*ZA, pny=2*ZA, pnz=2*ZA;
                
    /* next create the plans for FFTW*/
    
    /* create the forward and backward plans: */
    plan = fftw_plan_dft_r2c_3d(pnx, pny, pnz, dist, distft, FFTW_ESTIMATE);
      
          
    /* Now, compute the forward transform: */
    fftw_execute(plan);
    fftw_destroy_plan(plan);
               
    return;
}
