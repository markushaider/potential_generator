#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "variablen.h"
#include "proto.h"


int load_snapshot(int kk)
{
  FILE *fd;
  char   buf[200];
  int    i,k,dummy,ntot_withmasses, files=1;
  int    n,pc,pc_new,pc_sph;
  
  int kkroller=0, kkroller1=0;

  #define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if (kk < 10) sprintf(buf, "snapshot_00%i", kk);
      if (kk > 9 && kk < 100) sprintf(buf, "snapshot_0%i", kk);
      if (kk > 99) sprintf(buf, "snapshot_%i", kk);

      if(!(fd=fopen(buf,"r")))
    {
      printf("can't open file `%s`\n",buf);
      exit(0);
    }

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      for(k=0, ntot_withmasses=0; k<5; k++)
    {
      if(header1.mass[k]==0)
        ntot_withmasses+= header1.npart[k];
    }

      allocate_memory_leser();
   
      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<header1.npart[k];n++)
        {
          fread(&P2[pc_new].Pos[0], sizeof(float), 3, fd);
          POSITION_X[kkroller]=P2[pc_new].Pos[0]/0.71-(AllBoxSize/(2*0.71)+(float)SHIFT_X*(HYDROBOX/((float)ZA)));  // wir wollen echte kpc, daher division durch 0.71
          //printf("Po_x:%f\n", POSITION_X[kkroller]);
          POSITION_Y[kkroller]=P2[pc_new].Pos[1]/0.71-(AllBoxSize/(2*0.71)+(float)SHIFT_Y*(HYDROBOX/((float)ZA)));
          POSITION_Z[kkroller]=P2[pc_new].Pos[2]/0.71-(AllBoxSize/(2*0.71)+(float)SHIFT_Z*(HYDROBOX/((float)ZA)));
          kkroller++;
          pc_new++;
        }
    }
      SKIP;

      SKIP;
      kkroller=0;
      for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<header1.npart[k];n++)
        {
          fread(&P2[pc_new].Vel[0], sizeof(float), 3, fd);
          VELOCITY_X[kkroller]=P2[pc_new].Vel[0];
          VELOCITY_Y[kkroller]=P2[pc_new].Vel[1];
          VELOCITY_Z[kkroller]=P2[pc_new].Vel[2];
	  kkroller++;
          pc_new++;
        }
    }
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<header1.npart[k];n++)
        {
          fread(&Id[pc_new], sizeof(int), 1, fd);
          pc_new++;
        }
    }
      SKIP;

      if(ntot_withmasses>0)
    SKIP;
      for(k=0, pc_new=pc; k<6; k++)
    {
      for(n=0;n<header1.npart[k];n++)
        {
          P2[pc_new].Type=k;

          if(header1.mass[k]==0)
        fread(&P2[pc_new].Mass, sizeof(float), 1, fd);
          else
        P2[pc_new].Mass= header1.mass[k];
        MASSE[kkroller1]= P2[pc_new].Mass/0.71; //weil GADGET Einheit fuer Masse ist 10E10 M_sun/h
        kkroller1++;
          pc_new++;
        }
    }
      if(ntot_withmasses>0)
    SKIP;
      

      if(header1.npart[0]>0)
    {
      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
        {
          fread(&P2[pc_sph].U, sizeof(float), 1, fd);
          pc_sph++;
        }
      SKIP;

      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
        {
          fread(&P2[pc_sph].Rho, sizeof(float), 1, fd);
          pc_sph++;
        }
      SKIP;

      if(header1.flag_cooling)
        {
          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
        {
          fread(&P2[pc_sph].Ne, sizeof(float), 1, fd);
          pc_sph++;
        }
          SKIP;
        }
      else
        for(n=0, pc_sph=pc; n<header1.npart[0];n++)
          {
        P2[pc_sph].Ne= 1.0;
        pc_sph++;
          }
    }
  

     fclose(fd);
    }

  Time= header1.time;
  Redshift= header1.time;    
  
  return 0;

}
/*------------------------------------------------------------------------------------------*/
/* this routine allocates the memory for the 
 * particle data.
 */
 /*------------------------------------------------------------------------------------------*/
int allocate_memory_leser(void)
{
  printf("allocating memory...\n");

  if(!(P2=malloc(TotNumPart*sizeof(struct particle_data2))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  //P2--;   /* start with offset 1 */

  
  if(!(Id=malloc(TotNumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }

  
  //Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");

  return 0;
}
