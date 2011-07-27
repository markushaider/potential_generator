#include <stdio.h>
#include <complex.h>
#include <fftw3.h>

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7


int TotNumPart;

/*ARRAYS fuer fftw*/
int N, N22;

 /*ARRAYS fuer Dichte und Potentiale*/
float posgitter[ZA*ZA*ZA*4][3];
float *TABLE;


struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

int     NumPart1, Ngas1;

struct particle_data2
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *P2;

int *Id;
double  Time, Redshift;


/*K&K scale factor*/
float a_now;
/*K&K ENDE*/

/*K&K fuer Potentiale */

double Lx, Ly, Lz; /*die Kantenlaenge des GRID1 fuer Hydro*/
float *POSITION_X, *POSITION_Y, *POSITION_Z;
float *VELOCITY_X, *VELOCITY_Y, *VELOCITY_Z;

double *Rho, *Phi;

/*K&K Arrays fuer fftw*/
double *backpot;
fftw_complex *distft, *disttmp;     
fftw_complex *rhoft;	 
double *rho;
double *dist;

float kkanow;

/*velocity arrays*/

double *velox, *veloy, *veloz;
float *gwichterl;
   
float *MASSE;

/*#############################################################################################*/
int ThisTask;		/*!< the number of the local processor  */
int NTask;               /*!< number of processors */
int PTask;	        /*!< note: NTask = 2^PTask */

