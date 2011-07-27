/*K&K potnester routine*/

int dichterechner(int);
int unit_conversion(void);
int load_snapshot(int);
int allocate_memory_leser(void);
int reordering_leser(void);
/*K&K Potreader*/
void CICmethod(int);                     /* CIC method               */
void getWeights(int, int, float *,float);	/* calculates mass assignement function for CIC*/
void getPoten(void);    				/* calculates potential on nested grid*/ 
void FFTW_rho(int, int);
void FFTW_dist(void);
void pot_slab_collector(int);
void getbackPoten(int, int);


void speicher_pot_fits(int, int);
