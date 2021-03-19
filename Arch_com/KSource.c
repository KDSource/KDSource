#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<time.h>

#include "ksource.h"
#include "metrics.h"
#include "plists.h"
#include "aux.h"


static int initialized = 0;
static long int N = 0;
static double t_sample = 0;

void source(int *ipt, double *x, double *y, double *z, double *dx, double *dy, double *dz, double *e, double *we, double *param){

	clock_t start = clock();

/************************************************* Input *****************************************************/
	
	char filename[] = "/home/inti/Documents/Maestria/Simulaciones/1_guia_n_knn/D_tracks_source.txt";

/*********************************************** Fin Input ***************************************************/

   // *********************************** Declaracion variables globales **************************************

	static long int N_simul;

	static KSource *ksource;

	// **************************************** Inicializacion ************************************************

	int i;
	if(initialized == 0){
		printf("\nCargando fuentes...  ");

		ksource = KS_open(filename);

		N_simul = param[0];

		srand(time(NULL));

		initialized = 1;
		printf("Hecho\n");
	}

	// ********************************************** Sorteo ***********************************************************

	Part part;
	double w;
	char pt;

	KS_sample(ksource, &pt, &part, &w, 1e-2, NULL);

	if(pt == 'n') *ipt = 1;
	else if(pt == 'p') *ipt = 2;
	else{
		printf("Error: Particula no reconocida. Se tomara como neutron.\n");
		*ipt = 1;
	}
	*x = part.pos[0];
	*y = part.pos[1];
	*z = part.pos[2];
	*dx = part.dir[0];
	*dy = part.dir[1];
	*dz = part.dir[2];
	*e = part.E;
	*we = w;

	N++;

	// *********************************************** Finalizacion ***************************************************

	if(N%N_simul == 0){
		printf("\nDestruyendo fuentes...  ");
		
		KS_destroy(ksource);

		initialized = 0;
		printf("Hecho\n");
		printf("Tiempo de muestreo: %lf s\n", t_sample);
	}

	clock_t end = clock();
	t_sample += (float)(end - start) / CLOCKS_PER_SEC;

	return;
}