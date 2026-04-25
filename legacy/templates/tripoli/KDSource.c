#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>

#include<math.h>

#include "mcpl.h"

#include "kdsource.h"
#include "metrics.h"
#include "plists.h"
#include "aux.h"

#define EPSILON_GEO 1E-4


static int initialized = 0;
static long int N = 0;
static double I = 0, p2 = 0;
static double t_sample = 0;

void source(int *ipt, double *x, double *y, double *z, double *dx, double *dy, double *dz, double *e, double *we, double *param){

	clock_t start = clock();

/************************************************* Input *****************************************************/
	
	#define len 1                                          // Number of sources
	const char* filenames[len] = {"../simul1/source.xml"}; // Source XML files
	double ws[len] = {1};                                  // Weight of each source
	int use_kde = 1;                                       // Whether to use KDE for sampling
	int loop = 1;                                          // Whether to loop over list when reach end

	WeightFun bias = NULL;                                 // Bias function

/*********************************************** End Input ***************************************************/


   // *********************************** Global variables declaration ***************************************

	static long int N_simul;

	static MultiSource *msource;
	static double w_crit;

	// **************************************** Initialization ************************************************

	int i;
	if(initialized == 0){
		printf("\nLoading sources...\n");

		msource = MS_open(len, filenames, ws);
		w_crit = use_kde ? MS_w_mean(msource, 1000, bias) : -1;

		N_simul = param[0]*param[1];

		srand(time(NULL));

		initialized = 1;
		printf("Done.\n");
	}

	// ********************************************** Sampling ************************************************

	mcpl_particle_t part;
	double w;
	char pt;

	MS_sample2(msource, &part, use_kde, w_crit, bias, loop);

	if(part.pdgcode == 2112) *ipt = 1;
	else if(part.pdgcode == 22) *ipt = 2;
	else if(part.pdgcode == 11) *ipt = 3;
	else if(part.pdgcode == -11) *ipt = 4;
	else{
		printf("Warning: PDG code %d unknown. Taking as neutron.\n", part.pdgcode);
		*ipt = 1;
	}
	*x = part.position[0];
	*y = part.position[1];
	*z = part.position[2];
	*dx = part.direction[0];
	*dy = part.direction[1];
	*dz = part.direction[2];
	*e = part.ekin;
	*we = part.weight;

	*x += *dx * EPSILON_GEO;
	*y += *dy * EPSILON_GEO;
	*z += *dz * EPSILON_GEO;

	N++;
	I += *we;
	p2 += *we**we;

	// *********************************************** Finalizacion ***************************************************

		if(N%N_simul == 0){
		printf("\nDestroying sources...\n");
		
		MS_destroy(msource);

		initialized = 0;
		printf("Done.\n");
		printf("Sampling time: %lf s\n", t_sample);
		printf("Produced particles: I err N %lf %lf %ld\n", I, sqrt(p2), N);
	}

	clock_t end = clock();
	t_sample += (float)(end - start) / CLOCKS_PER_SEC;

	return;
}
