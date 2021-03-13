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

void source(int *ipt, double *x, double *y, double *z, double *dx, double *dy, double *dz, double *e, double *we, double *param){

/************************************************* Input *****************************************************/
	
	#define ord 3
	char pt = 'n';
    char* filenames[ord] = {"../Decay/Fe59.csv", "../2_bunker_n_bw0/activ_fe.ssv", NULL};
    ReadFun readfuns[ord] = {Decay_read, SSVtally_read, NULL};
    double trasl_plist[3] = {0, 0, 0};
    double rot_plist[3] = {0, 0, 0};
    int switch_x2z = 0;

    int dims[3] = {1, 3, 3};
    double bw_E[1] = {0};
    double bw_pos[3] = {1.831, 2.362, 126.326};
    double bw_dir[3] = {INFINITY,INFINITY,INFINITY};
    char* bwfilename = NULL; // "../1_guia_n_knn/D_tracks_bw_knn.txt"; // 
    int variable_bw = 0;
    PerturbFun perturb[] = {Let_perturb, Vol_perturb, Isotrop_perturb};
    int n_gp[3] = {0, 6, 0};
    double gp_E[] = {};
    double gp_pos[] = {-10,30,-15,15,0,2360}; // {7, 20, -92900};
    double gp_dir[] = {};
    double trasl_metric[3] = {0, 0, 0};
    double rot_metric[3] = {0, 0, 0};

/*********************************************** Fin Input ***************************************************/

   // *********************************** Declaracion variables globales **************************************

	static long int N_simul;

    static PList* plist;
    static MetricSepVar* metric;
	static KSource *ksource;

	// **************************************** Inicializacion ************************************************

	int i;
	if(initialized == 0){
		printf("\nCargando fuentes...  ");
		
		plist = PList_create(pt, ord, filenames, readfuns, trasl_plist, rot_plist, switch_x2z);
    	Metric* metrics[3] = {Metric_create(dims[0], bw_E, perturb[0], n_gp[0], gp_E),
    	                      Metric_create(dims[1], bw_pos, perturb[1], n_gp[1], gp_pos),
    	                      Metric_create(dims[2], bw_dir, perturb[2], n_gp[2], gp_dir)};
    	metric = MetricSepVar_create(3, metrics, bwfilename, variable_bw, trasl_metric, rot_metric);
		ksource = KS_create(1, plist, metric);

		N_simul = param[0];

		srand(time(NULL));

		initialized = 1;
		printf("Hecho\n");
	}

	// ********************************************** Sorteo ***********************************************************

	Part part;
	double w;

	KS_sample(ksource, &pt, &part, &w, 1);

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
	}

	return;
}