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
	
	char pt = 'n';
    char* filename = "../1_guia_n_mlcv/D_tracks.ssv";
    ReadFun readfun = SSV_read;
    double trasl_plist[3] = {0, 0, 300};
    double rot_plist[3] = {0, 0, 0};
    int switch_x2z = 0;

    int dims[3] = {1, 2, 3};
    double bw_E[1] = {0};
    double bw_pos[2] = {0,0};
    double bw_dir[3] = {0,0,0};
    char* bwfilename = NULL;
    int variable_bw = 0;
    PerturbFun perturb[] = {Let_perturb, Guide_perturb, Isotrop_perturb};
    int n_gp[3] = {0, 3, 0};
    double gp_E[] = {};
    double gp_pos[] = {7, 20, 92900};
    double gp_dir[] = {};
    double trasl_metric[3] = {0, 0, 300};
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
		
		plist = PListSimple_create(pt, filename, readfun, trasl_plist, rot_plist, switch_x2z);
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