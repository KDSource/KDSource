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
	
	static char pt = 'n';
    static double trasl_plist[3] = {0, 0, 300};
    static double rot_plist[3] = {0, 0, 0};
    static int switch_x2z = 0;
    static char* filename = "../1_guia_n_mlcv/D_tracks.ssv";
    static ReadFun readfun = SSV_read;

    static double trasl_metric[3] = {0, 0, 300};
    static double rot_metric[3] = {0, 0, 0};
    static double gp_E[] = {};
    static double gp_pos[] = {7, 20, 92900};
    static double gp_dir[] = {};
    static double bw_E[1] = {0};
    static double bw_pos[2] = {0,0};
    static double bw_dir[3] = {0,0,0};
    static char* bwfilename = NULL;
    static int variable_bw = 0;
    static PerturbFun perturb[] = {Let_perturb, Guide_perturb, Isotrop_perturb};

/*********************************************** Fin Input ***************************************************/

   // *********************************** Declaracion variables globales **************************************

	static Part part;
	static double w;
	static long int N_simul;

    static int dims[3] = {1, 2, 3};
    static double* bw[3] = {bw_E, bw_pos, bw_dir};
    static double* geom_par[3] = {gp_E, gp_pos, gp_dir};

    static PList* plist;
    static Metric* metric;
	static KSource *ksource;

	// **************************************** Inicializacion ************************************************

	int i;
	if(initialized == 0){
		printf("\nCargando fuentes...  ");
		
		plist = PListSimple_create(pt, trasl_plist, rot_plist, switch_x2z, filename, readfun);
		metric = MetricSepVar_create(dims, bw, bwfilename, variable_bw, perturb, trasl_metric, rot_metric, geom_par);
		ksource = KS_create(1, plist, metric);

		N_simul = param[0];

		srand(time(NULL));

		initialized = 1;
		printf("Hecho\n");
	}

	// ********************************************** Sorteo ***********************************************************

	KS_sample(ksource, &pt, &part, &w, 1);

	if(pt == 'n') *ipt = 1;
	else if(pt == 'p') *ipt = 2;
	else *ipt = 1;
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