#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "ksource.h"


KSource* KS_create(double J, PList* plist, Geometry* geom){
	KSource* ks = (KSource*)malloc(sizeof(KSource));
	ks->J = J;
	ks->plist = plist;
	ks->geom = geom;
	return ks;
}

KSource* KS_open(const char* filename, int bw_null){
	int i, j;
	char buffer[LINE_MAX_LEN], *pbuffer;
	double J;

	char pt;	
	char mcplfile[NAME_MAX_LEN];
	double *trasl_plist=NULL, *rot_plist=NULL;

	int ord_geom;
	double *trasl_metric=NULL, *rot_metric=NULL;
	int switch_x2z, variable_bw;
	char* bwfilename=NULL;

	// Leer archivo
	FILE* file;
	if((file=fopen(filename, "r")) == 0){
		printf("Error en KS_open: No se pudo abrir archivo %s\n", filename);
		return NULL;
	}
	fgets(buffer, LINE_MAX_LEN, file); // # J
	if(strcmp(buffer, "# J [1/s]:\n") != 0){
		printf("Error en KS_open: Formato de archivo de fuente %s invalido\n", filename);
		return NULL;
	}
	fscanf(file, "%le\n", &J); // Leer J
	// PList
	fgets(buffer, LINE_MAX_LEN, file); // # PList
	fscanf(file, "%c\n", &pt); // Leer pt
	fgets(mcplfile, NAME_MAX_LEN, file); // Leer nombre de archivo mcpl
	mcplfile[strcspn(mcplfile, "\n")] = 0;
	fgets(buffer, LINE_MAX_LEN, file); // Leer traslacion de PList
	if(strlen(buffer) > 1){
		trasl_plist = (double*)malloc(3 * sizeof(double));
		sscanf(buffer, "%lf %lf %lf", &trasl_plist[0], &trasl_plist[1], &trasl_plist[2]);
	}
	fgets(buffer, LINE_MAX_LEN, file); // Leer rotacion de PList
	if(strlen(buffer) > 1){
		rot_plist = (double*)malloc(3 * sizeof(double));
		sscanf(buffer, "%lf %lf %lf", &rot_plist[0], &rot_plist[1], &rot_plist[2]);
	}
	fscanf(file, "%d\n", &switch_x2z); // Leer switch_x2z
	// Metric
	fgets(buffer, LINE_MAX_LEN, file); // # Metric
	fscanf(file, "%d\n", &ord_geom); // Leer ord_geom
	int dims[ord_geom], ngps[ord_geom];
	char metricnames[ord_geom][NAME_MAX_LEN];
	double *gps[ord_geom];
	double *bws[ord_geom];
	PerturbFun perturbs[ord_geom];
	for(i=0; i<ord_geom; i++){
		fgets(metricnames[i], NAME_MAX_LEN, file); // Leer metricname
		fscanf(file, "%d\n", &dims[i]); // Leer dim
		fscanf(file, "%d", &ngps[i]); // Leer ngp
		gps[i] = (double*)malloc(ngps[i]*sizeof(double));
		for(j=0; j<ngps[i]; j++) fscanf(file, "%lf", &gps[i][j]); // Leer gp
		fgets(buffer, LINE_MAX_LEN, file);
	}
	fgets(buffer, LINE_MAX_LEN, file); // Leer traslacion de Metric
	if(strlen(buffer) > 1){
		trasl_metric = (double*)malloc(3 * sizeof(double));
		sscanf(buffer, "%lf %lf %lf", &trasl_metric[0], &trasl_metric[1], &trasl_metric[2]);
	}
	fgets(buffer, LINE_MAX_LEN, file); // Leer rotacion de Metric
	if(strlen(buffer) > 1){
		rot_metric = (double*)malloc(3 * sizeof(double));
		sscanf(buffer, "%lf %lf %lf", &rot_metric[0], &rot_metric[1], &rot_metric[2]);
	}
	fscanf(file, "%d\n", &variable_bw); // Leer variable_bw
	if(bw_null) variable_bw = 0;
	if(variable_bw){
		bwfilename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
		fgets(bwfilename, NAME_MAX_LEN, file);
		bwfilename[strcspn(bwfilename, "\n")] = 0;
		for(i=0; i<ord_geom; i++) bws[i] = NULL;
	}
	else{
		for(i=0; i<ord_geom; i++){
			bws[i] = (double*)malloc(dims[i]*sizeof(double));
			for(j=0; j<dims[i]; j++) bw_null ? bws[i][j]=0 : fscanf(file, "%lf", &bws[i][j]);
		}
	}
	// Crear PList
	PList* plist = PList_create(pt, mcplfile, trasl_plist, rot_plist, switch_x2z);
	// Crear Metric
	for(i=0; i<ord_geom; i++){
		if(strcmp(metricnames[i], "Energy\n") == 0) perturbs[i] = E_perturb;
		else if(strcmp(metricnames[i], "Lethargy\n") == 0) perturbs[i] = Let_perturb;
		else if(strcmp(metricnames[i], "SurfXY\n") == 0) perturbs[i] = SurfXY_perturb;
		else if(strcmp(metricnames[i], "Vol\n") == 0) perturbs[i] = Vol_perturb;
		else if(strcmp(metricnames[i], "Guide\n") == 0) perturbs[i] = Guide_perturb;
		else if(strcmp(metricnames[i], "Polar\n") == 0) perturbs[i] = Polar_perturb;
		else if(strcmp(metricnames[i], "Isotrop\n") == 0) perturbs[i] = Isotrop_perturb;
		else{
			printf("Error: Metrica %s invalida\n", metricnames[i]);
			return NULL;
		}
	}
	Metric* metrics[ord_geom];
	for(i=0; i<ord_geom; i++) metrics[i] = Metric_create(dims[i], bws[i], perturbs[i], ngps[i], gps[i]);
	Geometry* metric = Geom_create(ord_geom, metrics, bwfilename, variable_bw, trasl_metric, rot_metric);
	// Crear KSource
	KSource* s = KS_create(J, plist, metric);
	// Liberar variables alocadas
	free(trasl_plist); free(rot_plist);
	free(trasl_metric); free(rot_metric);
	free(bwfilename);
	for(i=0; i<ord_geom; i++){ free(gps[i]); free(bws[i]); }
	fclose(file);

	return s;
}

int KS_sample(KSource* ks, mcpl_particle_t* part, double w_crit, WeightFun bias){
	if(w_crit <= 0){
		PList_get(ks->plist, part);
		PList_next(ks->plist);
		Geom_next(ks->geom);
	}
	else{ // Normalizo w a 1
		double bs;
		int resamples = 0;
		while(1){
			PList_get(ks->plist, part);
			if(bias) bs = bias(part);
			else bs = 1;
			if(part->weight*bs > w_crit){ // Si w*bs>w_crit, uso w_crit/w*bs como prob de avanzar en la lista
				if(rand() < w_crit/(part->weight*bs)*RAND_MAX){
					PList_next(ks->plist);
					Geom_next(ks->geom);
				}
				break;
			}
			else{ // Si w*bs<w_crit, uso w*bs/w_crit como prob de tomar la particula
				PList_next(ks->plist);
				Geom_next(ks->geom);
				if(rand() < (part->weight*bs)/w_crit*RAND_MAX) break;
			}
			if(resamples++ > MAX_RESAMPLES) return 1;
		}
		part->weight = 1/bs;
	}
	Geom_perturb(ks->geom, part);
	return 0;
}

double KS_w_mean(KSource* ks, int N, WeightFun bias){
	int i;
	char pt;
	mcpl_particle_t part;
	double w_mean=0;
	for(i=0; i<N; i++){
		KS_sample(ks, &part, -1, NULL);
		if(bias) w_mean += part.weight * bias(&part);
		else w_mean += part.weight;
	}
	return w_mean / N;
}

void KS_destroy(KSource* ks){
	PList_destroy(ks->plist);
	Geom_destroy(ks->geom);
	free(ks);
}


MultiSource* MS_create(int len, KSource** s, const double* ws){
	MultiSource* ms = (MultiSource*)malloc(sizeof(MultiSource));
	ms->len = len;
	ms->s = (KSource**)malloc(len*sizeof(KSource*));
	ms->ws = (double*)malloc(len*sizeof(double));
	ms->J = 0;
	int i;
	for(i=0; i<len; i++){
		ms->s[i] = s[i];
		ms->ws[i] = ws[i];
		ms->J += s[i]->J;
	}
	ms->cdf = (double*)malloc(ms->len*sizeof(double));
	for(i=0; i<ms->len; i++) ms->cdf[i] = ms->ws[i];
	for(i=1; i<ms->len; i++) ms->cdf[i] += ms->cdf[i-1];
	return ms;
}

MultiSource* MS_open(int len, const char** filenames, const double* ws, int bw_null){
	KSource* s[len];
	int i;
	for(i=0; i<len; i++) s[i] = KS_open(filenames[i], bw_null);
	return MS_create(len, s, ws);
}

int MS_sample(MultiSource* ms, mcpl_particle_t* part, double w_crit, WeightFun bias){
	double y = rand() / ((double)RAND_MAX+1);
	int i, ret;
	if(ms->cdf[ms->len-1] <= 0) i = (int)(y*ms->len);
	else for(i=0; y*ms->cdf[ms->len-1]>ms->cdf[i]; i++);
	ret = KS_sample(ms->s[i], part, w_crit, bias);
	if(ms->cdf[ms->len-1] > 0) part->weight *= (ms->s[i]->J/ms->J) / (ms->ws[i]/ms->cdf[ms->len-1]);
	else part->weight *= (ms->s[i]->J/ms->J) * ms->len;
	return ret;
}

double MS_w_mean(MultiSource* ms, int N, WeightFun bias){
	double w_mean=0;
	int i;
	for(i=0; i<ms->len; i++) w_mean += KS_w_mean(ms->s[i], N, bias);
	return w_mean / ms->len;
}

void MS_destroy(MultiSource* ms){
	int i;
	for(i=0; i<ms->len; i++) KS_destroy(ms->s[i]);
	free(ms->s);
	free(ms->ws);
	free(ms->cdf);
	free(ms);
}
