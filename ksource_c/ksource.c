#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "ksource.h"


void Part_print(Part* part){
	printf("(%.3le, (%.3lf, %.3lf, %.3lf), (%.3lf, %.3lf, %.3lf))\n",
		part->E, part->pos[0], part->pos[1], part->pos[2], part->dir[0], part->dir[1], part->dir[2]);
}

KSource* KS_create(double J, PList* plist, MetricSepVar* metric){
	KSource* ks = (KSource*)malloc(sizeof(KSource));
	ks->J = J;
	ks->plist = plist;
	ks->metric = metric;
	return ks;
}

KSource* KS_open(char* filename){
	int i;
	char buffer[LINE_MAX_LEN], *pbuffer;
	double J;

	char pt;
	int ord;
	//char *format[ord], *tracksfiles[ord];
	//ReadFun reads[ord];
	double *trasl_plist=NULL, *rot_plist=NULL;

	int dim_E, dim_pos, dim_dir;
	char metricEname[NAME_MAX_LEN], metricposname[NAME_MAX_LEN], metricdirname[NAME_MAX_LEN];
	int ngp_E=0, ngp_pos=0, ngp_dir=0;
	double *gp_E=NULL, *gp_pos=NULL, *gp_dir=NULL;
	double *trasl_metric=NULL, *rot_metric=NULL;
	int switch_x2z, variable_bw;
	double *bw_E=NULL, *bw_pos=NULL, *bw_dir=NULL;
	char* bwfilename=NULL;
	PerturbFun perturb_E, perturb_pos, perturb_dir;

	// Leer archivo
	FILE* file;
	if((file=fopen(filename, "r")) == 0){
		printf("Error en KS_open: No se pudo abrir archivo %s\n", filename);
		return NULL;
	}
	fgets(buffer, LINE_MAX_LEN, file); // # J
	fscanf(file, "%le\n", &J); // Leer J
	// PList
	fgets(buffer, LINE_MAX_LEN, file); // # PList
	fscanf(file, "%c\n", &pt); // Leer pt
	fscanf(file, "%d\n", &ord); // Leer ord
	char formats[ord][NAME_MAX_LEN], *tracksfiles[ord];
	ReadFun reads[ord];
	for(i=0; i<ord; i++) fgets(formats[i], NAME_MAX_LEN, file); // Leer formatos
	for(i=0; i<ord; i++){
		tracksfiles[i] = (char*)malloc(NAME_MAX_LEN*sizeof(char));
		fgets(tracksfiles[i], NAME_MAX_LEN, file); // Leer nombres de archivos de tracks
		tracksfiles[i][strcspn(tracksfiles[i], "\n")] = 0;
	}
	fgets(buffer, LINE_MAX_LEN, file); // Leer traslacion de PList
	if(strlen(buffer)){
		trasl_plist = (double*)malloc(3 * sizeof(double));
		sscanf(buffer, "%lf %lf %lf", &trasl_plist[0], &trasl_plist[1], &trasl_plist[2]);
	}
	fgets(buffer, LINE_MAX_LEN, file); // Leer rotacion de PList
	if(strlen(buffer)){
		rot_plist = (double*)malloc(3 * sizeof(double));
		sscanf(buffer, "%lf %lf %lf", &rot_plist[0], &rot_plist[1], &rot_plist[2]);
	}
	fscanf(file, "%d\n", &switch_x2z); // Leer switch_x2z
	// Metric
	fgets(buffer, LINE_MAX_LEN, file); // # Metric
	fgets(buffer, LINE_MAX_LEN, file); // # E
	fgets(metricEname, NAME_MAX_LEN, file); // Leer metricEname
	fscanf(file, "%d\n", &dim_E); // Leer dim
	fscanf(file, "%d", &ngp_E); // Leer ngp_E
	gp_E = (double*)malloc(ngp_E*sizeof(double));
	if(ngp_E) for(i=0; i<ngp_E; i++) fscanf(file, "%lf\n", &gp_E[i]); // Leer gp_E
	else fgets(buffer, LINE_MAX_LEN, file);
	fgets(buffer, LINE_MAX_LEN, file); // # pos
	fgets(metricposname, NAME_MAX_LEN, file); // Leer metricposname
	fscanf(file, "%d\n", &dim_pos); // Leer dim
	fscanf(file, "%d", &ngp_pos); // Leer ngp_pos
	gp_pos = (double*)malloc(ngp_pos*sizeof(double));
	if(ngp_pos) for(i=0; i<ngp_pos; i++) fscanf(file, "%lf\n", &gp_pos[i]); // Leer gp_pos
	else fgets(buffer, LINE_MAX_LEN, file);
	fgets(buffer, LINE_MAX_LEN, file); // # dir
	fgets(metricdirname, NAME_MAX_LEN, file); // Leer metricdirname
	fscanf(file, "%d\n", &dim_dir); // Leer dim
	fscanf(file, "%d", &ngp_dir); // Leer ngp_dir
	gp_dir = (double*)malloc(ngp_dir*sizeof(double));
	if(ngp_dir) for(i=0; i<ngp_dir; i++) fscanf(file, "%lf\n", &gp_dir[i]); // Leer gp_dir
	else fgets(buffer, LINE_MAX_LEN, file);
	fgets(buffer, LINE_MAX_LEN, file); // Leer traslacion de Metric
	if(strlen(buffer)){
		trasl_metric = (double*)malloc(3 * sizeof(double));
		sscanf(buffer, "%lf %lf %lf", &trasl_metric[0], &trasl_metric[1], &trasl_metric[2]);
	}
	fgets(buffer, LINE_MAX_LEN, file); // Leer rotacion de Metric
	if(strlen(buffer)){
		rot_metric = (double*)malloc(3 * sizeof(double));
		sscanf(buffer, "%lf %lf %lf", &rot_metric[0], &rot_metric[1], &rot_metric[2]);
	}
	fscanf(file, "%d\n", &variable_bw); // Leer variable_bw
	if(variable_bw){
		bwfilename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
		fgets(bwfilename, NAME_MAX_LEN, file);
		bwfilename[strcspn(bwfilename, "\n")] = 0;
	}
	else{
		bw_E = (double*)malloc(dim_E*sizeof(double));
		for(i=0; i<dim_E; i++) fscanf(file, "%lf", &bw_E[i]);
		bw_pos = (double*)malloc(dim_pos*sizeof(double));
		for(i=0; i<dim_pos; i++) fscanf(file, "%lf", &bw_pos[i]);
		bw_dir = (double*)malloc(dim_dir*sizeof(double));
		for(i=0; i<dim_dir; i++) fscanf(file, "%lf", &bw_dir[i]);
	}
	// Crear PList
	for(i=0; i<ord; i++){
		if(strcmp(formats[i], "PTRAC\n") == 0) reads[i] = PTRAC_read;
		else if(strcmp(formats[i], "T4stock\n") == 0) reads[i] = T4stock_read;
		else if(strcmp(formats[i], "SSV\n") == 0) reads[i] = SSV_read;
		else if(strcmp(formats[i], "Decay\n") == 0) reads[i] = Decay_read;
		else if(strcmp(formats[i], "SSVtally\n") == 0) reads[i] = SSVtally_read;
		else if(strcmp(formats[i], "Isotrop\n") == 0) reads[i] = Isotrop_read;
		else{
			printf("Error: Formato %s invalido\n", formats[i]);
			return NULL;
		}
	}
	PList* plist = PList_create(pt, ord, tracksfiles, reads, trasl_plist, rot_plist, switch_x2z);
	// Crear Metric
	if(strcmp(metricEname, "Energy\n") == 0) perturb_E = E_perturb;
	else if(strcmp(metricEname, "Lethargy\n") == 0) perturb_E = Let_perturb;
	else{
		printf("Error: Metrica %s invalida\n", metricEname);
		return NULL;
	}
	if(strcmp(metricposname, "SurfXY\n") == 0) perturb_pos = SurfXY_perturb;
	else if(strcmp(metricposname, "Vol\n") == 0) perturb_pos = Vol_perturb;
	else if(strcmp(metricposname, "Guide\n") == 0) perturb_pos = Guide_perturb;
	else{
		printf("Error: Metrica %s invalida\n", metricposname);
		return NULL;
	}
	if(strcmp(metricdirname, "Polar\n") == 0) perturb_dir = Polar_perturb;
	else if(strcmp(metricdirname, "Isotrop\n") == 0) perturb_dir = Isotrop_perturb;
	else{
		printf("Error: Metrica %s invalida\n", metricdirname);
		return NULL;
	}
	Metric* metric_E = Metric_create(dim_E, bw_E, perturb_E, ngp_E, gp_E);
	Metric* metric_pos = Metric_create(dim_pos, bw_pos, perturb_pos, ngp_pos, gp_pos);
	Metric* metric_dir = Metric_create(dim_dir, bw_dir, perturb_dir, ngp_dir, gp_dir);
	Metric* metrics[3] = {metric_E, metric_pos, metric_dir};
	MetricSepVar* metric = MSV_create(3, metrics, bwfilename, variable_bw, trasl_metric, rot_metric);
	// Crear KSource
	KSource* s = KS_create(J, plist, metric);
	// Liberar variables alocadas
	free(trasl_plist); free(rot_plist);
	free(trasl_metric); free(rot_metric);
	free(gp_E); free(gp_pos); free(gp_dir);
	free(bwfilename);
	free(bw_E); free(bw_pos); free(bw_dir);
	for(i=0; i<ord; i++) free(tracksfiles[i]);
	fclose(file);

	return s;
}

int KS_sample(KSource* ks, char* pt, Part* part, double* w, double w_crit, WeightFun bias){
	*pt = ks->plist->pt;
	if(w_crit <= 0){
		PList_get(ks->plist, part, w);
		PList_next(ks->plist);
		MSV_next(ks->metric);
	}
	else{ // Normalizo w a 1
		double bs;
		int resamples = 0;
		while(1){
			PList_get(ks->plist, part, w);
			if(bias) bs = bias(part);
			else bs = 1;
			if(*w*bs > 1){ // Si w*bs>w_crit, uso w_crit/w*bs como prob de avanzar en la lista
				if(rand() < w_crit/(*w*bs)*RAND_MAX){
					PList_next(ks->plist);
					MSV_next(ks->metric);
				}
				break;
			}
			else{ // Si w*bs<w_crit, uso w*bs/w_crit como prob de tomar la particula
				PList_next(ks->plist);
				MSV_next(ks->metric);
				if(rand() < (*w*bs)/w_crit*RAND_MAX) break;
			}
			if(resamples++ > MAX_RESAMPLES) return 1;
		}
		*w = 1/bs;
	}
	MSV_perturb(ks->metric, part);
	return 0;
}

double KS_w_mean(KSource* ks, int N){
	int i;
	char pt;
	Part part;
	double w, w_mean=0;
	for(i=0; i<N; i++){
		KS_sample(ks, &pt, &part, &w, -1, NULL);
		w_mean += w;
	}
	return w_mean / N;
}

void KS_destroy(KSource* ks){
	PList_destroy(ks->plist);
	MSV_destroy(ks->metric);
	free(ks);
}


MultiSource* MS_create(int len, KSource** s, double* ws){
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

MultiSource* MS_open(int len, char** filenames, double* ws){
	KSource* s[len];
	int i;
	for(i=0; i<len; i++) s[i] = KS_open(filenames[i]);
	return MS_create(len, s, ws);
}

int MS_sample(MultiSource* ms, char* pt, Part* part, double* w, double w_crit, WeightFun bias){
	double y = rand() / ((double)RAND_MAX+1);
	int i, ret;
	if(ms->cdf[ms->len-1] <= 0) i = (int)(y*ms->len);
	else for(i=0; y*ms->cdf[ms->len-1]>ms->cdf[i]; i++);
	ret = KS_sample(ms->s[i], pt, part, w, w_crit, bias);
	*w *= ms->s[i]->J / (ms->J/ms->len);
	return ret;
}

double MS_w_mean(MultiSource* ms, int N){
	double w_mean=0;
	int i;
	for(i=0; i<ms->len; i++) w_mean += KS_w_mean(ms->s[i], N);
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
