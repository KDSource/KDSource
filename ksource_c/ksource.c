#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "ksource.h"


void Part_print(Part* part){
	printf("(%.3le, (%.3lf, %.3lf, %.3lf), (%.3lf, %.3lf, %.3lf))\n",
		part->E, part->pos[0], part->pos[1], part->pos[2], part->dir[0], part->dir[1], part->dir[2]);
}

KSource* KS_create(double J, PList* plist, Geometry* geom){
	KSource* ks = (KSource*)malloc(sizeof(KSource));
	ks->J = J;
	ks->plist = plist;
	ks->geom = geom;
	return ks;
}

KSource* KS_open(char* filename){
	int i, j;
	char buffer[LINE_MAX_LEN], *pbuffer;
	double J;

	char pt;
	int ord_plist;
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
	fscanf(file, "%le\n", &J); // Leer J
	// PList
	fgets(buffer, LINE_MAX_LEN, file); // # PList
	fscanf(file, "%c\n", &pt); // Leer pt
	fscanf(file, "%d\n", &ord_plist); // Leer ord_plist
	char formats[ord_plist][NAME_MAX_LEN], *tracksfiles[ord_plist];
	ReadFun reads[ord_plist];
	for(i=0; i<ord_plist; i++) fgets(formats[i], NAME_MAX_LEN, file); // Leer formatos
	for(i=0; i<ord_plist; i++){
		tracksfiles[i] = (char*)malloc(NAME_MAX_LEN*sizeof(char));
		fgets(tracksfiles[i], NAME_MAX_LEN, file); // Leer nombres de archivos de tracks
		tracksfiles[i][strcspn(tracksfiles[i], "\n")] = 0;
	}
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
	if(variable_bw){
		bwfilename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
		fgets(bwfilename, NAME_MAX_LEN, file);
		bwfilename[strcspn(bwfilename, "\n")] = 0;
		for(i=0; i<ord_geom; i++) bws[i] = NULL;
	}
	else{
		for(i=0; i<ord_geom; i++){
			bws[i] = (double*)malloc(dims[i]*sizeof(double));
			for(j=0; j<dims[i]; j++) fscanf(file, "%lf", &bws[i][j]);
		}
	}
	// Crear PList
	for(i=0; i<ord_plist; i++){
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
	PList* plist = PList_create(pt, ord_plist, tracksfiles, reads, trasl_plist, rot_plist, switch_x2z);
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
	for(i=0; i<ord_plist; i++) free(tracksfiles[i]);
	fclose(file);

	return s;
}

int KS_sample(KSource* ks, char* pt, Part* part, double* w, double w_crit, WeightFun bias){
	*pt = ks->plist->pt;
	if(w_crit <= 0){
		PList_get(ks->plist, part, w);
		PList_next(ks->plist);
		Geom_next(ks->geom);
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
					Geom_next(ks->geom);
				}
				break;
			}
			else{ // Si w*bs<w_crit, uso w*bs/w_crit como prob de tomar la particula
				PList_next(ks->plist);
				Geom_next(ks->geom);
				if(rand() < (*w*bs)/w_crit*RAND_MAX) break;
			}
			if(resamples++ > MAX_RESAMPLES) return 1;
		}
		*w = 1/bs;
	}
	Geom_perturb(ks->geom, part);
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
	Geom_destroy(ks->geom);
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
