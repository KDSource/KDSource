#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "ksource.h"

KSource* KS_create(double J_, PList* plist_, Metric* metric_){
	KSource* ks = (KSource*)malloc(sizeof(KSource));
	ks->J = J_;
	ks->plist = plist_;
	ks->metric = metric_;
	return ks;
}

int KS_sample(KSource* ks, char* pt, Part* part, double* w, int normalize_w){
	*pt = ks->plist->pt;
	if(!normalize_w){
		PList_get(ks->plist, part, w);
		PList_next(ks->plist);
		Metric_next(ks->metric);
	}
	else{ // Normalizo w a 1
		int resamples = 0;
		while(1){
			PList_get(ks->plist, part, w);
			if(*w > 1){ // Si w>1, uso 1/w como prob de avanzar en la lista
				if(rand() < 1/(*w)*RAND_MAX){
					PList_next(ks->plist);
					Metric_next(ks->metric);
				}
				break;
			}
			else{ // Si w<1, uso w como prob de tomar la particula
				PList_next(ks->plist);
				Metric_next(ks->metric);
				if(rand() < (*w)*RAND_MAX) break;
			}
			if(resamples++ > MAX_RESAMPLES) return 1;
		}
		*w = 1;
	}
	Metric_perturb(ks->metric, part);
	return 0;
}

void KS_destroy(KSource* ks){
	free(ks);
}


MultiSource* MS_create(int len_, KSource** s_, double* ws_){
	MultiSource* ms = (MultiSource*)malloc(sizeof(MultiSource));
	ms->len = len_;
	ms->s = s_;
	ms->ws = ws_;
	ms->cdf = (double*)malloc(ms->len*sizeof(double));
	int i;
	for(i=0; i<ms->len; i++) ms->cdf[i] = ms->ws[i];
	for(i=1; i<ms->len; i++) ms->cdf[i] += ms->cdf[i-1];
	return ms;
}

double MS_J(MultiSource* ms){
	double J = 0;
	int i;
	for(i=0; i<ms->len; i++) J += ms->s[i]->J;
	return J;
}

int MS_sample(MultiSource* ms, char* pt, Part* part, double* w, int normalize_w){
	double y = rand() / ((double)RAND_MAX+1);
	int i;
	if(ms->cdf[ms->len-1] <= 0) i = (int)(y*ms->len);
	else for(i=0; y*ms->cdf[ms->len-1]>ms->cdf[i]; i++);
	return KS_sample(ms->s[i], pt, part, w, normalize_w);
}

void MS_destroy(MultiSource* ms){
	free(ms->cdf);
	free(ms);
}


double rand_norm(){
	double y1 = (double)(rand()+1) / ((double)RAND_MAX+1), y2 = rand() / (double)RAND_MAX;
	return sqrt(-2 * log(y1)) * cos(2 * M_PI * y2);
}