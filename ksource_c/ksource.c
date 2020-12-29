#include<stdio.h>
#include<stdlib.h>

#include "ksource.h"


KSource* KS_create(double J, PList* plist, Metric* metric){
	KSource* ks = (KSource*)malloc(sizeof(KSource));
	ks->J = J_;
	ks->plist = plist_;
	ks->metric = metric_;
	return ks;
}

int KS_sample(KSource* ks, char* pt, Part* part, double* w, int normalize_w){
	pt = plist->pt;
	if(!normalize_w) PList_get(ks->plist, part, w);
	else{ // Normalizo w a 1
		int resamples = 0;
		// Si w<1, uso w como prob de tomar la particula
		while(true){
			PList_get(ks->plist, part, w);
			if(*w > 1) break;
			if(distr(gen) < *w) break;
			if(resamples++ > MAX_RESAMPLES) throw SamplingError();
		}
		// Si w>1, uso 1/w como prob de avanzar en la lista
		if(w<1){
			PList_next(ks->plist);
			Metric_next(ks->metric);
		}
		else if(distr(gen) < 1/w){
			PList_next(ks->plist);
			Metric_next(ks->metric);
		}
		w = 1;
	}
	ks->metric->perturb(ks->metric, part);
}

void KS_destroy(KSource* ks){
	free(ks);
}