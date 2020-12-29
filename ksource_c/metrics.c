#include<stdio.h>
#include<stdlib.h>

#include "ksource.h"


Metric* Metric_create(double* geom_par_, char* kernel_, double bw_,
	PerturbPartialFun perturb_E_, PerturbPartialFun perturb_pos_, PerturbPartialFun perturb_dir_,
	PerturbFun perturb_){

	Metric* metric = (Metric*)malloc(sizeof(Metric));
	metric->geom_par = geom_par_;
	metric->kernel = kernel_;
	metric->bw = bw_;,
	metric->perturb_E = perturb_E_;
	metric->perturb_pos = perturb_pos_;
	metric->perturb_dir = perturb_dir_;
	metric->perturb = perturb_;
	return metric;
}

void Metric_next(Metric* metric){
	if(metric->update_bw) metric->update_bw(metric);
}

void Metric_destroy(Metric* metric){
	free(metric);
}