#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifndef METRICS_H
#define METRICS_H


void (*MetricVoidFun)(Metric* metric);
typedef int (*PerturbFun)(Metric* m, Part* part);

struct Metric{
	double* geom_par;

	char* kernel;
	double* bw;
	MetricVoidFun update_bw;

	PerturbFun perturb_E;
	PerturbFun perturb_pos;
	PerturbFun perturb_dir;
	PerturbFun perturb;
};

Metric* Metric_create(double* geom_par_, char* kernel_, double bw_,
	PerturbFun perturb_E_, PerturbFun perturb_pos_, PerturbFun perturb_dir_, PerturbFun perturb_);
void Metric_next(Metric* metric);
void Metric_destroy(Metric* metric);


#endif