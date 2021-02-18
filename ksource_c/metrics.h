#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifndef METRICS_H
#define METRICS_H

#define BW_MIN 1e-5


typedef struct Part Part;
typedef struct Metric Metric;

typedef void (*MetricVoidFun)(Metric* metric);
typedef int (*PerturbFun)(Metric* metric, Part* part);

typedef struct Metric{
	double** bw;
	MetricVoidFun update_bw;

	double* trasl;
	double* rot;
	double** geom_par;

	int n; // Cantidad de funciones de perturbacion
	PerturbFun* perturb;
} Metric;

Metric* Metric_create(double** bw, int n, PerturbFun* perturb, double trasl[3], double rot[3], double** geom_par);
int Metric_perturb(Metric* metric, Part* part);
void Metric_next(Metric* metric);
void Metric_destroy(Metric* metric);

Metric* MetricSimple_create(double* bw, PerturbFun perturb, double trasl[3], double rot[3], double* geom_par);
void MetricSimple_destroy(Metric* metric);

Metric* MetricSepVar_create(double* bw[3], PerturbFun perturb[3], double trasl[3], double rot[3], double* geom_par[3]);
void MetricSepVar_destroy(Metric* metric);

int E_perturb(Metric* metric, Part* part);
int Let_perturb(Metric* metric, Part* part);

int Vol_perturb(Metric* metric, Part* part);
int SurfXY_perturb(Metric* metric, Part* part);
int Guide_perturb(Metric* metric, Part* part);

int Isotrop_perturb(Metric* metric, Part* part);
int Polar_perturb(Metric* metric, Part* part);


#endif
