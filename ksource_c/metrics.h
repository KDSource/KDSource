#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifndef METRICS_H
#define METRICS_H

#define BW_MIN 1e-5
#define E_MIN 1e-11


typedef struct Part Part;
typedef struct Metric Metric;

typedef int (*PerturbFun)(Metric* metric, Part* part);

typedef struct Metric{
	int dim; // Dimensiones de cada submetrica
	double* bw; // Anchos de banda
	PerturbFun perturb; // Funcion de perturbacion
	int n_gp; // Cantidad de parametros geometricos
	double* geom_par; // Parametros geometricos de cada submetrica
} Metric;

Metric* Metric_create(int dim, double* bw, PerturbFun perturb, int n_gp, double* geom_par);
void Metric_destroy(Metric* metric);

typedef struct MetricSepVar{
	int ord; // Cantidad de submetricas
	Metric** ms; // Submetricas
	FILE* file_bw; // Archivo con anchos de banda

	double* trasl; // Traslacion de la metrica
	double* rot; // Rotacion de la metrica
} MetricSepVar;

MetricSepVar* MetricSepVar_create(int ord, Metric** metrics, char* bwfilename, int variable_bw, double* trasl, double* rot);
int MetricSepVar_perturb(MetricSepVar* metric, Part* part);
int MetricSepVar_next(MetricSepVar* metric);
void MetricSepVar_destroy(MetricSepVar* metric);

int E_perturb(Metric* metric, Part* part);
int Let_perturb(Metric* metric, Part* part);

int Vol_perturb(Metric* metric, Part* part);
int SurfXY_perturb(Metric* metric, Part* part);
int Guide_perturb(Metric* metric, Part* part);

int Isotrop_perturb(Metric* metric, Part* part);
int Polar_perturb(Metric* metric, Part* part);


#endif
