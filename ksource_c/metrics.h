#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifndef METRICS_H
#define METRICS_H

#define BW_MIN 1e-5


typedef struct Part Part;
typedef struct Metric Metric;

typedef int (*PerturbFun)(Metric* metric, Part* part);

typedef struct Metric{
	int n; // Multiplicidad de la metrica
	int* dims; // Dimensiones de cada submetrica
	double** bw; // Anchos de banda
	FILE* file_bw; // Archivo con anchos de banda
	int variable_bw; // Si es true, se releen anchos de banda en cada sampleo

	PerturbFun* perturb; // Funciones de perturbacion

	double* trasl; // Traslacion de la metrica
	double* rot; // Rotacion de la metrica
	double** geom_par; // Parametros geometricos de cada submetrica
} Metric;

Metric* Metric_create(int n, int* dims, double** bw, char* bwfilename, int variable_bw, PerturbFun* perturb,
	double trasl[3], double rot[3], double** geom_par);
int Metric_perturb(Metric* metric, Part* part);
int Metric_next(Metric* metric);
void Metric_destroy(Metric* metric);

Metric* MetricSimple_create(int dim, double* bw, char* bwfilename, int variable_bw, PerturbFun perturb,
	double trasl[3], double rot[3], double* geom_par);
void MetricSimple_destroy(Metric* metric);

Metric* MetricSepVar_create(int dims[3], double* bw[3], char* bwfilename, int variable_bw, PerturbFun perturb[3],
	double trasl[3], double rot[3], double** geom_par);
void MetricSepVar_destroy(Metric* metric);

int E_perturb(Metric* metric, Part* part);
int Let_perturb(Metric* metric, Part* part);

int Vol_perturb(Metric* metric, Part* part);
int SurfXY_perturb(Metric* metric, Part* part);
int Guide_perturb(Metric* metric, Part* part);

int Isotrop_perturb(Metric* metric, Part* part);
int Polar_perturb(Metric* metric, Part* part);


#endif
