#include<stdio.h>
#include<stdlib.h>

#ifndef METRICS_H
#define METRICS_H

#define MAX_RESAMPLES 1000
#define E_MIN 1e-11
#define NAME_MAX_LEN 96


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
Metric* Metric_copy(Metric* from);
void Metric_destroy(Metric* metric);

typedef struct Geometry{
	int ord; // Cantidad de submetricas
	Metric** ms; // Submetricas
	char* bwfilename; // Nombre de archivo con anchos de banda
	FILE* bwfile; // Archivo con anchos de banda

	double* trasl; // Traslacion de la metrica
	double* rot; // Rotacion de la metrica
} Geometry;

Geometry* Geom_create(int ord, Metric** metrics, char* bwfilename, int variable_bw, double* trasl, double* rot);
Geometry* Geom_copy(Geometry* from);
int Geom_perturb(Geometry* geom, Part* part);
int Geom_next(Geometry* geom);
void Geom_destroy(Geometry* geom);

int E_perturb(Metric* metric, Part* part);
int Let_perturb(Metric* metric, Part* part);

int Vol_perturb(Metric* metric, Part* part);
int SurfXY_perturb(Metric* metric, Part* part);
int Guide_perturb(Metric* metric, Part* part);

int Isotrop_perturb(Metric* metric, Part* part);
int Polar_perturb(Metric* metric, Part* part);


#endif
