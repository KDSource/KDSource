#include<stdio.h>
#include<stdlib.h>

#ifndef METRICS_H
#define METRICS_H

#include "mcpl.h"

#define MAX_RESAMPLES 1000
#define E_MIN 1e-11
#define E_MAX 20
#define NAME_MAX_LEN 256


typedef struct Metric Metric;

typedef int (*PerturbFun)(const Metric* metric, mcpl_particle_t* part, double bw);

struct Metric{
	int dim; // Dimensiones de cada submetrica
	float* scaling; // Escaleos de variables
	PerturbFun perturb; // Funcion de perturbacion
	int nps; // Cantidad de parametros geometricos
	double* params; // Parametros geometricos de cada submetrica
};

Metric* Metric_create(int dim, const double* scaling, PerturbFun perturb, int nps, const double* params);
Metric* Metric_copy(const Metric* from);
void Metric_destroy(Metric* metric);

typedef struct Geometry{
	int ord; // Cantidad de submetricas
	Metric** ms; // Submetricas
	char* bwfilename; // Nombre de archivo con anchos de banda
	FILE* bwfile; // Archivo con anchos de banda
	double bw; // Ancho de banda normalizado

	double* trasl; // Traslacion de la metrica
	double* rot; // Rotacion de la metrica
} Geometry;

Geometry* Geom_create(int ord, Metric** metrics, double bw, const char* bwfilename,
	const double* trasl, const double* rot);
Geometry* Geom_copy(const Geometry* from);
int Geom_perturb(const Geometry* geom, mcpl_particle_t* part);
int Geom_next(Geometry* geom);
void Geom_destroy(Geometry* geom);

int E_perturb(const Metric* metric, mcpl_particle_t* part, double bw);
int Let_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

int Vol_perturb(const Metric* metric, mcpl_particle_t* part, double bw);
int SurfXY_perturb(const Metric* metric, mcpl_particle_t* part, double bw);
int Guide_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

int Isotrop_perturb(const Metric* metric, mcpl_particle_t* part, double bw);
int Polar_perturb(const Metric* metric, mcpl_particle_t* part, double bw);


#endif
