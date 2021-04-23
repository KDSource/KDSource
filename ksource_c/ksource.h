#include<stdio.h>
#include<stdlib.h>

#ifndef KSOURCE_H
#define KSOURCE_H

#include "mcpl.h"

#include "aux.h"
#include "metrics.h"
#include "plists.h"


typedef struct KSource KSource;
typedef struct MultiSource MultiSource;
typedef struct PList PList;
typedef struct Metric Metric;

typedef double (*WeightFun)(const mcpl_particle_t* part);

typedef struct KSource{
	double J;
	PList* plist;
	Geometry* geom;
} KSource;

KSource* KS_create(double J, PList* plist, Geometry* geom);
KSource* KS_open(const char* filename);
int KS_sample(KSource* ks, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias);
double KS_w_mean(KSource* ks, int N, WeightFun bias);
void KS_destroy(KSource* ks);

typedef struct MultiSource{
	int len;
	KSource** s;
	double J; // Corriente total
	double* ws; // Pesos de cada fuente
	double* cdf; // cdf de los pesos de fuentes
} MultiSource;

MultiSource* MS_create(int len, KSource** s, const double* ws);
MultiSource* MS_open(int len, const char** filenames, const double* ws);
int MS_sample(MultiSource* ms, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias);
double MS_w_mean(MultiSource* ms, int N, WeightFun bias);
void MS_destroy(MultiSource* ms);


#endif
