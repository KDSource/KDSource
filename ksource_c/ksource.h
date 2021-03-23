#include<stdio.h>
#include<stdlib.h>

#ifndef KSOURCE_H
#define KSOURCE_H

#include "aux.h"
#include "metrics.h"
#include "plists.h"


typedef struct KSource KSource;
typedef struct MultiSource MultiSource;
typedef struct PList PList;
typedef struct Metric Metric;

typedef struct Part Part;

void Part_print(Part* part);

typedef int (*KSSampleFun)(KSource* ks, char* pt, Part* part, double* w, double w_crit);
typedef int (*MSSampleFun)(MultiSource* ms, char* pt, Part* part, double* w, double w_crit);
typedef double (*WeightFun)(Part* part);

typedef struct KSource{
	double J;
	PList* plist;
	Geometry* geom;
} KSource;

KSource* KS_create(double J, PList* plist, Geometry* geom);
KSource* KS_open(char* filename);
int KS_sample(KSource* ks, char* pt, Part* part, double* w, double w_crit, WeightFun bias);
double KS_w_mean(KSource* ks, int N);
void KS_destroy(KSource* ks);

typedef struct MultiSource{
	int len;
	KSource** s;
	double J; // Corriente total
	double* ws; // Pesos de cada fuente
	double* cdf; // cdf de los pesos de fuentes
} MultiSource;

MultiSource* MS_create(int len, KSource** s, double* ws);
MultiSource* MS_open(int len, char** filenames, double* ws);
int MS_sample(MultiSource* ms, char* pt, Part* part, double* w, double w_crit, WeightFun bias);
double MS_w_mean(MultiSource* ms, int N);
void MS_destroy(MultiSource* ms);


#endif
