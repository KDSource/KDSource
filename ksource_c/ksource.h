#include<stdio.h>
#include<stdlib.h>

#ifndef KSOURCE_H
#define KSOURCE_H

#include "aux.h"
#include "metrics.h"
#include "plists.h"

#define MAX_RESAMPLES 1000
#define LINE_MAX_LEN 256
#define NAME_MAX_LEN 96


typedef struct KSource KSource;
typedef struct MultiSource MultiSource;
typedef struct PList PList;
typedef struct Metric Metric;

typedef struct Part{
	double E;
	double pos[3];
	double dir[3];
} Part;

typedef int (*KSSampleFun)(KSource* ks, char* pt, Part* part, double* w, int normalize_w);
typedef int (*MSSampleFun)(MultiSource* ms, char* pt, Part* part, double* w, int normalize_w);

typedef struct KSource{
	double J;
	PList* plist;
	MetricSepVar* metric;
} KSource;

KSource* KS_create(double J, PList* plist, MetricSepVar* metric);
int KS_sample(KSource* ks, char* pt, Part* part, double* w, int normalize_w);
void KS_destroy(KSource* ks);

typedef struct MultiSource{
	int len;
	KSource** s;
	double* ws; // Pesos de cada fuente
	double* cdf; // cdf de los pesos de fuentes
} MultiSource;

MultiSource* MS_create(int len, KSource** s, double* ws);
double MS_J(MultiSource* ms);
int MS_sample(MultiSource* ms, char* pt, Part* part, double* w, int normalize_w);
void MS_destroy(MultiSource* ms);


#endif
