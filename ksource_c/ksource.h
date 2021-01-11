#include<stdio.h>
#include<stdlib.h>

#ifndef KSOURCE_H
#define KSOURCE_H

#include "aux.h"
#include "metrics.h"
#include "plists.h"

#define MAX_RESAMPLES 1000


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
	Metric* metric;
} KSource;

KSource* KS_create(double J_, PList* plist_, Metric* metric_);
int KS_sample(KSource* ks, char* pt, Part* part, double* w, int normalize_w);
void KS_destroy(KSource* ks);

typedef struct MultiSource{
	int len;
	KSource** s;
	double* ws;
	double* cdf;
} MultiSource;

MultiSource* MS_create(int len_, KSource** s_, double* ws_);
double MS_J(MultiSource* ms);
int MS_sample(MultiSource* ms, char* pt, Part* part, double* w, int normalize_w);
void MS_destroy(MultiSource* ms);


#endif
