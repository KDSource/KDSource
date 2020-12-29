#include<stdio.h>
#include<stdlib.h>

#include "plists.h"
#include "metrics.h"

#ifndef KSOURCE_H
#define KSOURCE_H


struct KSource;
struct MultiSource;
struct PList;
struct Metric;

struct Part{
	double E;
	double pos[3];
	double dir[3];
};

typedef int (*KSSampleFun)(KSource* ks, char* pt, Part* part, double* w, int normalize_w);
typedef int (*MSSampleFun)(MultiSource* ms, char* pt, Part* part, double* w, int normalize_w);

struct KSource{
	double J;
	PList* plist;
	Metric* metric;
};

KSource* KS_create(double J_, PList* plist_, Metric* metric_);
int KS_sample(KSource* ks, char* pt, Part* part, double* w);
void KS_destroy(KSource* ks);

struct MultiSource{
	KSource* s;
	double* ws;
};

MultiSource* MS_create(KSource* s_, double* ws_);
double MS_J(MultiSource* ms);
int MS_sample(MultiSource* ms, char* pt, Part* part, double* w);
void MS_destroy(MultiSource* ms);


#endif