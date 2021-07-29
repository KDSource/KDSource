#ifndef KSOURCE_H
#define KSOURCE_H

#include<stdio.h>
#include<stdlib.h>

/***********************************************************************************/
/*                                                                                 */
/*  KSource: KDE particle sources                                                  */
/*                                                                                 */
/*  Utilities for sampling particles from a KDE source. KSource sources use        */
/*  particle lists in MCPL format, and apply on them the Kernel Density Estimation */
/*  (KDE) method. This allows coupling different Monte Carlo radiation transport   */
/*  simulations, and gives variance reduction.                                     */
/*                                                                                 */
/*  Find more information and updates at https://github.com/inti-abbate/KSource    */
/*                                                                                 */
/*  This file can be freely used as per the terms in the LICENSE file.             */
/*                                                                                 */
/*  Written by Osiris Inti Abbate, 2021.                                           */
/*                                                                                 */
/***********************************************************************************/

#include "mcpl.h"

#include "KSourceConfig.h"
#include "utils.h"
#include "geom.h"
#include "plist.h"


typedef struct KSource KSource;
typedef struct MultiSource MultiSource;

typedef double (*WeightFun)(const mcpl_particle_t* part);

struct KSource{
	double J;
	PList* plist;
	Geometry* geom;
};

KSource* KS_create(double J, PList* plist, Geometry* geom);
KSource* KS_open(const char* filename);
int KS_sample2(KSource* ks, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias, int loop);
int KS_sample(KSource* ks, mcpl_particle_t* part);
double KS_w_mean(KSource* ks, int N, WeightFun bias);
void KS_destroy(KSource* ks);

struct MultiSource{
	int len;
	KSource** s;
	double J;    // Total current [1/s]
	double* ws;  // Weights of each source
	double* cdf; // cdf of sources weights
};

MultiSource* MS_create(int len, KSource** s, const double* ws);
MultiSource* MS_open(int len, const char** filenames, const double* ws);
int MS_sample2(MultiSource* ms, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias, int loop);
int MS_sample(MultiSource* ms, mcpl_particle_t* part);
double MS_w_mean(MultiSource* ms, int N, WeightFun bias);
void MS_destroy(MultiSource* ms);


#endif
