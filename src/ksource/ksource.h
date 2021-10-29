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


#define MAX_RESAMPLES 1000000
#define NAME_MAX_LEN 256


typedef double (*WeightFun)(const mcpl_particle_t* part);

typedef struct KSource{
	double J;       // Total current [1/s]
	PList* plist;   // Particle list
	Geometry* geom; // Geometry defining variable treatment
} KSource;

KSource* KS_create(double J, PList* plist, Geometry* geom);
KSource* KS_open(const char* xmlfilename);
int KS_sample2(KSource* ks, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias, int loop);
int KS_sample(KSource* ks, mcpl_particle_t* part);
double KS_w_mean(KSource* ks, int N, WeightFun bias);
void KS_destroy(KSource* ks);

typedef struct MultiSource{
	int len;     // Number of sources
	KSource** s; // Array of sources
	double J;    // Total current [1/s]
	double* ws;  // Frequency weights of sources
	double* cdf; // cdf of sources weights
} MultiSource;

MultiSource* MS_create(int len, KSource** s, const double* ws);
MultiSource* MS_open(int len, const char** xmlfilenames, const double* ws);
int MS_sample2(MultiSource* ms, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias, int loop);
int MS_sample(MultiSource* ms, mcpl_particle_t* part);
double MS_w_mean(MultiSource* ms, int N, WeightFun bias);
void MS_destroy(MultiSource* ms);


#endif
