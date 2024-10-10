#ifndef KDSOURCE_H
#define KDSOURCE_H

#include <stdio.h>
#include <stdlib.h>

/***********************************************************************************/
/*                                                                                 */
/*  KDSource: KDE particle sources */
/*                                                                                 */
/*  Utilities for sampling particles from a KDE source. KDSource sources use */
/*  particle lists in MCPL format, and apply on them the Kernel Density
 * Estimation */
/*  (KDE) method. This allows coupling different Monte Carlo radiation transport
 */
/*  simulations, and gives variance reduction. */
/*                                                                                 */
/*  Find more information and updates at https://github.com/KDSource/KDSource */
/*                                                                                 */
/*  This file can be freely used as per the terms in the LICENSE file. */
/*                                                                                 */
/*  Written by Osiris Inti Abbate, 2021. */
/*                                                                                 */
/***********************************************************************************/

#include "mcpl.h"

#include "KDSourceConfig.h"
#include "geom.h"
#include "plist.h"
#include "utils.h"

#define MAX_RESAMPLES 1000000
#define NAME_MAX_LEN 256

typedef double (*WeightFun)(const mcpl_particle_t *part);

typedef struct KDSource {
  double J;       // Total current [1/s]
  char kernel;    // Kernel
  PList *plist;   // Particle list
  Geometry *geom; // Geometry defining variable treatment
} KDSource;

KDSource *KDS_create(double J, char kernel, PList *plist, Geometry *geom);
KDSource *KDS_open(const char *xmlfilename);
int KDS_sample2(KDSource *kds, mcpl_particle_t *part, int perturb,
                double w_crit, WeightFun bias, int loop);
int KDS_sample(KDSource *kds, mcpl_particle_t *part);
double KDS_w_mean(KDSource *kds, int N, WeightFun bias);
void KDS_destroy(KDSource *kds);

typedef struct MultiSource {
  int len;      // Number of sources
  KDSource **s; // Array of sources
  double J;     // Total current [1/s]
  double *ws;   // Frequency weights of sources
  double *cdf;  // cdf of sources weights
} MultiSource;

MultiSource *MS_create(int len, KDSource **s, const double *ws);
MultiSource *MS_open(int len, const char **xmlfilenames, const double *ws);
int MS_sample2(MultiSource *ms, mcpl_particle_t *part, int perturb,
               double w_crit, WeightFun bias, int loop);
int MS_sample(MultiSource *ms, mcpl_particle_t *part);
double MS_w_mean(MultiSource *ms, int N, WeightFun bias);
void MS_destroy(MultiSource *ms);

#endif
