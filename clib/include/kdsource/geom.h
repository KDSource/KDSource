#ifndef KDSOURCE_GEOM_H
#define KDSOURCE_GEOM_H

#include <stdio.h>
#include <stdlib.h>

/***********************************************************************************/
/*                                                                                 */
/*  Utilities for geometry and metrics handling for KDSource. */
/*                                                                                 */
/*  This file can be freely used as per the terms in the LICENSE file. */
/*                                                                                 */
/*  Written by Osiris Inti Abbate, 2021. */
/*                                                                                 */
/***********************************************************************************/

#include "mcpl.h"

typedef struct Metric Metric;

#ifndef KDS_RNG_FCT
#define KDS_RNG_FCT
typedef double (*kds_rng_fct_t)(void);
#endif

typedef int (*PerturbFun)(kds_rng_fct_t, const Metric *metric,
                          mcpl_particle_t *part, double bw, char kernel);

struct Metric {
  int dim;            // Dimension
  float *scaling;     // Variables scaling
  PerturbFun perturb; // Perturbation function
  int nps;            // Number of metric parameters
  double *params;     // Metric parameters
  kds_rng_fct_t rng;  // Random number generator
};

Metric *Metric_create(int dim, const double *scaling, PerturbFun perturb,
                      int nps, const double *params);
Metric *Metric_copy(const Metric *from);
void Metric_destroy(Metric *metric);

typedef struct Geometry {
  int ord;          // Number of submetrics
  Metric **ms;      // Submetrics
  char *bwfilename; // Bandwidth file name
  FILE *bwfile;     // Bandwidth file
  double bw;        // Normalized bandwidth
  char kernel;      // Kernel

  double *trasl; // Geometry translation
  double *rot;   // Geometry rotation
} Geometry;

Geometry *Geom_create(int ord, Metric **metrics, double bw,
                      const char *bwfilename, char kernel, const double *trasl,
                      const double *rot);
Geometry *Geom_copy(const Geometry *from);
int Geom_perturb(kds_rng_fct_t, const Geometry *geom, mcpl_particle_t *part);
int Geom_next(Geometry *geom, int loop);
void Geom_destroy(Geometry *geom);

int E_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
              double bw, char kernel);
int Let_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                double bw, char kernel);
int wl_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
               double bw, char kernel);

int t_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
              double bw, char kernel);
int Dec_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                double bw, char kernel);

int Vol_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                double bw, char kernel);
int SurfXY_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                   double bw, char kernel);
int SurfR_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                  double bw, char kernel);
int SurfR2_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                   double bw, char kernel);
int SurfCircle_perturb(kds_rng_fct_t, const Metric *metric,
                       mcpl_particle_t *part, double bw, char kernel);
int Guide_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                  double bw, char kernel);

int Isotrop_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                    double bw, char kernel);
int Polar_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                  double bw, char kernel);
int PolarMu_perturb(kds_rng_fct_t, const Metric *metric, mcpl_particle_t *part,
                    double bw, char kernel);

extern const int _n_metrics;
extern const char *_metric_names[];
extern const PerturbFun _metric_perturbs[];

#endif
