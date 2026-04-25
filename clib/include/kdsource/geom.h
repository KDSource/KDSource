#ifndef METRICS_H
#define METRICS_H

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

typedef int (*PerturbFun)(const Metric *metric, mcpl_particle_t *part,
                          double bw, char kernel);

struct Metric {
  int dim;            // Dimension
  float *scaling;     // Variables scaling
  PerturbFun perturb; // Perturbation function
  int nps;            // Number of metric parameters
  double *params;     // Metric parameters
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
int Geom_perturb(const Geometry *geom, mcpl_particle_t *part);
int Geom_next(Geometry *geom, int loop);
void Geom_destroy(Geometry *geom);

int E_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
              char kernel);
int Let_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                char kernel);
int wl_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
               char kernel);

int t_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
              char kernel);
int Dec_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                char kernel);

int Vol_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                char kernel);
int SurfXY_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                   char kernel);
int SurfR_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                  char kernel);
int SurfR2_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                   char kernel);
int SurfCircle_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                       char kernel);
int Guide_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                  char kernel);

int Isotrop_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                    char kernel);
int Polar_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                  char kernel);
int PolarMu_perturb(const Metric *metric, mcpl_particle_t *part, double bw,
                    char kernel);

static const int _n_metrics = 14;
static const char *_metric_names[] = {
    "Energy", "Lethargy", "Wavelength", "Vol",   "SurfXY",
    "SurfR",  "SurfR2",   "SurfCircle", "Guide", "Isotrop",
    "Polar",  "PolarMu",  "Time",       "Decade"};
static const PerturbFun _metric_perturbs[] = {
    E_perturb,      Let_perturb,     wl_perturb,     Vol_perturb,
    SurfXY_perturb, SurfR_perturb,   SurfR2_perturb, SurfCircle_perturb,
    Guide_perturb,  Isotrop_perturb, Polar_perturb,  PolarMu_perturb,
    t_perturb,      Dec_perturb};

#endif
