#ifndef PLISTS_H
#define PLISTS_H

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

/***********************************************************************************/
/*                                                                                 */
/*  Utilities for particle lists handling for KDSource. */
/*                                                                                 */
/*  This file can be freely used as per the terms in the LICENSE file. */
/*                                                                                 */
/*  Written by Osiris Inti Abbate, 2021. */
/*                                                                                 */
/***********************************************************************************/

#include "mcpl.h"

typedef struct PList {
  char pt;     // Particle type ("n", "p", "e", ...)
  int pdgcode; // PDG code for particle type

  char *filename;   // Name of MCPL file
  mcpl_file_t file; // MCPL file

  double *trasl; // PList translation
  double *rot;   // PList rotation
  int x2z;       // If true, apply permutation x,y,z -> y,z,x

  const mcpl_particle_t *part; // Pointer to selected particle
  long long npts;
} PList;

PList *PList_create(char pt, const char *filename, const double *trasl,
                    const double *rot, int switch_x2z);
PList *PList_copy(const PList *from);
int PList_get(const PList *plist, mcpl_particle_t *part);
int PList_next(PList *plist, int loop);
void PList_destroy(PList *plist);

#endif
