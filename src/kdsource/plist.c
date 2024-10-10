#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "kdsource.h"

void KDS_error(const char *msg);
void KDS_end(const char *msg);

PList *PList_create(char pt, const char *filename, const double *trasl,
                    const double *rot, int switch_x2z) {
  PList *plist = (PList *)malloc(sizeof(PList));
  plist->pt = pt;
  plist->pdgcode = pt2pdg(pt);
  mcpl_file_t file = mcpl_open_file(filename);
  plist->npts = mcpl_hdr_nparticles(file);
  plist->filename = (char *)malloc(NAME_MAX_LEN * sizeof(char));
  strcpy(plist->filename, filename);
  plist->file = file;
  int i;
  if (trasl) {
    plist->trasl = (double *)malloc(3 * sizeof(double));
    for (i = 0; i < 3; i++)
      plist->trasl[i] = trasl[i];
  } else
    plist->trasl = NULL;
  if (rot) {
    plist->rot = (double *)malloc(3 * sizeof(double));
    for (i = 0; i < 3; i++)
      plist->rot[i] = rot[i];
  } else
    plist->rot = NULL;
  plist->x2z = switch_x2z;
  plist->part = NULL;
  PList_next(plist, 1);
  return plist;
}

PList *PList_copy(const PList *from) {
  PList *plist = (PList *)malloc(sizeof(PList));
  *plist = *from;
  plist->filename = (char *)malloc(NAME_MAX_LEN * sizeof(char));
  strcpy(plist->filename, from->filename);
  plist->file = mcpl_open_file(plist->filename);
  int i;
  if (from->trasl) {
    plist->trasl = (double *)malloc(3 * sizeof(double));
    for (i = 0; i < 3; i++)
      plist->trasl[i] = from->trasl[i];
  }
  if (from->rot) {
    plist->rot = (double *)malloc(3 * sizeof(double));
    for (i = 0; i < 3; i++)
      plist->rot[i] = from->rot[i];
  }
  plist->part = NULL;
  PList_next(plist, 1);
  return plist;
}

int PList_get(const PList *plist, mcpl_particle_t *part) {
  *part = *plist->part;
  if (plist->trasl)
    traslv(part->position, plist->trasl, 0);
  if (plist->rot) {
    rotv(part->position, plist->rot, 0);
    rotv(part->direction, plist->rot, 0);
  }
  if (plist->x2z) {
    double x = part->position[0], y = part->position[1], z = part->position[2];
    part->position[0] = y;
    part->position[1] = z;
    part->position[2] = x;
    double dx = part->direction[0], dy = part->direction[1],
           dz = part->direction[2];
    part->direction[0] = dy;
    part->direction[1] = dz;
    part->direction[2] = dx;
  }
  return 0;
}

int PList_next(PList *plist, int loop) {
  int ret = 0, resamples = 0;
  while (resamples++ < MAX_RESAMPLES) {
    plist->part = mcpl_read(plist->file);
    if (plist->part == NULL) { // After 1st failed try, rewind
      if (!loop)
        KDS_end("PList_next: End of particle list reached.");
      ret = 1;
      mcpl_rewind(plist->file);
      plist->part = mcpl_read(plist->file);
      if (plist->part == NULL)
        KDS_error("Error in PList_next: Could not get particle");
    }
    if (plist->part->weight > 0 && plist->part->pdgcode == plist->pdgcode)
      return ret;
  }
  KDS_error("Error in PList_next: MAX_RESAMPLES reached.");
  return 1;
}

void PList_destroy(PList *plist) {
  free(plist->filename);
  mcpl_close_file(plist->file);
  free(plist->trasl);
  free(plist->rot);
  free(plist);
}
