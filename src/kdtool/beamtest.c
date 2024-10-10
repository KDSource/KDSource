#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "kdsource.h"

void display_usage() {
  printf("Usage: kdtool beamtest sourcefile [options]\n\n");
  printf("Executes a simple simulation with source defined in XML file "
         "sourcefile, in\n");
  printf("which calculates the number of particles passing thru a rectangular "
         "collimator.\n");
  printf("The simulation is repeated using the particle list directly as "
         "source, to\n");
  printf("compare the results.\n\n");
  printf("Results are computed for 4 energy groups, and stored in a results "
         "file which\n");
  printf("can be imported from a spreadsheet.\n\n");
  printf("This tool is designed to be used with flat neutron sources with "
         "particles\n");
  printf("propagating towards z direction.\n\n");
  printf("Options:\n");
  printf("\t-n N:           Number of source particles (default: 1E6).\n");
  printf("\t-o results:     Name of file to store results.\n");
  printf("\t-xwidth value:  Width of the collimator, in cm (default: 7).\n");
  printf("\t-yheight value: Height of the collimator, in cm (default: 20).\n");
  printf("\t-z value:       Position of the collimator along z axis, in cm\n");
  printf("\t                (default: 500).\n");
  printf(
      "\t-xshift:        Horizontal shift of the center of the collimator,\n");
  printf("\t                in cm (default: 0)\n");
  printf("\t-yshift:        Vertical shift of the center of the collimator,\n");
  printf("\t                in cm (default: 0)\n");
  printf("\t-h, --help:     Display usage instructions.\n");
}

int beamtest_parse_args(int argc, char **argv, const char **sourcefile, long *N,
                        const char **results, double *xwidth, double *yheight,
                        double *z, double *xshift, double *yshift) {
  *sourcefile = 0;
  *N = 1E6;
  *results = "results_beamtest.txt";
  *xwidth = 7;
  *yheight = 20;
  *z = 500;
  *xshift = 0;
  *yshift = 0;
  int i;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '\0')
      continue;
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      display_usage();
      exit(0);
    }
    if (strcmp(argv[i], "-n") == 0) {
      *N = atof(argv[++i]);
      continue;
    }
    if (strcmp(argv[i], "-o") == 0) {
      *results = argv[++i];
      continue;
    }
    if (strcmp(argv[i], "-xwidth") == 0) {
      *xwidth = atof(argv[++i]);
      continue;
    }
    if (strcmp(argv[i], "-yheight") == 0) {
      *yheight = atof(argv[++i]);
      continue;
    }
    if (strcmp(argv[i], "-z") == 0) {
      *z = atof(argv[++i]);
      continue;
    }
    if (strcmp(argv[i], "-xshift") == 0) {
      *xshift = atof(argv[++i]);
      continue;
    }
    if (strcmp(argv[i], "-yshift") == 0) {
      *yshift = atof(argv[++i]);
      continue;
    }
    if (argv[i][0] == '-') {
      printf("Error: Invalid argument: %s\n", argv[i]);
      exit(1);
    }
    if (!*sourcefile) {
      *sourcefile = argv[i];
      continue;
    }
    printf("Too many arguments. Use -h or --help for help.\n");
    exit(1);
  }
  if (!*sourcefile) {
    printf("No source file. Use -h or --help for help.\n");
    exit(1);
  }
  return 0;
}

// Determine if particle passes thru collimator
int propagate_collimator(mcpl_particle_t *part, double xwidth, double yheight,
                         double z, double xshift, double yshift) {
  if ((part->position[2] - z) * part->direction[2] > 0)
    return 0;
  if (part->ekin <= 0)
    return 0;
  double t = (z - part->position[2]) /
             part->direction[2]; // Time to reach collimator assuming v=1
  double x = part->position[0] + t * part->direction[0];
  double y = part->position[1] + t * part->direction[1];
  return (fabs(x - xshift) < xwidth / 2 && fabs(y - yshift) < yheight / 2);
}

int main(int argc, char *argv[]) {
  // User defined options
  const char *sourcefile, *results;
  double xwidth, yheight, z, xshift, yshift;
  long int N;
  if (beamtest_parse_args(argc, argv, &sourcefile, &N, &results, &xwidth,
                          &yheight, &z, &xshift, &yshift))
    return 1;

  // KDSource object
  KDSource *kds = KDS_open(sourcefile);
  mcpl_particle_t part;
  double w_crit = KDS_w_mean(kds, 1000, NULL);

  // Variables to store results
  long int N_source_KDE = 0, N_source_tracks = 0;
  double I_source_KDE = 0, I_source_tracks = 0, I_KDE_tot = 0, I_tracks_tot = 0,
         err_KDE_tot = 0, err_tracks_tot = 0;
  int ngroups = 4;
  double Ecuts[] = {0, 1E-8, 3E-7, 0.5, 100}; // Energy groups cuts
  double I_KDE[] = {0, 0, 0, 0}, I_tracks[] = {0, 0, 0, 0},
         err_KDE[] = {0, 0, 0, 0}, err_tracks[] = {0, 0, 0, 0};
  long int i;
  int j;

  // KDE simulation
  mcpl_rewind(kds->plist->file);
  if (kds->geom->bwfile)
    rewind(kds->geom->bwfile);
  printf("Executing simulation with KDE source...\n");
  for (i = 0; i < N; i++) {
    KDS_sample2(kds, &part, 1, w_crit, NULL, 1);
    N_source_KDE++;
    I_source_KDE += part.weight;
    if (propagate_collimator(&part, xwidth, yheight, z, xshift, yshift)) {
      for (j = 0; j < ngroups; j++) {
        if (part.ekin > Ecuts[j] && part.ekin < Ecuts[j + 1]) {
          I_KDE[j] += part.weight;
          err_KDE[j] += part.weight * part.weight;
        }
      }
    }
  }
  printf("Finished. Particles simulated: N_source = %ld, I_source = %lf\n\n",
         N_source_KDE, I_source_KDE);

  // Tracks simulation
  mcpl_rewind(kds->plist->file);
  if (kds->geom->bwfile)
    rewind(kds->geom->bwfile);
  printf("Executing simulation with tracks source...\n");
  for (i = 0; i < N; i++) {
    if (KDS_sample2(kds, &part, 0, w_crit, NULL, 1))
      break;
    N_source_tracks++;
    I_source_tracks += part.weight;
    if (propagate_collimator(&part, xwidth, yheight, z, xshift, yshift)) {
      for (j = 0; j < ngroups; j++) {
        if (part.ekin > Ecuts[j] && part.ekin < Ecuts[j + 1]) {
          I_tracks[j] += part.weight;
          err_tracks[j] += part.weight * part.weight;
        }
      }
    }
  }
  printf("Finished. Particles simulated: N_source = %ld, I_source = %lf\n\n",
         N_source_tracks, I_source_tracks);

  // Normalize results
  for (j = 0; j < ngroups; j++) {
    I_KDE[j] /= I_source_KDE;
    I_tracks[j] /= I_source_tracks;
    I_KDE_tot += I_KDE[j];
    I_tracks_tot += I_tracks[j];
    err_KDE_tot += err_KDE[j];
    err_tracks_tot += err_tracks[j];
    err_KDE[j] = sqrt(err_KDE[j]) / I_source_KDE;
    err_tracks[j] = sqrt(err_tracks[j]) / I_source_tracks;
  }
  err_KDE_tot = sqrt(err_KDE_tot) / I_source_KDE;
  err_tracks_tot = sqrt(err_tracks_tot) / I_source_tracks;

  // Save results
  FILE *resfile = fopen(results, "w");
  fprintf(resfile, "Results are given for source intensity = 1.\n");
  fprintf(resfile, "Energy group cuts: [ ");
  for (j = 0; j < ngroups; j++)
    fprintf(resfile, "%.3le, ", Ecuts[j]);
  fprintf(resfile, "%.3le ]\n\n", Ecuts[ngroups]);
  fprintf(
      resfile,
      "E group\tJ KDE\terr\tJ tracks\terr\tDiff.\terr\tRel. diff. [%%]\terr\n");
  double I_dif, err_dif, reldif, err_rd;
  for (j = 0; j < ngroups; j++) {
    I_dif = I_KDE[j] - I_tracks[j];
    err_dif = sqrt(err_KDE[j] * err_KDE[j] + err_tracks[j] * err_tracks[j]);
    reldif = I_dif / I_tracks[j];
    err_rd = reldif * sqrt((err_dif / I_dif) * (err_dif / I_dif) +
                           (err_tracks[j] / I_tracks[j]) *
                               (err_tracks[j] / I_tracks[j]));
    fprintf(resfile,
            "%d\t%.3le\t%.1le\t%.3le\t%.1le\t%.3le\t%.1le\t%.2lf\t%.1lf\n",
            j + 1, I_KDE[j], err_KDE[j], I_tracks[j], err_tracks[j], I_dif,
            err_dif, 100 * reldif, 100 * err_rd);
  }
  I_dif = I_KDE_tot - I_tracks_tot;
  err_dif = sqrt(err_KDE_tot * err_KDE_tot + err_tracks_tot * err_tracks_tot);
  reldif = I_dif / I_tracks_tot;
  err_rd = fabs(reldif * sqrt((err_dif / I_dif) * (err_dif / I_dif) +
                              (err_tracks_tot / I_tracks_tot) *
                                  (err_tracks_tot / I_tracks_tot)));
  fprintf(resfile,
          "Tot.\t%.3le\t%.1le\t%.3le\t%.1le\t%.3le\t%.1le\t%.2lf\t%.1lf\n",
          I_KDE_tot, err_KDE_tot, I_tracks_tot, err_tracks_tot, I_dif, err_dif,
          100 * reldif, 100 * err_rd);
  printf("Results stored in text file %s\n", results);

  return 0;
}