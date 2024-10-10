#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "kdsource.h"

void display_usage() {
  printf("Usage: kdtool resample sourcefile [options]\n\n");
  printf("Resample particles from source defined in XML file sourcefile, and "
         "save them in\n");
  printf("a MCPL file.\n\n");
  printf("Options:\n");
  printf("\t-o outfile: Name of MCPL file with new samples\n");
  printf("\t            (default: \"resampled.mcpl\").\n");
  printf("\t-n N:       Number of new samples (default: 1E5).\n");
  printf("\t-h, --help: Display usage instructions.\n");
}

int resample_parse_args(int argc, char **argv, const char **filename,
                        const char **outfilename, long int *N) {
  *filename = 0;
  *outfilename = 0;
  *N = 1E5;
  int i;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '\0')
      continue;
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      display_usage();
      exit(0);
    }
    if (strcmp(argv[i], "-o") == 0) {
      *outfilename = argv[++i];
      continue;
    }
    if (strcmp(argv[i], "-n") == 0) {
      *N = atof(argv[++i]);
      continue;
    }
    if (argv[i][0] == '-') {
      printf("Error: Invalid argument: %s\n", argv[i]);
      exit(1);
    }
    if (!*filename) {
      *filename = argv[i];
      continue;
    }
    printf("Too many arguments. Use -h or --help for help.\n");
    exit(1);
  }
  if (!*filename) {
    printf("No source file. Use -h or --help for help.\n");
    exit(1);
  }
  if (!*outfilename)
    *outfilename = "resampled.mcpl";
  return 0;
}

int main(int argc, char *argv[]) {
  const char *filename;
  const char *outfilename;
  long int N;

  if (resample_parse_args(argc, argv, &filename, &outfilename, &N))
    return 1;

  KDSource *kds = KDS_open(filename);
  mcpl_particle_t part;

  mcpl_outfile_t file = mcpl_create_outfile(outfilename);
  mcpl_hdr_set_srcname(file, "KDSource resample");

  double w_crit = KDS_w_mean(kds, 1000, NULL);

  printf("Resampling...\n");
  long int i;
  for (i = 0; i < N; i++) {
    KDS_sample2(kds, &part, 1, w_crit, NULL, 1);
    mcpl_add_particle(file, &part);
  }
  mcpl_closeandgzip_outfile(file);
  printf("Successfully sampled %ld particles.\n", N);

  return 0;
}
