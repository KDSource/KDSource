
#include "kdsource/kdsource.h"
#include <stdint.h>

typedef double (*kds_stateless_rng_fct_t)(void);

typedef struct {
  kds_stateless_rng_fct_t rngfct;
} kds_stateless_rng_fct_data_t;

double kdsource_stateless_rng_wrapper(void *thefct_as_state) {
  kds_stateless_rng_fct_data_t *data =
      (kds_stateless_rng_fct_data_t *)(thefct_as_state);
  return data->rngfct();
}

void kdsource_resample_to_mcpl(kds_stateless_rng_fct_t rng_stateless,
                               const char *kds_sourcefile,
                               const char *destination_mcpl, uint64_t nout) {
  const char *filename = kds_sourcefile;
  const char *outfilename = destination_mcpl;

  kds_rng_fct_t rng = kdsource_stateless_rng_wrapper;
  kds_stateless_rng_fct_data_t rng_data;
  rng_data.rngfct = rng_stateless;
  void *rngstate = (void *)(&rng_data);

  KDSource *kds = KDS_open(filename);

  char *realoutfilename = mcpl_name_helper(outfilename, 'B');
  mcpl_outfile_t file = mcpl_create_outfile(realoutfilename);
  free(realoutfilename);
  mcpl_hdr_set_srcname(file, "KDSource resample");

  double w_crit = KDS_w_mean(kds, 1000, NULL);

  printf("Resampling...\n");
  uint64_t i;
  mcpl_particle_t *part = mcpl_get_empty_particle(file);
  for (i = 0; i < nout; ++i) {
    KDS_rand_sample2(rng, rngstate, kds, part, 1, w_crit, NULL, 1);
    mcpl_add_particle(file, part);
  }
  mcpl_closeandgzip_outfile(file);
  KDS_destroy(kds);
  printf("Successfully sampled %llu particles.\n", (unsigned long long)nout);
}
