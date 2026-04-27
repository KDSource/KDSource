
#include "kdsource/kdsource.h"
#include <stdint.h>

void kdsource_resample_to_mcpl( const char * kds_sourcefile,
                                const char * destination_mcpl,
                                uint64_t nout )
{
  const char *filename = kds_sourcefile;
  const char *outfilename = destination_mcpl;

  KDSource* kds = KDS_open(filename);

  char * realoutfilename = mcpl_name_helper( outfilename, 'B' );
  mcpl_outfile_t file = mcpl_create_outfile(realoutfilename);
  free(realoutfilename);
  mcpl_hdr_set_srcname(file, "KDSource resample");

  double w_crit = KDS_w_mean(kds, 1000, NULL);

  printf("Resampling...\n");
  uint64_t i;
  mcpl_particle_t* part = mcpl_get_empty_particle(file);
  for(i=0; i<nout; ++i){
    KDS_sample2(kds, part, 1, w_crit, NULL, 1);
    mcpl_add_particle(file, part);
  }
  mcpl_closeandgzip_outfile(file);
  KDS_destroy(kds);
  printf("Successfully sampled %llu particles.\n", (unsigned long long)nout);
}
