
#include "mcpl.h"
#include <stdlib.h>

//flags encoding:
//   1*ENABLE_DOUBLE_PRECISION
// + 2*ENABLE_UNIVERSAL_PDGCODE
// + 4*ENABLE_UNIVERSAL_WEIGHT
//
// TODO: Accept polx/poly/polz/userflags vars (NULL indicates the corresponding
// features should not be enabled).

void kdsource_write_mcpl( const char * filename,
                          uint64_t nparticles,
                          unsigned flags,
                          const double * ekin,
                          const double * x,
                          const double * y,
                          const double * z,
                          const double * ux,
                          const double * uy,
                          const double * uz,
                          const double * time,
                          const double * weight,
                          const int32_t * pdgcode )
{
  //decode flags:
  const unsigned enable_doubleprec = flags % 2;
  flags /= 2;
  const unsigned enable_universal_pdgcode = flags % 2;
  flags /= 2;
  const unsigned enable_universal_weight = flags % 2;

  char * realoutfilename = mcpl_name_helper( filename, 'B' );
  mcpl_outfile_t file = mcpl_create_outfile(realoutfilename);
  free(realoutfilename);
  mcpl_hdr_set_srcname(file, "KDSource MCPL writer");
  mcpl_particle_t* part = mcpl_get_empty_particle(file);
  if ( enable_doubleprec )
    mcpl_enable_doubleprec(file);
  if ( enable_universal_pdgcode )
    mcpl_enable_universal_pdgcode(file,pdgcode[0]);
  if ( enable_universal_weight )
    mcpl_enable_universal_weight(file,weight[0]);

  for ( uint64_t i = 0; i < nparticles; ++i ) {
    //NB: Leaving userflags and polarisation 0.
    part->ekin = ekin[i];
    part->position[0] = x[i];
    part->position[1] = y[i];
    part->position[2] = z[i];
    part->direction[0] = ux[i];
    part->direction[1] = uy[i];
    part->direction[2] = uz[i];
    part->time = time[i];
    if ( !enable_universal_weight )
      part->weight = weight[i];
    if ( !enable_universal_pdgcode )
      part->pdgcode = pdgcode[i];
    mcpl_add_particle(file, part);
  }
  mcpl_closeandgzip_outfile(file);
}
