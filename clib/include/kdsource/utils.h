#ifndef KDSOURCE_UTILS_H
#define KDSOURCE_UTILS_H

#include <stdio.h>
#include <stdlib.h>

/***********************************************************************************/
/*                                                                                 */
/*  General utilities for KDSource. */
/*                                                                                 */
/*  This file can be freely used as per the terms in the LICENSE file. */
/*                                                                                 */
/*  Written by Osiris Inti Abbate, 2021. */
/*                                                                                 */
/***********************************************************************************/

// random sampling:
#ifndef KDS_RNG_FCT
#define KDS_RNG_FCT
typedef double (*kds_rng_fct_t)(void *);
#endif
double kds_rand_norm(kds_rng_fct_t, void *rngstate);
double kds_rand_epan(kds_rng_fct_t, void *rngstate);
double kds_rand_box(kds_rng_fct_t, void *rngstate);
double kds_rand_type(kds_rng_fct_t, void *rngstate, char kernel);

double *traslv(double *vect, const double *trasl, int inverse);
double *rotv(double *vect, const double *rotvec, int inverse);

int pt2pdg(char pt);
char pdg2pt(int pdgcode);

double interp(double x, const double *xs, const double *ys, int N);

double H10_n_ARN(double E);
double H10_p_ARN(double E);
double H10_n_ICRP(double E);
double H10_p_ICRP(double E);

#endif
