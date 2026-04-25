#ifndef AUX_H
#define AUX_H

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

double rand_norm();
double rand_epan();
double rand_box();
double rand_type(char kernel);

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
