#include<stdio.h>
#include<stdlib.h>

#ifndef AUX_H
#define AUX_H


double rand_norm();

double *traslv(double *vect, const double *trasl, int inverse);

double *rotv(double *vect, const double *rotvec, int inverse);

long pt2pdg(char pt);
char pdg2pt(long pdgcode);

double interp(double x, const double *xs, const double *ys, int N);

double H10_n_ARN(double E);
double H10_p_ARN(double E);
double H10_n_ICRP(double E);
double H10_p_ICRP(double E);


#endif
