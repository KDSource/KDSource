#include<stdio.h>
#include<stdlib.h>

#ifndef AUX_H
#define AUX_H


double rand_norm();

double *traslv(double *vect, double *trasl, int inverse);

double *rotv(double *vect, double *rotvec, int inverse);

long pt2pdg(char pt);
char pdg2pt(long pdgcode);


#endif
