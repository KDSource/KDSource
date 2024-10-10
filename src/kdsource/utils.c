#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "kdsource.h"

// Sample with normal distribution (x0=0, s=1)
double rand_norm() {
  double y1 = (double)(rand() + 1) / ((double)RAND_MAX + 1),
         y2 = rand() / (double)RAND_MAX;
  return sqrt(-2 * log(y1)) * cos(2 * M_PI * y2);
}

// Sample with Epanechnikov distribution (x0=0, s=1)
double rand_epan() {
  double x1 = (double)rand() / ((double)RAND_MAX / 2.0) - 1.0;
  double x2 = (double)rand() / ((double)RAND_MAX / 2.0) - 1.0;
  double x3 = (double)rand() / ((double)RAND_MAX / 2.0) - 1.0;
  if (abs(x3) >= abs(x2) && abs(x3) >= abs(x1))
    return x2;
  else
    return x3;
}

// Sample with Tophat distribution (x0=0, s=1)
double rand_box() {
  double y1 = (double)rand() / ((double)RAND_MAX / 2.0) - 1.0;
  return y1;
}

// Sample depending on kernel selected
double rand_type(char kernel) {
  if (kernel == 'g')
    return rand_norm();
  if (kernel == 'e')
    return rand_epan();
  if (kernel == 'b')
    return rand_box();
  else {
    printf("Cannot perturbate with current kernel. \n");
    // KDS_error("Cannot perturbate with current kernel.");
    return 0;
  }
  return 0;
}

// Translation
double *traslv(double *vect, const double *trasl, int inverse) {
  int i;
  if (!inverse)
    for (i = 0; i < 3; i++)
      vect[i] += trasl[i];
  else
    for (i = 0; i < 3; i++)
      vect[i] -= trasl[i];
  return vect;
}

// Rotation (rotvec in axis-angle format)
double *rotv(double *vect, const double *rotvec, int inverse) {
  double theta, kdotv, kcrossv[3], vrot[3];
  int i;
  theta = sqrt(rotvec[0] * rotvec[0] + rotvec[1] * rotvec[1] +
               rotvec[2] * rotvec[2]);
  if (theta == 0)
    return vect;
  if (inverse)
    theta = -theta;

  kdotv = rotvec[0] * vect[0] + rotvec[1] * vect[1] + rotvec[2] * vect[2];
  kcrossv[0] = rotvec[1] * vect[2] - rotvec[2] * vect[1];
  kcrossv[1] = rotvec[2] * vect[0] - rotvec[0] * vect[2];
  kcrossv[2] = rotvec[0] * vect[1] - rotvec[1] * vect[0];

  for (i = 0; i < 3; i++)
    vrot[i] = vect[i] * cos(theta) + (kcrossv[i] / fabs(theta)) * sin(theta) +
              rotvec[i] * kdotv / (theta * theta) * (1 - cos(theta));
  for (i = 0; i < 3; i++)
    vect[i] = vrot[i];
  return vect;
}

// Conversion between char pt and pdg code for particle type
int pt2pdg(char pt) {
  if (pt == 'n')
    return 2112;
  if (pt == 'p')
    return 22;
  if (pt == 'e')
    return 11;
  return 0;
}
char pdg2pt(int pdgcode) {
  if (pdgcode == 2112)
    return 'n';
  if (pdgcode == 22)
    return 'p';
  if (pdgcode == 11)
    return 'e';
  return '0';
}

// Interpolation
double interp(double x, const double *xs, const double *ys, int N) {
  if (x < xs[0] || x > xs[N - 1]) {
    printf("Error in interp: x outside range.\n");
    return 0;
  }
  int i = 0;
  while (x > xs[i + 1])
    i++;
  return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i]);
}

// Dosimetric factors [pSv cm2]

#define N_n_ARN 15
const double log_E_n_ARN[N_n_ARN] = {-1e3,      -1.75e+01, -6.21e+00, -3.69e+00,
                                     -1.94e+00, -1.39e+00, -5.71e-01, 1.82e-01,
                                     9.16e-01,  1.03e+00,  1.16e+00,  1.61e+00,
                                     2.69e+00,  2.94e+00,  3.91e+00};
const double log_f_n_ARN[N_n_ARN] = {-1e3,     2.36e+00, 2.04e+00, 2.96e+00,
                                     4.84e+00, 5.31e+00, 5.84e+00, 6.05e+00,
                                     6.03e+00, 6.02e+00, 6.02e+00, 6.00e+00,
                                     6.28e+00, 6.37e+00, 6.37e+00};

#define N_p_ARN 27
const double log_E_p_ARN[N_p_ARN] = {
    -1e3,      -4.61e+00, -4.20e+00, -3.91e+00, -3.51e+00, -3.22e+00, -3.00e+00,
    -2.81e+00, -2.53e+00, -2.30e+00, -1.90e+00, -1.61e+00, -1.20e+00, -9.16e-01,
    -6.93e-01, -5.11e-01, -2.23e-01, 0.00e+00,  4.05e-01,  6.93e-01,  1.10e+00,
    1.39e+00,  1.61e+00,  1.79e+00,  2.08e+00,  2.30e+00,  3.91e+00};
const double log_f_p_ARN[N_p_ARN] = {
    -1e3,      -2.80e+00, -1.86e-01, 4.88e-02,  -2.11e-01, -4.46e-01, -5.98e-01,
    -6.73e-01, -6.35e-01, -4.94e-01, -1.17e-01, 1.82e-01,  5.88e-01,  8.67e-01,
    1.08e+00,  1.24e+00,  1.48e+00,  1.65e+00,  1.93e+00,  2.15e+00,  2.41e+00,
    2.60e+00,  2.74e+00,  2.87e+00,  3.07e+00,  3.24e+00,  3.24e+00};

#define N_n_ICRP 69
const double log_E_n_ICRP[N_n_ICRP] = {
    -1e3,      -2.07e+01, -1.84e+01, -1.75e+01, -1.61e+01, -1.54e+01, -1.45e+01,
    -1.38e+01, -1.31e+01, -1.22e+01, -1.15e+01, -1.08e+01, -9.90e+00, -9.21e+00,
    -8.52e+00, -7.60e+00, -6.91e+00, -6.21e+00, -5.30e+00, -4.61e+00, -3.91e+00,
    -3.51e+00, -3.00e+00, -2.66e+00, -2.30e+00, -1.90e+00, -1.61e+00, -1.20e+00,
    -6.93e-01, -3.57e-01, -1.05e-01, 0.00e+00,  1.82e-01,  4.05e-01,  6.93e-01,
    1.10e+00,  1.39e+00,  1.61e+00,  1.79e+00,  1.95e+00,  2.08e+00,  2.20e+00,
    2.30e+00,  2.48e+00,  2.64e+00,  2.71e+00,  2.77e+00,  2.89e+00,  3.00e+00,
    3.04e+00,  3.40e+00,  3.91e+00,  4.32e+00,  4.61e+00,  4.87e+00,  5.01e+00,
    5.19e+00,  5.30e+00,  5.70e+00,  5.99e+00,  6.21e+00,  6.40e+00,  6.55e+00,
    6.68e+00,  6.80e+00,  6.91e+00,  7.60e+00,  8.52e+00,  9.21e+00};
const double log_f_n_ICRP[N_n_ICRP] = {
    -1e3,     1.13e+00, 1.27e+00, 1.39e+00, 1.65e+00, 1.77e+00, 1.89e+00,
    1.95e+00, 2.00e+00, 2.04e+00, 2.06e+00, 2.06e+00, 2.06e+00, 2.05e+00,
    2.05e+00, 2.02e+00, 2.02e+00, 2.03e+00, 2.08e+00, 2.21e+00, 2.50e+00,
    2.75e+00, 3.14e+00, 3.42e+00, 3.74e+00, 4.10e+00, 4.37e+00, 4.74e+00,
    5.18e+00, 5.45e+00, 5.63e+00, 5.71e+00, 5.80e+00, 5.90e+00, 6.01e+00,
    6.13e+00, 6.18e+00, 6.20e+00, 6.21e+00, 6.21e+00, 6.21e+00, 6.21e+00,
    6.21e+00, 6.21e+00, 6.20e+00, 6.20e+00, 6.19e+00, 6.18e+00, 6.17e+00,
    6.16e+00, 6.12e+00, 6.07e+00, 6.08e+00, 6.10e+00, 6.10e+00, 6.10e+00,
    6.10e+00, 6.10e+00, 6.16e+00, 6.24e+00, 6.28e+00, 6.34e+00, 6.44e+00,
    6.46e+00, 6.47e+00, 6.50e+00, 6.65e+00, 6.95e+00, 7.24e+00};

#define N_p_ICRP 65
const double log_E_p_ICRP[N_p_ICRP] = {
    -1e3,      -5.30e+00, -5.12e+00, -4.96e+00, -4.83e+00, -4.71e+00, -4.61e+00,
    -4.42e+00, -4.34e+00, -4.20e+00, -4.07e+00, -3.91e+00, -3.69e+00, -3.51e+00,
    -3.22e+00, -3.00e+00, -2.81e+00, -2.66e+00, -2.53e+00, -2.30e+00, -1.90e+00,
    -1.61e+00, -1.20e+00, -9.16e-01, -6.93e-01, -6.71e-01, -5.11e-01, -4.12e-01,
    -2.23e-01, 0.00e+00,  1.13e-01,  2.85e-01,  4.05e-01,  6.93e-01,  1.10e+00,
    1.39e+00,  1.61e+00,  1.79e+00,  1.81e+00,  2.08e+00,  2.30e+00,  2.71e+00,
    3.00e+00,  3.40e+00,  3.69e+00,  3.91e+00,  4.09e+00,  4.38e+00,  4.61e+00,
    5.01e+00,  5.30e+00,  5.70e+00,  5.99e+00,  6.21e+00,  6.40e+00,  6.68e+00,
    6.91e+00,  7.31e+00,  7.60e+00,  8.01e+00,  8.29e+00,  8.52e+00,  8.70e+00,
    8.99e+00,  9.21e+00};
const double log_f_p_ICRP[N_p_ICRP] = {
    -1e3,      -4.31e+00, -4.10e+00, -3.79e+00, -3.40e+00, -3.02e+00, -2.68e+00,
    -2.25e+00, -2.10e+00, -1.86e+00, -1.71e+00, -1.49e+00, -1.29e+00, -1.16e+00,
    -1.05e+00, -9.97e-01, -9.44e-01, -8.89e-01, -8.14e-01, -6.58e-01, -2.92e-01,
    0.00e+00,  4.12e-01,  6.93e-01,  9.04e-01,  9.24e-01,  1.07e+00,  1.15e+00,
    1.32e+00,  1.50e+00,  1.59e+00,  1.72e+00,  1.81e+00,  2.01e+00,  2.28e+00,
    2.46e+00,  2.60e+00,  2.71e+00,  2.72e+00,  2.92e+00,  3.10e+00,  3.41e+00,
    3.64e+00,  3.94e+00,  4.12e+00,  4.28e+00,  4.41e+00,  4.59e+00,  4.70e+00,
    4.87e+00,  4.97e+00,  5.08e+00,  5.15e+00,  5.20e+00,  5.23e+00,  5.28e+00,
    5.33e+00,  5.36e+00,  5.46e+00,  5.53e+00,  5.59e+00,  5.62e+00,  5.65e+00,
    5.70e+00,  5.73e+00};

double H10_n_ARN(double E) {
  double log_E = log(E);
  double log_f = interp(log_E, log_E_n_ARN, log_f_n_ARN, N_n_ARN);
  return exp(log_f);
}
double H10_p_ARN(double E) {
  double log_E = log(E);
  double log_f = interp(log_E, log_E_p_ARN, log_f_p_ARN, N_p_ARN);
  return exp(log_f);
}
double H10_n_ICRP(double E) {
  double log_E = log(E);
  double log_f = interp(log_E, log_E_n_ICRP, log_f_n_ICRP, N_n_ICRP);
  return exp(log_f);
}
double H10_p_ICRP(double E) {
  double log_E = log(E);
  double log_f = interp(log_E, log_E_p_ICRP, log_f_p_ICRP, N_p_ICRP);
  return exp(log_f);
}
