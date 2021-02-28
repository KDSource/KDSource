#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "ksource.h"


double rand_norm(){
	double y1 = (double)(rand()+1) / ((double)RAND_MAX+1), y2 = rand() / (double)RAND_MAX;
	return sqrt(-2 * log(y1)) * cos(2 * M_PI * y2);
}

double *traslv(double *vect, double *trasl, int inverse){
	int i;
	if(!inverse) for(i=0; i<3; i++) vect[i] += trasl[i];
	else for(i=0; i<3; i++) vect[i] -= trasl[i];
	return vect;
}

double *rotv(double *vect, double *rotvec, int inverse){
	double theta, kdotv, kcrossv[3], vrot[3];
	int i;
	theta = sqrt(rotvec[0]*rotvec[0]+rotvec[1]*rotvec[1]+rotvec[2]*rotvec[2]);
	if(theta == 0) return vect;
	if(inverse) theta = -theta;

	kdotv = rotvec[0]*vect[0]+rotvec[1]*vect[1]+rotvec[2]*vect[2];
	kcrossv[0] = rotvec[1]*vect[2]-rotvec[2]*vect[1];
	kcrossv[1] = rotvec[2]*vect[0]-rotvec[0]*vect[2];
	kcrossv[2] = rotvec[0]*vect[1]-rotvec[1]*vect[0];

	for(i=0; i<3; i++) vrot[i] = vect[i]*cos(theta) + (kcrossv[i]/fabs(theta))*sin(theta) + rotvec[i]*kdotv/(theta*theta)*(1-cos(theta));
	for(i=0; i<3; i++) vect[i] = vrot[i];
	return vect;
}
