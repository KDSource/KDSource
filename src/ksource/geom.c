#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<math.h>

#include "ksource.h"


void KS_error(const char* msg);
void KS_end(const char* msg);

Metric* Metric_create(int dim, const double* scaling, PerturbFun perturb, int nps, const double* params){
	Metric* metric = (Metric*)malloc(sizeof(Metric));
	int i;
	metric->dim = dim;
	metric->scaling = (float*)malloc(dim*sizeof(float));
	if(scaling) for(i=0; i<dim; i++) metric->scaling[i] = (float)scaling[i];
	else for(i=0; i<dim; i++) metric->scaling[i] = 0;
	metric->perturb = perturb;
	metric->nps = nps;
	metric->params = (double*)malloc(nps * sizeof(double));
	for(i=0; i<nps; i++) metric->params[i] = params[i];
	return metric;
}

Metric* Metric_copy(const Metric* from){
	Metric* metric = (Metric*)malloc(sizeof(Metric));
	*metric = *from;
	int i;
	metric->scaling = (float*)malloc(metric->dim*sizeof(double));
	for(i=0; i<metric->dim; i++) metric->scaling[i] = from->scaling[i];
	metric->params = (double*)malloc(metric->nps * sizeof(double));
	for(i=0; i<metric->nps; i++) metric->params[i] = from->params[i];
	return metric;
}

void Metric_destroy(Metric* metric){
	free(metric->scaling);
	free(metric);
}

Geometry* Geom_create(int ord, Metric** metrics, double bw, const char* bwfilename,
		const double* trasl, const double* rot){
	Geometry* geom = (Geometry*)malloc(sizeof(Geometry));
	geom->ord = ord;
	int i;
	geom->ms = (Metric**)malloc(ord * sizeof(Metric*));
	for(i=0; i<ord; i++) geom->ms[i] = metrics[i];
	geom->bw = bw;
	geom->bwfilename = NULL;
	geom->bwfile = NULL;
	if(bwfilename) if(strlen(bwfilename)){
		FILE* bwfile;
		if((bwfile=fopen(bwfilename, "r")) == 0){
			printf("Could not open file %s\n", bwfilename);
			KS_error("Error in Geom_create");
		}
		geom->bwfilename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
		strcpy(geom->bwfilename, bwfilename);
		geom->bwfile = bwfile;
		Geom_next(geom, 1);
	}
	if(trasl){
		geom->trasl = (double*)malloc(3 * sizeof(double));
		for(i=0; i<3; i++) geom->trasl[i] = trasl[i];
	}
	else geom->trasl = NULL;
	if(rot){
		geom->rot = (double*)malloc(3 * sizeof(double));
		for(i=0; i<3; i++) geom->rot[i] = rot[i];
	}
	else geom->rot = NULL;
	return geom;
}

Geometry* Geom_copy(const Geometry* from){
	Geometry* geom = (Geometry*)malloc(sizeof(Geometry));
	*geom = *from;
	int i;
	geom->ms = (Metric**)malloc(geom->ord * sizeof(Metric*));
	for(i=0; i<geom->ord; i++) geom->ms[i] = Metric_copy(from->ms[i]);
	geom->bwfilename = NULL;
	geom->bwfile = NULL;
	if(geom->bwfile){
		geom->bwfilename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
		strcpy(geom->bwfilename, from->bwfilename);
		geom->bwfile = fopen(geom->bwfilename, "r");
		Geom_next(geom, 1);
	}
	if(from->trasl){
		geom->trasl = (double*)malloc(3 * sizeof(double));
		for(i=0; i<3; i++) geom->trasl[i] = from->trasl[i];
	}
	else geom->trasl = NULL;
	if(from->rot){
		geom->rot = (double*)malloc(3 * sizeof(double));
		for(i=0; i<3; i++) geom->rot[i] = from->rot[i];
	}
	else geom->rot = NULL;
	return geom;
}

int Geom_perturb(const Geometry* geom, mcpl_particle_t* part){
	int i, ret=0;
	if(geom->trasl) traslv(part->position, geom->trasl, 1);
	if(geom->rot){ rotv(part->position, geom->rot, 1); rotv(part->direction, geom->rot, 1); }
	for(i=0; i<geom->ord; i++)
		ret += geom->ms[i]->perturb(geom->ms[i], part, geom->bw);
	if(geom->rot){ rotv(part->position, geom->rot, 0); rotv(part->direction, geom->rot, 0); }
	if(geom->trasl) traslv(part->position, geom->trasl, 0);
	return ret;
}

int Geom_next(Geometry* geom, int loop){
	if(geom->bwfile){
		float temp;
		if(fread(&temp, sizeof(float), 1, geom->bwfile) == 1){
			geom->bw = temp;
			return 0;
		}
		if(!loop) KS_end("Geom_next: End of BW file reached.");
		rewind(geom->bwfile); // After 1st failed try, rewind
		if(fread(&temp, sizeof(float), 1, geom->bwfile) == 1){
			geom->bw = temp;
			return 1;
		}
		KS_error("Error in Geom_next: Could not read BW.");
	}
	return 0;
}

void Geom_destroy(Geometry* geom){
	int i;
	for(i=0; i<geom->ord; i++) Metric_destroy(geom->ms[i]);
	free(geom->ms);
	free(geom->trasl);
	free(geom->rot);
	free(geom);
}


int E_perturb(const Metric* metric, mcpl_particle_t* part, double bw){
	int cont=0;
	double E0 = part->ekin, E=E0;
	E = E0 + bw*metric->scaling[0] * rand_norm();
	while(E < 0  && cont++ < MAX_RESAMPLES){ // Keep E > 0
		E = E0 + bw*metric->scaling[0] * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in E_perturb: MAX_RESAMPLES reached (E0 = %le).\n", E0);
	part->ekin = E;
	return 0;
}
int Let_perturb(const Metric* metric, mcpl_particle_t* part, double bw){
	part->ekin *= exp(bw*metric->scaling[0] * rand_norm());
	return 0;
}

int Vol_perturb(const Metric* metric, mcpl_particle_t* part, double bw){
	int cont=0;
	double x0=part->position[0], y0=part->position[1], z0=part->position[1], x=x0, y=y0, z=z0;
	double xmin=metric->params[0], xmax=metric->params[1], ymin=metric->params[2], ymax=metric->params[3], zmin=metric->params[4], zmax=metric->params[5];
	x = x0 + bw*metric->scaling[0] * rand_norm();
	while((x < xmin || x > xmax)  && cont++ < MAX_RESAMPLES){ // Keet x within range [xmin,xmax]
		x = x0 + bw*metric->scaling[0] * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in Vol_perturb: MAX_RESAMPLES reached (x0 = %lf).\n", x0);
	cont = 0;
	y = y0 + bw*metric->scaling[1] * rand_norm();
	while((y < ymin || y > ymax)  && cont++ < MAX_RESAMPLES){ // Keet y within range [ymin,ymax]
		y = y0 + bw*metric->scaling[1] * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in Vol_perturb: MAX_RESAMPLES reached (y0 = %lf).\n", y0);
	cont = 0;
	z = z0 + bw*metric->scaling[2] * rand_norm();
	while((z < zmin || z > zmax)  && cont++ < MAX_RESAMPLES){ // Keet z within range [zmin,zmax]
		z = z0 + bw*metric->scaling[2] * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in Vol_perturb: MAX_RESAMPLES reached (z0 = %lf).\n", z0);
	part->position[0] = x; part->position[1] = y; part->position[2] = z;
	return 0;
}
int SurfXY_perturb(const Metric* metric, mcpl_particle_t* part, double bw){
	int cont=0;
	double x0=part->position[0], y0=part->position[1], x=x0, y=y0;
	double xmin=metric->params[0], xmax=metric->params[1], ymin=metric->params[2], ymax=metric->params[3];
	double z = metric->params[4];
	x = x0 + bw*metric->scaling[0] * rand_norm();
	while((x < xmin || x > xmax)  && cont++ < MAX_RESAMPLES){ // Keet x within range [xmin,xmax]
		x = x0 + bw*metric->scaling[0] * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in SurfXY_perturb: MAX_RESAMPLES reached (x0 = %lf).\n", x0);
	cont = 0;
	y = y0 + bw*metric->scaling[1] * rand_norm();
	while((y < ymin || y > ymax)  && cont++ < MAX_RESAMPLES){ // Keet y within range [ymin,ymax]
		y = y0 + bw*metric->scaling[1] * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in SurfXY_perturb: MAX_RESAMPLES reached (y0 = %lf).\n", y0);
	part->position[0] = x; part->position[1] = y; part->position[2] = z;
	return 0;
}
int Guide_perturb(const Metric* metric, mcpl_particle_t* part, double bw){
	double x=part->position[0], y=part->position[1], z0=part->position[2], z=z0, dx=part->direction[0], dy=part->direction[1], dz=part->direction[2];
	double xwidth=metric->params[0], yheight=metric->params[1], zmax=metric->params[2], rcurv=metric->params[3];
	double t, theta, phi, theta0, dx2, dz2;
	int cont=0, mirror;
	if(rcurv != 0){ // Transform to curved guide variables
		double r = sqrt((rcurv+x)*(rcurv+x) + z*z);
		x = copysign(1, rcurv) * r - rcurv; z = fabs(rcurv) * asin(z / r);
		dx2 = dx; dz2 = dz; dx = dx2*cos(z/rcurv) + dz2*sin(z/rcurv); dz = -dx2*sin(z/rcurv) + dz2*cos(z/rcurv);
	}
	// Transform from (x,y,z,dx,dy,dz) to (z,t,theta,phi)
	if((y/yheight > -x/xwidth) && (y/yheight <  x/xwidth))      mirror=0; // mirror x pos
	else if((y/yheight >  x/xwidth) && (y/yheight > -x/xwidth)) mirror=1; // mirror y pos
	else if((y/yheight < -x/xwidth) && (y/yheight >  x/xwidth)) mirror=2; // mirror x neg
	else                                                        mirror=3; // mirror y neg
	switch(mirror){
		case 0:
			t = 0.5*yheight + y;
			theta0 = acos(dx); phi = atan2(-dy, dz);
			break;
		case 1:
			t = 1.0*yheight + 0.5*xwidth - x;
			theta0 = acos(dy); phi = atan2(dx, dz);
			break;
		case 2:
			t = 1.5*yheight + 1.0*xwidth - y;
			theta0 = acos(-dx); phi = atan2(dy, dz);
			break;
		case 3:
			t = 2.0*yheight + 1.5*xwidth + x;
			theta0 = acos(-dy); phi = atan2(-dx, dz);
			break;
	}
	// Perturb
	z = z0 + bw*metric->scaling[0] * rand_norm();
	while((z < 0 || z > zmax)  && cont++ < MAX_RESAMPLES){ // Keep z within range [0,zmax]
		z = z0 + bw*metric->scaling[0] * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in Guide_perturb: MAX_RESAMPLES reached (z0 = %lf).\n", z0);
	t += bw*metric->scaling[1] * rand_norm();
	while(t < 0) t += 2*(xwidth+yheight);
	while(t > 2*(xwidth+yheight)) t -= 2*(xwidth+yheight);
	theta = theta0 + bw*metric->scaling[2]*M_PI/180 * rand_norm();
	cont = 0;
	while(cos(theta0)*cos(theta) < 0 && cont++ < MAX_RESAMPLES){ // Avoid perturbation to change propagation direction
		theta = theta0 + bw*metric->scaling[2]*M_PI/180 * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in Guide_perturb: MAX_RESAMPLES reached (theta0 = %lf).\n", theta0);
	phi += bw*metric->scaling[3]*M_PI/180 * rand_norm();
	switch(mirror){
		case 0: if(t<  0)              t=  0;              else if(t>  yheight)          t=  yheight;          break;
		case 1: if(t<  yheight)        t=  yheight;        else if(t>  yheight+  xwidth) t=  yheight+  xwidth; break;
		case 2: if(t<  yheight+xwidth) t=  yheight+xwidth; else if(t>2*yheight+  xwidth) t=2*yheight+  xwidth; break;
		case 3: if(t<2*yheight+xwidth) t=2*yheight+xwidth; else if(t>2*yheight+2*xwidth) t=2*yheight+2*xwidth; break;
	}
	// Antitransform from (z,t,theta_n,theta_t) to (x,y,z,dx,dy,dz)
	switch(mirror){
		case 0:
			x =  xwidth/2; y =  t - 0.5*yheight;
			dx = cos(theta); dz = sin(theta)*cos(phi); dy = -sin(theta)*sin(phi);
			break;
		case 1:
			y =  yheight/2; x = -t + yheight + 0.5*xwidth;
			dy = cos(theta); dz = sin(theta)*cos(phi); dx = sin(theta)*sin(phi);
			break;
		case 2:
			x = -xwidth/2; y = -t + 1.5*yheight + xwidth;
			dx = -cos(theta); dz = sin(theta)*cos(phi); dy = sin(theta)*sin(phi);
			break;
		case 3:
			y = -yheight/2; x =  t - 2*yheight - 1.5*xwidth;
			dy = -cos(theta); dz = sin(theta)*cos(phi); dx = -sin(theta)*sin(phi);
			break;
	}
	if(rcurv != 0){ // Antitransform curved guide variables
		double r = (rcurv + x) * copysign(1, rcurv), ang = z / rcurv;
		x = copysign(1, rcurv) * r * cos(ang) - rcurv; z = r * sin(fabs(ang));
		dx2 = dx; dz2 = dz; dx = dx2*cos(ang) - dz2*sin(ang); dz = dx2*sin(ang) + dz2*cos(ang);
	}
	part->position[0] = x; part->position[1] = y; part->position[2] = z;
	part->direction[0] = dx; part->direction[1] = dy; part->direction[2] = dz;
	return 0;
}

int Isotrop_perturb(const Metric* metric, mcpl_particle_t* part, double bw){
	if(bw*metric->scaling[0] == INFINITY){
		part->direction[2] = -1 + 2.*rand()/RAND_MAX;
		double dxy = sqrt(1-part->direction[2]*part->direction[2]);
		double phi = 2.*M_PI*rand()/RAND_MAX;
		part->direction[0] = dxy*cos(phi);
		part->direction[1] = dxy*sin(phi);
	}
	else if(bw*metric->scaling[0] > 0){
		double xi = (double)rand()/RAND_MAX;
		double w = 1;
		w += bw*metric->scaling[0]*bw*metric->scaling[0] * log(xi+(1-xi)*exp(-2/(bw*metric->scaling[0]*bw*metric->scaling[0])));
		double phi = 2.*M_PI*rand()/RAND_MAX;
		double uv = sqrt(1-w*w), u = uv*cos(phi), v = uv*sin(phi);
		double x=part->direction[0], y=part->direction[1], z=part->direction[2];
		if(part->direction[2] > 0){
			part->direction[0] = u*z + w*x - (v*x-u*y)*y/(1+z);
			part->direction[1] = v*z + w*y + (v*x-u*y)*x/(1+z);
		}
		else{
			part->direction[0] = u*z + w*x + (v*x-u*y)*y/(1-z);
			part->direction[1] = v*z + w*y - (v*x-u*y)*x/(1-z);
		}
		part->direction[2] = w*z - u*x - v*y;
	}
	return 0;
}

int Polar_perturb(const Metric* metric, mcpl_particle_t* part, double bw){
	double theta, phi, theta0;
	int cont=0;
	theta0 = acos(part->direction[2]);
	phi   = atan2(part->direction[1], part->direction[0]);
	theta = theta0 + bw*metric->scaling[0]*M_PI/180 * rand_norm();
	while(cos(theta0)*cos(theta) < 0 && cont++ < MAX_RESAMPLES){ // Avoid perturbation to change propagation direction
		theta = theta0 + bw*metric->scaling[0]*M_PI/180 * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in Polar_perturb: MAX_RESAMPLES reached (theta0 = %lf).\n", theta0);
	phi   += bw*metric->scaling[1]*M_PI/180 * rand_norm();
	part->direction[0] = sin(theta) * cos(phi);
	part->direction[1] = sin(theta) * sin(phi);
	part->direction[2] = cos(theta);
	return 0;
}

int PolarMu_perturb(const Metric* metric, mcpl_particle_t* part, double bw){
	double mu, phi, mu0;
	int cont=0;
	mu0 = part->direction[2];
	phi   = atan2(part->direction[1], part->direction[0]);
	mu = mu0 + bw*metric->scaling[0] * rand_norm();
	while((mu*mu0 < 0 || fabs(mu) > 1) && cont++ < MAX_RESAMPLES){ // Avoid perturbation to change propagation direction
		mu = mu0 + bw*metric->scaling[0] * rand_norm();
	}
	if(cont > MAX_RESAMPLES) printf("Warning in PolarMu_perturb: MAX_RESAMPLES reached (mu0 = %lf).\n", mu0);
	phi   += bw*metric->scaling[1]*M_PI/180 * rand_norm();
	part->direction[0] = sqrt(1-mu*mu) * cos(phi);
	part->direction[1] = sqrt(1-mu*mu) * sin(phi);
	part->direction[2] = mu;
	return 0;
}
