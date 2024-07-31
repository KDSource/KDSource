#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include<math.h>

#include "kdsource.h"


void KDS_error(const char* msg);
void KDS_end(const char* msg);

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

Geometry* Geom_create(int ord, Metric** metrics, double bw, const char* bwfilename, char kernel,
		const double* trasl, const double* rot){
	Geometry* geom = (Geometry*)malloc(sizeof(Geometry));
	geom->ord = ord;
	int i;
	geom->ms = (Metric**)malloc(ord * sizeof(Metric*));
	for(i=0; i<ord; i++) geom->ms[i] = metrics[i];
	geom->bw = bw;
	geom->bwfilename = NULL;
	geom->bwfile = NULL;
	geom->kernel = kernel;
	if(bwfilename) if(strlen(bwfilename)){
		FILE* bwfile;
		if((bwfile=fopen(bwfilename, "rb")) == 0){
			printf("Could not open file %s\n", bwfilename);
			KDS_error("Error in Geom_create");
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
		geom->bwfile = fopen(geom->bwfilename, "rb");
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
		ret += geom->ms[i]->perturb(geom->ms[i], part, geom->bw, geom->kernel);
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
		if(!loop) KDS_end("Geom_next: End of BW file reached.");
		rewind(geom->bwfile); // After 1st failed try, rewind
		if(fread(&temp, sizeof(float), 1, geom->bwfile) == 1){
			geom->bw = temp;
			return 1;
		}
		KDS_error("Error in Geom_next: Could not read BW.");
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


int E_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	part->ekin += bw*metric->scaling[0] * rand_type(kernel);
	if(part->ekin < 0) part->ekin *= -1;
	return 0;
}
int Let_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	float E = part->ekin;
	E *= exp(bw*metric->scaling[0] * rand_type(kernel));
	while(E > metric->params[0]){
		E = part->ekin;
		E *= exp(bw*metric->scaling[0] * rand_type(kernel));
	}
	part->ekin = E;
	return 0;
}

int t_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	part->time += bw*metric->scaling[0] * rand_type(kernel);
	return 0;
}
int Dec_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	part->time *= pow(10, bw*metric->scaling[0] * rand_type(kernel));
	return 0;
}

int Vol_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	double xmin=metric->params[0], xmax=metric->params[1];
	double ymin=metric->params[2], ymax=metric->params[3];
	double zmin=metric->params[4], zmax=metric->params[5];
	part->position[0] += bw*metric->scaling[0] * rand_type(kernel);
	while((part->position[0] < xmin) || (part->position[0] > xmax)){
		if(part->position[0] < xmin) part->position[0] += 2 * (xmin - part->position[0]);
		else                         part->position[0] -= 2 * (part->position[0] - xmax);
	}
	part->position[1] += bw*metric->scaling[1] * rand_type(kernel);
	while((part->position[1] < ymin) || (part->position[1] > ymax)){
		if(part->position[1] < ymin) part->position[1] += 2 * (ymin - part->position[1]);
		else                         part->position[1] -= 2 * (part->position[1] - ymax);
	}
	part->position[2] += bw*metric->scaling[2] * rand_type(kernel);
	while((part->position[2] < zmin) || (part->position[2] > zmax)){
		if(part->position[2] < zmin) part->position[2] += 2 * (zmin - part->position[2]);
		else                         part->position[2] -= 2 * (part->position[2] - zmax);
	}
	return 0;
}
int SurfXY_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	double xmin=metric->params[0], xmax=metric->params[1];
	double ymin=metric->params[2], ymax=metric->params[3];
	part->position[0] += bw*metric->scaling[0] * rand_type(kernel);
	while((part->position[0] < xmin) || (part->position[0] > xmax)){
		if(part->position[0] < xmin) part->position[0] += 2 * (xmin - part->position[0]);
		else                         part->position[0] -= 2 * (part->position[0] - xmax);
	}
	part->position[1] += bw*metric->scaling[1] * rand_type(kernel);
	while((part->position[1] < ymin) || (part->position[1] > ymax)){
		if(part->position[1] < ymin) part->position[1] += 2 * (ymin - part->position[1]);
		else                         part->position[1] -= 2 * (part->position[1] - ymax);
	}
	return 0;
}
int SurfR_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	double rho_min=metric->params[0], rho_max=metric->params[1];
	double psi_min=metric->params[2], psi_max=metric->params[3];
	double rho, psi, rho2, psi2;

	rho = sqrt(part->position[0]*part->position[0]+part->position[1]*part->position[1]);
	psi = atan2(part->position[1], part->position[0]);

	rho2 = rho + bw*metric->scaling[0] * rand_type(kernel);
	psi2 = psi + bw*metric->scaling[1]*M_PI/180 * rand_type(kernel);
	
	while((rho2 < rho_min) || (rho2 > rho_max)){
		if(rho2 < rho_min) rho2 += 2 * (rho_min - rho2);
		else               rho2 -= 2 * (rho2 - rho_max);
	}

	while((psi2 < psi_min) || (psi2>psi_max)){
		if(psi2 < psi_min) psi2 += 2 * (psi_min - psi2);
		else               psi2 -= 2 * (psi2 - psi_max);
	}

	part->position[0] = rho2*cos(psi2);
	part->position[1] = rho2*sin(psi2);
	return 0;
}
int SurfR2_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	double rho_min=metric->params[0], rho_max=metric->params[1];
	double psi_min=metric->params[2], psi_max=metric->params[3];
	double rho, psi, rho2, psi2;

	rho_min *= rho_min;
	rho_max *= rho_max;

	rho = part->position[0]*part->position[0]+part->position[1]*part->position[1];
	psi = atan2(part->position[1], part->position[0]);

	rho2 = rho + bw*metric->scaling[0] * rand_type(kernel);
	psi2 = psi + bw*metric->scaling[1]*M_PI/180 * rand_type(kernel);
	
	while((rho2 < rho_min) || (rho2 > rho_max)){
		if(rho2 < rho_min) rho2 += 2 * (rho_min - rho2);
		else               rho2 -= 2 * (rho2 - rho_max);
	}

	while((psi2 < psi_min) || (psi2>psi_max)){
		if(psi2 < psi_min) psi2 += 2 * (psi_min - psi2);
		else               psi2 -= 2 * (psi2 - psi_max);
	}

	part->position[0] = sqrt(rho2)*cos(psi2);
	part->position[1] = sqrt(rho2)*sin(psi2);
	return 0;
}
int SurfCircle_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	double rho_min=metric->params[0], rho_max=metric->params[1];
	double psi_min=metric->params[2], psi_max=metric->params[3];
	double x, y, rho2=0, psi2=0;

	x = part->position[0] + bw*metric->scaling[0] * rand_type(kernel);
	y = part->position[1] + bw*metric->scaling[1] * rand_type(kernel);
	
	rho2 = sqrt(x * x + y * y);
	psi2 = atan2(y, x);

	while((rho2 < rho_min) || (rho2 > rho_max)){
		if(rho2 < rho_min) rho2 += 2 * (rho_min - rho2);
		else               rho2 -= 2 * (rho2 - rho_max);
	}

	while((psi2 < psi_min) || (psi2>psi_max)){
		if(psi2 < psi_min) psi2 += 2 * (psi_min - psi2);
		else               psi2 -= 2 * (psi2 - psi_max);
	}
	part->position[0] = rho2*cos(psi2);
	part->position[1] = rho2*sin(psi2);
	return 0;
}
int Guide_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	double x=part->position[0], y=part->position[1], z=part->position[2];
	double dx=part->direction[0], dy=part->direction[1], dz=part->direction[2], dx2, dz2;
	double xwidth=metric->params[0], yheight=metric->params[1], zmax=metric->params[2], rcurv=metric->params[3];
	double t, mu, mu2, phi;
	int cont=0, mirror;
	if(rcurv != 0){ // Transform to curved guide variables
		double r = sqrt((rcurv+x)*(rcurv+x) + z*z);
		x = copysign(1, rcurv) * r - rcurv; z = fabs(rcurv) * asin(z / r);
		dx2 = dx; dz2 = dz;
		dx = dx2*cos(z/rcurv) + dz2*sin(z/rcurv); dz = -dx2*sin(z/rcurv) + dz2*cos(z/rcurv);
	}
	// Transform from (x,y,z,dx,dy,dz) to (z,t,mu,phi)
	if     ((y/yheight > -x/xwidth) && (y/yheight <  x/xwidth)) mirror=0; // mirror x pos
	else if((y/yheight >  x/xwidth) && (y/yheight > -x/xwidth)) mirror=1; // mirror y pos
	else if((y/yheight < -x/xwidth) && (y/yheight >  x/xwidth)) mirror=2; // mirror x neg
	else                                                        mirror=3; // mirror y neg
	switch(mirror){
		case 0:
			t = 0.5*yheight + y;
			mu = dx; phi = atan2(-dy, dz);
			break;
		case 1:
			t = 1.0*yheight + 0.5*xwidth - x;
			mu = dy; phi = atan2(dx, dz);
			break;
		case 2:
			t = 1.5*yheight + 1.0*xwidth - y;
			mu = -dx; phi = atan2(dy, dz);
			break;
		case 3:
			t = 2.0*yheight + 1.5*xwidth + x;
			mu = -dy; phi = atan2(-dx, dz);
			break;
	}
	// Perturb
	z += bw*metric->scaling[0] * rand_type(kernel);
	while((z < 0) || (z > zmax)){
		if(z < 0) z *= -1;
		else      z -= 2 * (z - zmax);
	}
	t += bw*metric->scaling[1] * rand_type(kernel);
	switch(mirror){ // Avoid perturbation to change mirror
		case 0:
			while((t < 0) || (t > yheight)){
				if(t < 0) t *= -1;              
				else      t -= 2 * (t - yheight);          
			}
			break;
		case 1:
			while((t < yheight) || (t > yheight+xwidth)){
				if(t < yheight) t += 2 * (yheight - t);        
				else            t -= 2 * (t - (yheight+xwidth)); 
			}
			break;
		case 2:
			while((t < yheight+xwidth) || (t > 2*yheight+xwidth)){
				if(t < yheight+xwidth) t += 2 * (yheight+xwidth - t); 
				else                   t -= 2 * (t - (2*yheight+xwidth)); 
			}
			break;
		case 3:
			while((t < 2*yheight+xwidth) || (t > 2*yheight+2*xwidth)){
				if(t < 2*yheight+xwidth) t += 2 * (2*yheight+xwidth - t); 
				else                     t -= 2 * (t - (2*yheight+2*xwidth)); 
			}
			break;
	}
	if(isinf(metric->scaling[2])) mu = -1 + 2.*rand()/RAND_MAX;
	else{
		mu2 = mu + bw*metric->scaling[2] * rand_type(kernel);
		if(mu >= 0){
			while((mu2 < 0) || (mu2 > 1)){
				if(mu2 < 0) mu2 *= -1;
				else        mu2 -= 2 * (mu2 - 1);
			}
		}
		else{
			while((mu2 < -1) || (mu2 > 0)){
				if(mu2 < -1) mu2 += 2 * (-1 - mu2);
				else         mu2 *= -1;
			}
		}
		mu = mu2;
	}
	if(isinf(metric->scaling[3])) phi = 2.*M_PI*rand()/RAND_MAX;
	else phi += bw*metric->scaling[3]*M_PI/180 * rand_type(kernel);
	// Antitransform from (z,t,theta_n,theta_t) to (x,y,z,dx,dy,dz)
	switch(mirror){
		case 0:
			x =  xwidth/2; y =  t - 0.5*yheight;
			dx = mu; dz = sqrt(1-mu*mu)*cos(phi); dy = -sqrt(1-mu*mu)*sin(phi);
			break;
		case 1:
			y =  yheight/2; x = -t + yheight + 0.5*xwidth;
			dy = mu; dz = sqrt(1-mu*mu)*cos(phi); dx = sqrt(1-mu*mu)*sin(phi);
			break;
		case 2:
			x = -xwidth/2; y = -t + 1.5*yheight + xwidth;
			dx = -mu; dz = sqrt(1-mu*mu)*cos(phi); dy = sqrt(1-mu*mu)*sin(phi);
			break;
		case 3:
			y = -yheight/2; x =  t - 2*yheight - 1.5*xwidth;
			dy = -mu; dz = sqrt(1-mu*mu)*cos(phi); dx = -sqrt(1-mu*mu)*sin(phi);
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

void _vMF_perturb(double bw, double* dx, double* dy, double* dz){
	double xi = (double)rand()/RAND_MAX;
	double w = 1;
	w += bw*bw * log(xi+(1-xi)*exp(-2/(bw*bw)));
	double phi = 2.*M_PI*rand()/RAND_MAX;
	double uv = sqrt(1-w*w), u = uv*cos(phi), v = uv*sin(phi);
	double dx0=*dx, dy0=*dy, dz0=*dz;
	if(dz0 > 0){
		*dx = u*dz0 + w*dx0 - (v*dx0-u*dy0)*dy0/(1+dz0);
		*dy = v*dz0 + w*dy0 + (v*dx0-u*dy0)*dx0/(1+dz0);
	}
	else{
		*dx = u*dz0 + w*dx0 + (v*dx0-u*dy0)*dy0/(1-dz0);
		*dy = v*dz0 + w*dy0 - (v*dx0-u*dy0)*dx0/(1-dz0);
	}
	*dz = w*dz0 - u*dx0 - v*dy0;
}

int Isotrop_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	if(isinf(bw*metric->scaling[0])){ // Generate isotropic direction
		part->direction[2] = -1 + 2.*rand()/RAND_MAX;
		double dxy = sqrt(1-part->direction[2]*part->direction[2]);
		double phi = 2.*M_PI*rand()/RAND_MAX;
		part->direction[0] = dxy*cos(phi);
		part->direction[1] = dxy*sin(phi);
	}
	else if(bw*metric->scaling[0] > 0){ // Perturb following von Mises-Fischer distribution
		double dx=part->direction[0], dy=part->direction[1], dz=part->direction[2];
		_vMF_perturb(bw*metric->scaling[0], &part->direction[0], &part->direction[1], &part->direction[2]);
		int keepx = (int)metric->params[0], keepy = (int)metric->params[1], keepz = (int)metric->params[2];
		if(keepx && (dx*part->direction[0] < 0)) part->direction[0] *= -1;
		if(keepy && (dy*part->direction[1] < 0)) part->direction[1] *= -1;
		if(keepz && (dz*part->direction[2] < 0)) part->direction[2] *= -1;
	}
	return 0;
}

int Polar_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	double theta, theta2, phi;
	int cont=0;
	theta = acos(part->direction[2]);
	phi   = atan2(part->direction[1], part->direction[0]);
	theta2 = theta + bw*metric->scaling[0]*M_PI/180 * rand_type(kernel);
	if(cos(theta2)*cos(theta) < 0) theta += 2 * (M_PI/2 - theta);
	theta = theta2;
	phi   += bw*metric->scaling[1]*M_PI/180 * rand_type(kernel);
	part->direction[0] = sin(theta) * cos(phi);
	part->direction[1] = sin(theta) * sin(phi);
	part->direction[2] = cos(theta);
	return 0;
}

int PolarMu_perturb(const Metric* metric, mcpl_particle_t* part, double bw, char kernel){
	double mu, mu2, phi;
	int cont=0;
	mu = part->direction[2];
	phi   = atan2(part->direction[1], part->direction[0]);
	mu2 = mu + bw*metric->scaling[0] * rand_type(kernel);
	if(mu >= 0){
		while((mu2 < 0) || (mu2 > 1)){
			if(mu2 < 0)  mu2 *= -1;
			if(mu2 > 1)  mu2 -= 2 * (mu2 - 1);
		}
	}
	else{
		while((mu2 < -1) || (mu2 > 0)){
			if(mu2 < -1)  mu2 += 2 * (-1 - mu2);
			if(mu2 > 0)  mu2 *= -1;
		}
	}
	mu = mu2;
	phi   += bw*metric->scaling[1]*M_PI/180 * rand_type(kernel);
	part->direction[0] = sqrt(1-mu*mu) * cos(phi);
	part->direction[1] = sqrt(1-mu*mu) * sin(phi);
	part->direction[2] = mu;
	return 0;
}
