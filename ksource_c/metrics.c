#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "ksource.h"


Metric* Metric_create(int n, int* dims, double** bw, char* bwfilename, int variable_bw, PerturbFun* perturb,
		double trasl[3], double rot[3], double** geom_par){
	Metric* metric = (Metric*)malloc(sizeof(Metric));
	int i, j;
	metric->n = n;
	metric->dims = dims;
	metric->bw = (double**)malloc(n*sizeof(double*));
	for(i=0; i<n; i++) metric->bw[i] = (double*)malloc(dims[i]*sizeof(double));
	if(bw) for(i=0; i<n; i++) for(j=0; j<dims[i]; j++) metric->bw[i][j] = bw[i][j];
	if(bwfilename && strlen(bwfilename)){
		metric->file_bw = fopen(bwfilename, "r");
		Metric_next(metric);
		if(!variable_bw){
			fclose(metric->file_bw);
			metric->file_bw = NULL;
		}
		metric->variable_bw = variable_bw;
	}
	else{
		metric->file_bw = NULL;
		metric->variable_bw = 0;
	}
	metric->perturb = perturb;
	metric->trasl = trasl;
	metric->rot = rot;
	metric->geom_par = geom_par;
	return metric;
}

int Metric_perturb(Metric* metric, Part* part){
	int i, ret=1;
	if(metric->trasl) traslv(part->pos, metric->trasl, 1);
	if(metric->rot){ rotv(part->pos, metric->rot, 1); rotv(part->dir, metric->rot, 1); }
	for(i=0; i<metric->n; i++)
		if(metric->perturb[i](metric, part)) return 1;
	if(metric->rot){ rotv(part->pos, metric->rot, 0); rotv(part->dir, metric->rot, 0); }
	if(metric->trasl) traslv(part->pos, metric->trasl, 0);
	return 0;
}

int Metric_next(Metric* metric){
	if(metric->variable_bw){
		int i, j, ret=1;
		for(i=0; i<metric->n; i++)
			for(j=0; j<metric->dims[i]; j++)
				if(!fscanf(metric->file_bw, "%lf", &metric->bw[i][j])){
					rewind(metric->file_bw);
					if(!fscanf(metric->file_bw, "%lf", &metric->bw[i][j])) return 1;
				}
	}
	return 0;
}

void Metric_destroy(Metric* metric){
	int i;
	for(i=0; i<metric->n; i++) free(metric->bw[i]);
	free(metric->bw);
	if(metric->file_bw) fclose(metric->file_bw);
	free(metric);
}


Metric* MetricSimple_create(int dim, double* bw, char* bwfilename, int variable_bw, PerturbFun perturb,
		double trasl[3], double rot[3], double* geom_par){
	int* pdim = (int*)malloc(sizeof(int));
	double** pbw = (double**)malloc(sizeof(double*));
	*pbw = bw;
	PerturbFun* pperturb = (PerturbFun*)malloc(sizeof(PerturbFun));
	*pperturb = perturb;
	double** pgeom_par = (double**)malloc(sizeof(double*));
	*pgeom_par = geom_par;
	return Metric_create(1, pdim, pbw, bwfilename, variable_bw, pperturb, trasl, rot, pgeom_par);
}

void MetricSimple_destroy(Metric* metric){
	free(metric->dims);
	free(metric->bw);
	free(metric->perturb);
	free(metric->geom_par);
}

Metric* MetricSepVar_create(int dims[3], double* bw[3], char* bwfilename, int variable_bw, PerturbFun perturb[3],
		double trasl[3], double rot[3], double** geom_par){
	int* pdims = (int*)malloc(3*sizeof(int));
	pdims[0] = dims[0];
	pdims[1] = dims[1];
	pdims[2] = dims[2];
	double** pbw = (double**)malloc(3*sizeof(double*));
	pbw[0] = bw[0];
	pbw[1] = bw[1];
	pbw[2] = bw[2];
	PerturbFun* pperturb = (PerturbFun*)malloc(3*sizeof(PerturbFun));
	pperturb[0] = perturb[0];
	pperturb[1] = perturb[1];
	pperturb[2] = perturb[2];
	double** pgeom_par = (double**)malloc(3*sizeof(double*));
	if(geom_par){
		pgeom_par[0] = geom_par[0];
		pgeom_par[1] = geom_par[1];
		pgeom_par[2] = geom_par[2];
	}
	else pgeom_par[0] = pgeom_par[1] = pgeom_par[2] = NULL;
	return Metric_create(3, dims, pbw, bwfilename, variable_bw, pperturb, trasl, rot, pgeom_par);
}

void MetricSepVar_destroy(Metric* metric){
	free(metric->dims);
	free(metric->bw);
	free(metric->perturb);
	free(metric->geom_par);
}


int E_perturb(Metric* metric, Part* part){
	part->E += *metric->bw[0] * rand_norm();
	return 0;
}
int Let_perturb(Metric* metric, Part* part){
	part->E *= exp(*metric->bw[0] * rand_norm());
	return 0;
}

int Vol_perturb(Metric* metric, Part* part){
	part->pos[0] += metric->bw[1][0] * rand_norm();
	part->pos[1] += metric->bw[1][1] * rand_norm();
	part->pos[2] += metric->bw[1][2] * rand_norm();
	return 0;
}
int SurfXY_perturb(Metric* metric, Part* part){
	part->pos[0] += metric->bw[1][0] * rand_norm();
	part->pos[1] += metric->bw[1][1] * rand_norm();
	return 0;
}
int Guide_perturb(Metric* metric, Part* part){
	double x=part->pos[0], y=part->pos[1], z=part->pos[2];
	double xwidth=metric->geom_par[1][0], yheight=metric->geom_par[1][1], rcurv=metric->geom_par[1][2];
	if(rcurv != 0){
		double r = sqrt((rcurv+x)*(rcurv+x) + z*z);
		x = r - rcurv;
		z = rcurv * atan2(z, rcurv+x);
	}
	z += metric->bw[1][0] * rand_norm();
	if( (y/yheight > -x/xwidth) && (y/yheight <  x/xwidth) ){ // espejo x pos
		y += metric->bw[1][1] * rand_norm();
		if(y >  yheight/2){ x -=  y-yheight/2; y =  yheight/2; }
		if(y < -yheight/2){ x -= -y-yheight/2; y = -yheight/2; }
	}
	if( (y/yheight >  x/xwidth) && (y/yheight > -x/xwidth) ){ // espejo y pos
		x += metric->bw[1][1] * rand_norm();
		if(x >  xwidth/2){ y -=  x-xwidth/2; x =  xwidth/2; }
		if(x < -xwidth/2){ y -= -x-xwidth/2; x = -xwidth/2; }
	}
	if( (y/yheight < -x/xwidth) && (y/yheight >  x/xwidth) ){ // espejo x neg
		y += metric->bw[1][1] * rand_norm();
		if(y >  yheight/2){ x +=  y-yheight/2; y =  yheight/2; }
		if(y < -yheight/2){ x += -y-yheight/2; y = -yheight/2; }
	}
	if( (y/yheight <  x/xwidth) && (y/yheight < -x/xwidth) ){ // espejo y neg
		x += metric->bw[1][1] * rand_norm();
		if(x >  xwidth/2){ y +=  x-xwidth/2; x =  xwidth/2; }
		if(x < -xwidth/2){ y += -x-xwidth/2; x = -xwidth/2; }
	}
	if(rcurv != 0){
		double r = rcurv + x;
		double ang = tan(z / rcurv);
		x = r * cos(ang) - rcurv;
		z = r * sin(ang);
	}
	part->pos[0] = x;
	part->pos[1] = y;
	part->pos[2] = z;
	return 0;
}

int Isotrop_perturb(Metric* metric, Part* part){
	double xi = (double)rand()/RAND_MAX;
	double w = 1;
	if(*metric->bw[2] > BW_MIN) w += *metric->bw[2]**metric->bw[2] * log(xi+(1-xi)*exp(-2/(*metric->bw[2]**metric->bw[2])));
	double phi = 2.*M_PI*rand()/RAND_MAX;
	double uv = sqrt(1-w*w), u = uv*cos(phi), v = uv*sin(phi);
	double x=part->dir[0], y=part->dir[1], z=part->dir[2];
	if(part->dir[2] > 0){
		part->dir[0] = u*z + w*x - (v*x-u*y)*y/(1+z);
		part->dir[1] = v*z + w*y + (v*x-u*y)*x/(1+z);
	}
	else{
		part->dir[0] = u*z + w*x + (v*x-u*y)*y/(1-z);
		part->dir[1] = v*z + w*y - (v*x-u*y)*x/(1-z);
	}
	part->dir[2] = w*z - u*x - v*y;
	return 0;
}

int Polar_perturb(Metric* metric, Part* part){
	double theta, phi;
	theta = acos(part->dir[2]);
	phi   = atan2(part->dir[1], part->dir[0]);
	theta += metric->bw[2][0] * rand_norm();
	phi   += metric->bw[2][1] * rand_norm();
	part->dir[0] = sin(theta) * cos(phi);
	part->dir[1] = sin(theta) * sin(phi);
	part->dir[2] = cos(theta);
	return 0;
}
