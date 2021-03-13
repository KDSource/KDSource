#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "ksource.h"


Metric* Metric_create(int dim, double* bw, PerturbFun perturb, int n_gp, double* geom_par){
	Metric* metric = (Metric*)malloc(sizeof(Metric));
	int i;
	metric->dim = dim;
	metric->bw = (double*)malloc(dim*sizeof(double));
	if(bw) for(i=0; i<dim; i++) metric->bw[i] = bw[i];
	else for(i=0; i<dim; i++) metric->bw[i] = 0;
	metric->perturb = perturb;
	metric->n_gp = n_gp;
	metric->geom_par = (double*)malloc(n_gp * sizeof(double));
	for(i=0; i<n_gp; i++) metric->geom_par[i] = geom_par[i];
	return metric;
}

void Metric_destroy(Metric* metric){
	free(metric->bw);
	free(metric);
}

MetricSepVar* MetricSepVar_create(int ord, Metric** metrics, char* bwfilename, int variable_bw, double* trasl, double* rot){
	MetricSepVar* metric = (MetricSepVar*)malloc(sizeof(MetricSepVar));
	metric->ord = ord;
	int i;
	metric->ms = (Metric**)malloc(ord * sizeof(Metric*));
	for(i=0; i<ord; i++) metric->ms[i] = metrics[i];
	if(bwfilename && strlen(bwfilename)){
		FILE* file_bw;
		if((file_bw=fopen(bwfilename, "r")) == 0){
			printf("Error: No se pudo abrir archivo %s\n", bwfilename);
			return NULL;
		}
		metric->file_bw = file_bw;
		MetricSepVar_next(metric);
		if(!variable_bw){
			fclose(metric->file_bw);
			metric->file_bw = NULL;
		}
	}
	else metric->file_bw = NULL;
	if(trasl){
		metric->trasl = (double*)malloc(3 * sizeof(double));
		for(int i=0; i<3; i++) metric->trasl[i] = trasl[i];
	}
	else metric->trasl = NULL;
	if(rot){
		metric->rot = (double*)malloc(3 * sizeof(double));
		for(int i=0; i<3; i++) metric->rot[i] = rot[i];
	}
	else metric->rot = NULL;
	return metric;
}

int MetricSepVar_perturb(MetricSepVar* metric, Part* part){
	int i, ret=0;
	if(metric->trasl) traslv(part->pos, metric->trasl, 1);
	if(metric->rot){ rotv(part->pos, metric->rot, 1); rotv(part->dir, metric->rot, 1); }
	for(i=0; i<metric->ord; i++) ret += metric->ms[i]->perturb(metric->ms[i], part);
	if(metric->rot){ rotv(part->pos, metric->rot, 0); rotv(part->dir, metric->rot, 0); }
	if(metric->trasl) traslv(part->pos, metric->trasl, 0);
	return ret;
}

int MetricSepVar_next(MetricSepVar* metric){
	if(metric->file_bw){
		int i, j, dim=0, readed=0;
		for(i=0; i<metric->ord; i++) dim += metric->ms[i]->dim;
		for(i=0; i<metric->ord; i++)
			for(j=0; j<metric->ms[i]->dim; j++)
				readed += fscanf(metric->file_bw, "%lf", &metric->ms[i]->bw[j]);
		if(readed < dim){ // Si llego al final del archivo, vuelvo al inicio y releo
			if(readed != -dim) printf("Warning: Archivo de anchos de banda no tiene el formato esperado\n");
			rewind(metric->file_bw);
			readed = 0;
			for(i=0; i<metric->ord; i++)
				for(j=0; j<metric->ms[i]->dim; j++)
					readed += fscanf(metric->file_bw, "%lf", &metric->ms[i]->bw[j]);
		}
		if(readed < dim){
			printf("Error: No se pudo leer ancho de banda\n");
			return 1;
		}
	}
	return 0;
}

void MetricSepVar_destroy(MetricSepVar* metric){
	free(metric->ms);
	free(metric->trasl);
	free(metric->rot);
	free(metric);
}


int E_perturb(Metric* metric, Part* part){
	part->E += metric->bw[0] * rand_norm();
	if(part->E < E_MIN) part->E = E_MIN;
	return 0;
}
int Let_perturb(Metric* metric, Part* part){
	part->E *= exp(metric->bw[0] * rand_norm());
	if(part->E < E_MIN) part->E = E_MIN;
	return 0;
}

int Vol_perturb(Metric* metric, Part* part){
	part->pos[0] += metric->bw[0] * rand_norm();
	part->pos[1] += metric->bw[1] * rand_norm();
	part->pos[2] += metric->bw[2] * rand_norm();
	if(part->pos[0] < metric->geom_par[0]) part->pos[0] = metric->geom_par[0];
	else if(part->pos[0] > metric->geom_par[1]) part->pos[0] = metric->geom_par[1];
	if(part->pos[1] < metric->geom_par[2]) part->pos[1] = metric->geom_par[2];
	else if(part->pos[1] > metric->geom_par[3]) part->pos[1] = metric->geom_par[3];
	if(part->pos[2] < metric->geom_par[4]) part->pos[2] = metric->geom_par[4];
	else if(part->pos[2] > metric->geom_par[5]) part->pos[2] = metric->geom_par[5];
	return 0;
}
int SurfXY_perturb(Metric* metric, Part* part){
	part->pos[0] += metric->bw[0] * rand_norm();
	part->pos[1] += metric->bw[1] * rand_norm();
	if(part->pos[0] < metric->geom_par[0]) part->pos[0] = metric->geom_par[0];
	else if(part->pos[0] > metric->geom_par[1]) part->pos[0] = metric->geom_par[1];
	if(part->pos[1] < metric->geom_par[2]) part->pos[1] = metric->geom_par[2];
	else if(part->pos[1] > metric->geom_par[3]) part->pos[1] = metric->geom_par[3];
	return 0;
}
int Guide_perturb(Metric* metric, Part* part){
	double t, x=part->pos[0], y=part->pos[1], z=part->pos[2];
	double xwidth=metric->geom_par[0], yheight=metric->geom_par[1], zlength=metric->geom_par[2], rcurv=metric->geom_par[3];
	if(rcurv != 0){
		double r = sqrt((rcurv+x)*(rcurv+x) + z*z);
		x = copysign(1, rcurv) * r - rcurv;
		z = fabs(rcurv) * asin(z / r);
	}
	if     ((y/yheight > -x/xwidth) && (y/yheight <  x/xwidth)) t = 0.5*yheight + y;              // espejo x pos
	else if((y/yheight >  x/xwidth) && (y/yheight > -x/xwidth)) t = 1.0*yheight + 0.5*xwidth - x; // espejo y pos
	else if((y/yheight < -x/xwidth) && (y/yheight >  x/xwidth)) t = 1.5*yheight + 1.0*xwidth - y; // espejo x neg
	else                                                        t = 2.0*yheight + 1.5*xwidth + x; // espejo y neg
	z += metric->bw[0] * rand_norm();
	t += metric->bw[1] * rand_norm();
	if(z < 0) z = 0;
	else if(z > zlength) z = zlength;
	while(t < 0) t += 2*(xwidth+yheight);
	while(t > 2*(xwidth+yheight)) t -= 2*(xwidth+yheight);
	if     (t < yheight)                                    {x =  xwidth /2; y =  t - 0.5*yheight;             } // espejo x pos
	else if((t > yheight)        && (t <   yheight+xwidth)) {y =  yheight/2; x = -t + 1.0*yheight + 0.5*xwidth;} // espejo y pos
	else if((t > yheight+xwidth) && (t < 2*yheight+xwidth)) {x = -xwidth /2; y = -t + 1.5*yheight + 1.0*xwidth;} // espejo x neg
	else                                                    {y = -yheight/2; x =  t - 2.0*yheight - 1.5*xwidth;} // espejo y neg
	if(rcurv != 0){
		double r = (rcurv + x) * copysign(1, rcurv);
		double ang = z / fabs(rcurv);
		x = copysign(1, rcurv) * r * cos(ang) - rcurv;
		z = r * sin(ang);
	}
	part->pos[0] = x;
	part->pos[1] = y;
	part->pos[2] = z;
	return 0;
}

int Isotrop_perturb(Metric* metric, Part* part){
	if(metric->bw[0] == INFINITY){
		part->dir[2] = -1 + 2.*rand()/RAND_MAX;
		double dxy = sqrt(1-part->dir[2]*part->dir[2]);
		double phi = 2.*M_PI*rand()/RAND_MAX;
		part->dir[0] = dxy*cos(phi);
		part->dir[1] = dxy*sin(phi);
	}
	else{
		double xi = (double)rand()/RAND_MAX;
		double w = 1;
		if(metric->bw[2] > BW_MIN) w += metric->bw[2]*metric->bw[2] * log(xi+(1-xi)*exp(-2/(metric->bw[2]*metric->bw[2])));
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
	}
	return 0;
}

int Polar_perturb(Metric* metric, Part* part){
	double theta, phi;
	theta = acos(part->dir[2]);
	phi   = atan2(part->dir[1], part->dir[0]);
	theta += metric->bw[0] * rand_norm();
	phi   += metric->bw[1] * rand_norm();
	part->dir[0] = sin(theta) * cos(phi);
	part->dir[1] = sin(theta) * sin(phi);
	part->dir[2] = cos(theta);
	return 0;
}
