#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "ksource.h"


PList* PList_create(char pt, int ord, char** filenames, ReadFun* read, double* trasl, double* rot, int switch_x2z){
	PList* plist = (PList*)malloc(sizeof(PList));
	plist->pt = pt;
	plist->ord = ord;
	int i;
	FILE* files[ord];
	for(i=0; i<ord; i++){
		files[i] = NULL;
		if(filenames[i]) if(strlen(filenames[i])){
			if((files[i]=fopen(filenames[i], "r")) == 0){
				printf("Error en PList_create: No se pudo abrir archivo %s\n", filenames[i]);
				return NULL;
			}
		}
	}
	plist->filenames = (char**)malloc(ord * sizeof(char*));
	plist->files = (FILE**)malloc(ord * sizeof(FILE*));
	plist->read = (ReadFun*)malloc(ord * sizeof(ReadFun));
	for(i=0; i<ord; i++){
		plist->filenames[i] = NULL;
		if(files[i]){
			plist->filenames[i] = (char*)malloc(NAME_MAX_LEN*sizeof(char));
			strcpy(plist->filenames[i], filenames[i]);
		}
		plist->files[i] = files[i];
		plist->read[i] = read[i];
	}
	if(trasl){
		plist->trasl = (double*)malloc(3 * sizeof(double));
		for(int i=0; i<3; i++) plist->trasl[i] = trasl[i];
	}
	else plist->trasl = NULL;
	if(rot){
		plist->rot = (double*)malloc(3 * sizeof(double));
		for(int i=0; i<3; i++) plist->rot[i] = rot[i];
	}
	else plist->rot = NULL;
	plist->x2z = switch_x2z;
	PList_next(plist);
	return plist;
}

PList* PList_copy(PList* from){
	PList* plist = (PList*)malloc(sizeof(PList));
	*plist = *from;
	plist->filenames = (char**)malloc(plist->ord * sizeof(char*));
	plist->files = (FILE**)malloc(plist->ord * sizeof(FILE*));
	plist->read = (ReadFun*)malloc(plist->ord * sizeof(ReadFun));
	int i;
	for(i=0; i<plist->ord; i++){
		plist->filenames[i] = NULL;
		if(from->filenames[i]) if(strlen(from->filenames[i])){
			plist->filenames[i] = (char*)malloc(NAME_MAX_LEN*sizeof(char));
			strcpy(plist->filenames[i], from->filenames[i]);
		}
		plist->files[i] = fopen(from->filenames[i], "r");
		plist->read[i] = from->read[i];
	}
	if(from->trasl){
		plist->trasl = (double*)malloc(3 * sizeof(double));
		for(int i=0; i<3; i++) plist->trasl[i] = from->trasl[i];
	}
	if(from->rot){
		plist->rot = (double*)malloc(3 * sizeof(double));
		for(int i=0; i<3; i++) plist->rot[i] = from->rot[i];
	}
	return plist;
}

int PList_get(PList* plist, Part* part, double* w){
	*part = plist->part;
	*w = plist->w;
	if(plist->trasl)
		traslv(part->pos, plist->trasl, 0);
	if(plist->rot){
		rotv(part->pos, plist->rot, 0);
		rotv(part->dir, plist->rot, 0);
	}
	if(plist->x2z){
		double x=part->pos[0], y=part->pos[1], z=part->pos[2];
		part->pos[0] = y;
		part->pos[1] = z;
		part->pos[2] = x;
		double dx=part->dir[0], dy=part->dir[1], dz=part->dir[2];
		part->dir[0] = dy;
		part->dir[1] = dz;
		part->dir[2] = dx;
	}
	return 0;
}

int PList_next(PList* plist){
	double w=1, wi;
	Part part = {1, {0,0,0}, {0,0,1}};
	int i, cont=0, ret=0;
	for(i=0; i<plist->ord; i++){
		wi = 1;
		while(cont++ < MAX_SEARCH){
			if(plist->files[i]){  // Si hay archivo, leo una linea
				if(!fgets(plist->line, LINE_MAX_LEN, plist->files[i])){ // Si llego al final, vuelvo al principio
					rewind(plist->files[i]);
					fgets(plist->line, LINE_MAX_LEN, plist->files[i]);
					ret = 1;
				}
			}
			if(plist->read[i](plist->line, &part, &wi)) break; // Extraigo linea
		}
		if(cont > MAX_SEARCH) printf("Warning en PList_next: No se encontro particula\n");
		w *= wi;
	}
	plist->part = part;
	plist->w = w;
	return ret;
}

void PList_destroy(PList* plist){
	int i;
	for(i=0; i<plist->ord; i++){
		free(plist->filenames[i]);
		fclose(plist->files[i]);
	}
	free(plist->filenames);
	free(plist->files);
	free(plist->read);
	free(plist->trasl);
	free(plist->rot);
	free(plist);
}

PList* PListSimple_create(char pt, char* filename, ReadFun read, double* trasl, double* rot, int switch_x2z){
	char* pfilename[] = {filename};
	ReadFun pread[] = {read};
	return PList_create(pt, 1, pfilename, pread, trasl, rot, switch_x2z);
}


PList* PListSepVar_create(char pt, char* filenames[3], ReadFun read[3], double* trasl, double* rot, int switch_x2z){
	return PList_create(pt, 3, filenames, read, trasl, rot, switch_x2z);
}


int PTRAC_read(char* line, Part* part, double* w){
	double aux;
	int nreaded = sscanf(line, "%le %le %le %le %le %le %le %le %le %le",
		&part->pos[0], &part->pos[1], &part->pos[2], &part->dir[0], &part->dir[1], &part->dir[2], &part->E, w, &aux, &aux);
	return (nreaded == 9);
}

int T4stock_read(char* line, Part* part, double* w){
	if(strncmp(line,"NEUTRON",7)==0 || strncmp(line,"PHOTON",6)==0){
		sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf",
			&part->E, &part->pos[0], &part->pos[1], &part->pos[2], &part->dir[0], &part->dir[1], &part->dir[2], w);
		return 1;
	}
	return 0;
}

int SSV_read(char* line, Part* part, double* w){
	double aux;
	int nreaded = sscanf(line, "%le %le %le %le %le %le %le %le %le",
		&part->E, &part->pos[0], &part->pos[1], &part->pos[2], &part->dir[0], &part->dir[1], &part->dir[2], w, &aux);
	return (nreaded == 8);
}

int Decay_read(char* line, Part* part, double* w){
	double aux;
	int nreaded = sscanf(line, "%lf,%lf,%lf", &part->E, &aux, w);
	part->E /= 1000;
	return (nreaded == 3);
}
int SSVtally_read(char* line, Part* part, double* w){
	double aux;
	int nreaded = sscanf(line, "%le %le %le %le %le",
		&part->pos[0], &part->pos[1], &part->pos[2], w, &aux);
	return (nreaded == 4);
}
int Isotrop_read(char* line, Part* part, double* w){
	part->dir[2] = -1 + 2.*rand()/RAND_MAX;
	double dxy = sqrt(1-part->dir[2]*part->dir[2]);
	double phi = 2.*M_PI*rand()/RAND_MAX;
	part->dir[0] = dxy*cos(phi);
	part->dir[1] = dxy*sin(phi);
	return 1;
}
