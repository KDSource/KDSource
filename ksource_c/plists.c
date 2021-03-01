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
	for(i=0; i<ord; i++)
		if((files[i]=fopen(filenames[i], "r")) == 0){
			printf("Error: No se pudo abrir archivo %s\n", filenames[i]);
			return NULL;
		}
	plist->files = (FILE**)malloc(ord * sizeof(FILE*));
	plist->read = (ReadFun*)malloc(ord * sizeof(ReadFun));
	for(i=0; i<ord; i++){
		plist->files[i] = files[i];
		plist->read[i] = read[i];
	}
	if(trasl){
		plist->trasl = (double*)malloc(3 * sizeof(double));
		for(int i=0; i<3; i++) plist->trasl[i] = trasl[i];
	}
	if(rot){
		plist->rot = (double*)malloc(3 * sizeof(double));
		for(int i=0; i<3; i++) plist->rot[i] = rot[i];
	}
	plist->x2z = switch_x2z;
	plist->part = (Part*)malloc(sizeof(Part));
	PList_next(plist);
	return plist;
}

int PList_get(PList* plist, Part* part, double* w){
	*part = (*plist->part);
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
	Part part;
	int i;
	int part_updated;
	for(i=0; i<plist->ord; i++){
		part_updated = 0;
		while(!part_updated){
			if(plist->files[i]){ // Si hay archivo, leo una linea
				if(!fgets(plist->line, LINE_MAX_LEN, plist->files[i])){ // Si llego al final, vuelvo al principio
					rewind(plist->files[i]);
					fgets(plist->line, LINE_MAX_LEN, plist->files[i]);
				}
			}
			if(plist->read[i](plist->line, &part, &wi)) part_updated = 1; // Extraigo linea
		}
		w *= wi;
	}
	*plist->part = part;
	plist->w = w;
	return 0;
}

void PList_destroy(PList* plist){
	int i;
	for(i=0; i<plist->ord; i++) fclose(plist->files[i]);
	free(plist->files);
	free(plist->read);
	free(plist->trasl);
	free(plist->rot);
	free(plist->part);
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

int Tripoli_read_part(char* line, Part* part, double* w){
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
	sscanf(line, "%lf,%lf,%lf", &part->E, &aux, w);
	return 1;
}

int Tripoli_read_pos(char* line, Part* part, double* w){
	double aux;
	if(strncmp(line,"NEUTRON",7)==0 || strncmp(line,"PHOTON",6)==0){
		sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf",
			&aux, &part->pos[0], &part->pos[1], &part->pos[2], &aux, &aux, &aux, w);
		return 1;
	}
	return 0;
}

int Isotrop_read(char* line, Part* part, double* w){
	*w = 1;
	part->dir[2] = -1 + 2.*rand()/RAND_MAX;
	double dxy = sqrt(1-part->dir[2]*part->dir[2]);
	double phi = 2.*M_PI*rand()/RAND_MAX;
	part->dir[0] = dxy*cos(phi);
	part->dir[1] = dxy*sin(phi);
}
