#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "ksource.h"


PList* PList_create(char pt, double* trasl, double* rot, int switch_x2z, int n, FILE** files, ReadFun* read){
	PList* plist = (PList*)malloc(sizeof(PList));
	plist->pt = pt;
	plist->trasl = trasl;
	plist->rot = rot;
	plist->x2z = switch_x2z;
	plist->line = (char*)malloc(LINE_MAX_LEN*sizeof(char));
	plist->n = n;
	plist->files = files;
	plist->part = (Part*)malloc(sizeof(Part));
	plist->read = read;
	PList_next(plist);
	return plist;
}

int PList_get(PList* plist, Part* part, double* w){
	*part = (*plist->part);
	*w = plist->w;
	if(plist->trasl)
		traslv(part->pos, plist->trasl);
	if(plist->rot){
		rotv(part->pos, plist->rot);
		rotv(part->dir, plist->rot);
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
	for(i=0; i<plist->n; i++){
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
	// Aplicar trasl, rot
	return 0;
}

void PList_destroy(PList* plist){
	free(plist->line);
	free(plist->part);
	free(plist);
}


PList* PListSimple_create(char pt, double* trasl, double* rot, int switch_x2z, char* filename, ReadFun read){
	FILE** files = (FILE**)malloc(sizeof(FILE*));
	*files = fopen(filename, "r");
	ReadFun* pread = (ReadFun*)malloc(sizeof(ReadFun));
	*pread = read;
	return PList_create(pt, trasl, rot, switch_x2z, 1, files, pread);
}

void PListSimple_destroy(PList* plist){
	fclose(*plist->files);
	free(plist->files);
	free(plist->read);
	PList_destroy(plist);
}


PList* PListSepVar_create(char pt, double* trasl, double* rot, int switch_x2z, char* filename[3], ReadFun read[3]){
	FILE** files = (FILE**)malloc(3*sizeof(FILE*));
	files[0] = fopen(filename[0], "r");
	files[1] = fopen(filename[1], "r");
	files[2] = fopen(filename[2], "r");
	ReadFun* pread = (ReadFun*)malloc(3*sizeof(ReadFun));
	pread[0] = read[0];
	pread[1] = read[1];
	pread[2] = read[2];
	return PList_create(pt, trasl, rot, switch_x2z, 3, files, pread);
}

void PListSepVar_destroy(PList* plist){
	int i;
	for(i=0; i<3; i++) fclose(plist->files[i]);
	free(plist->files);
	free(plist->read);
	PList_destroy(plist);
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
