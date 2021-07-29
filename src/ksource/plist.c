#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "ksource.h"


void KS_error(const char* msg);
void KS_end(const char* msg);

PList* PList_create(char pt, const char* filename, const double* trasl, const double* rot, int switch_x2z){
	PList* plist = (PList*)malloc(sizeof(PList));
	plist->pt = pt;
	mcpl_file_t file = mcpl_open_file(filename);
	plist->filename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
	strcpy(plist->filename, filename);
	plist->file = file;
	int i;
	if(trasl){
		plist->trasl = (double*)malloc(3 * sizeof(double));
		for(i=0; i<3; i++) plist->trasl[i] = trasl[i];
	}
	else plist->trasl = NULL;
	if(rot){
		plist->rot = (double*)malloc(3 * sizeof(double));
		for(i=0; i<3; i++) plist->rot[i] = rot[i];
	}
	else plist->rot = NULL;
	plist->x2z = switch_x2z;
	plist->part = NULL;
	PList_next(plist, 1);
	return plist;
}

PList* PList_copy(const PList* from){
	PList* plist = (PList*)malloc(sizeof(PList));
	*plist = *from;
	plist->filename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
	strcpy(plist->filename, from->filename);
	plist->file = mcpl_open_file(plist->filename);
	int i;
	if(from->trasl){
		plist->trasl = (double*)malloc(3 * sizeof(double));
		for(i=0; i<3; i++) plist->trasl[i] = from->trasl[i];
	}
	if(from->rot){
		plist->rot = (double*)malloc(3 * sizeof(double));
		for(i=0; i<3; i++) plist->rot[i] = from->rot[i];
	}
	plist->part = NULL;
	PList_next(plist, 1);
	return plist;
}

int PList_get(const PList* plist, mcpl_particle_t* part){
	*part = *plist->part;
	if(plist->trasl)
		traslv(part->position, plist->trasl, 0);
	if(plist->rot){
		rotv(part->position, plist->rot, 0);
		rotv(part->direction, plist->rot, 0);
	}
	if(plist->x2z){
		double x=part->position[0], y=part->position[1], z=part->position[2];
		part->position[0] = y; part->position[1] = z; part->position[2] = x;
		double dx=part->direction[0], dy=part->direction[1], dz=part->direction[2];
		part->direction[0] = dy; part->direction[1] = dz; part->direction[2] = dx;
	}
	return 0;
}

int PList_next(PList* plist, int loop){
	int ret=0, resamples=0;
	while(resamples++ < MAX_RESAMPLES){
		plist->part = mcpl_read(plist->file);
		if(plist->part == NULL){ // Luego del 1er intento fallido rebobino
			if(!loop) KS_end("PList_next: Fin de lista de particulas alcanzado");
			ret = 1;
			mcpl_rewind(plist->file);
			plist->part = mcpl_read(plist->file);
			if(plist->part == NULL)
				KS_error("Error en PList_next: No se pudo obtener particula");
		}
		if(plist->part->weight > 0) return ret;
	}
	KS_error("Error en PList_next: MAX_RESAMPLES alcanzado");
}

void PList_destroy(PList* plist){
	free(plist->filename);
	mcpl_close_file(plist->file);
	free(plist->trasl);
	free(plist->rot);
	free(plist);
}
