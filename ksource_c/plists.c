#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "ksource.h"


PList* PList_create(char pt, const char* filename, const double* trasl, const double* rot, int switch_x2z){
	PList* plist = (PList*)malloc(sizeof(PList));
	plist->pt = pt;
	mcpl_file_t file = mcpl_open_file(filename);
	plist->filename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
	strcpy(plist->filename, filename);
	plist->file = file;
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

PList* PList_copy(const PList* from){
	PList* plist = (PList*)malloc(sizeof(PList));
	*plist = *from;
	plist->filename = (char*)malloc(NAME_MAX_LEN*sizeof(char));
	strcpy(plist->filename, from->filename);
	plist->file = mcpl_open_file(plist->filename);
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

int PList_get(const PList* plist, mcpl_particle_t* part){
	*part = plist->part;
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

int PList_next(PList* plist){
	const mcpl_particle_t* ppart = mcpl_read(plist->file);
	if(!ppart) return 1;
	plist->part = *ppart;
	return 0;
}

void PList_destroy(PList* plist){
	free(plist->filename);
	mcpl_close_file(plist->file);
	free(plist->trasl);
	free(plist->rot);
	free(plist);
}
