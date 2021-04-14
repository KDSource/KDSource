#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "mcpl.h" // "/opt/mcpl/include/mcpl.h"

#ifndef PLISTS_H
#define PLISTS_H


#define MAX_SEARCH 1E5
#define LINE_MAX_LEN 256
#define NAME_MAX_LEN 96

typedef struct PList{
	char pt; // Tipo de particula (n, p, ...)

	char* filename; // Nombres de archivos de tracks (unico o uno por variable)
	mcpl_file_t file; // Archivos de tracks (unico o uno por variable)

	double* trasl; // Traslacion a aplicar
	double* rot; // Rotacion a aplicar
	int x2z; // Si es true, se aplica permutacion x,y,z -> y,z,x

	mcpl_particle_t part; // Buffer para particula
} PList;

PList* PList_create(char pt, const char* filename, const double* trasl, const double* rot, int switch_x2z);
PList* PList_copy(const PList* from);
int PList_get(const PList* plist, mcpl_particle_t* part);
int PList_next(PList* plist);
void PList_destroy(PList* plist);


#endif
