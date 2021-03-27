#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifndef PLISTS_H
#define PLISTS_H


#define MAX_SEARCH 1E5
#define LINE_MAX_LEN 256
#define NAME_MAX_LEN 96

typedef struct Part{
	double E;
	double pos[3];
	double dir[3];
} Part;

typedef int (*ReadFun)(const char* line, Part* part, double* w);

typedef struct PList{
	char pt; // Tipo de particula (n, p, ...)

	int ord; // Cantidad de archivos de tracks
	char** filenames; // Nombres de archivos de tracks (unico o uno por variable)
	FILE** files; // Archivos de tracks (unico o uno por variable)
	ReadFun* read; // Funciones de lectura

	double* trasl; // Traslacion a aplicar
	double* rot; // Rotacion a aplicar
	int x2z; // Si es true, se aplica permutacion x,y,z -> y,z,x

	char line[LINE_MAX_LEN]; // Buffer para linea de texto
	Part part; // Buffer para particula
	double w; // Buffer para peso
} PList;

PList* PList_create(char pt, int ord, const char** filenames, const ReadFun* reads,
	const double* trasl, const double* rot, int switch_x2z);
PList* PList_copy(const PList* from);
int PList_get(const PList* plist, Part* part, double* w);
int PList_next(PList* plist);
void PList_destroy(PList* plist);

PList* PListSimple_create(char pt, const char* filename, ReadFun read,
	const double* trasl, const double* rot, int switch_x2z);
PList* PListSepVar_create(char pt, const char* filenames[3], ReadFun reads[3],
	const double* trasl, const double* rot, int switch_x2z);


int PTRAC_read(const char* line, Part* part, double* w);

int T4stock_read(const char* line, Part* part, double* w);

int SSV_read(const char* line, Part* part, double* w);

int Decay_read(const char* line, Part* part, double* w);
int SSVtally_read(const char* line, Part* part, double* w);
int Isotrop_read(const char* line, Part* part, double* w);


#endif
